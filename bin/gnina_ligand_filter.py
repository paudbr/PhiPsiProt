#!/usr/bin/env python3

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
import math
import matplotlib
matplotlib.use('Agg')  # sin display
import matplotlib.pyplot as plt
from scipy import stats


WATER_NAMES = {"HOH", "WAT", "DOD"}

def clean_receptor(reference_pdb, out_file):
    """Escribe un receptor solo-proteina: conserva ATOM/TER, elimina
    HETATM (ligando, aguas, iones). Necesario para que el pocket este
    libre antes de dockear."""
    kept = 0
    with open(reference_pdb) as src, open(out_file, "w") as dst:
        for line in src:
            if line.startswith("ATOM") or line.startswith("TER"):
                dst.write(line)
                if line.startswith("ATOM"):
                    kept += 1
        dst.write("END\n")
    if kept == 0:
        raise ValueError(f"No protein ATOM records found in {reference_pdb}")
    return out_file

def parse_residues(text):
    residues = []
    for token in re.split(r"[,\s]+", text.strip()):
        if not token:
            continue
        match = re.fullmatch(r"([A-Za-z]?)(-?\d+)", token)
        if not match:
            raise ValueError(f"Invalid pocket residue '{token}'. Use forms like A123 or 123.")
        chain = match.group(1) or "A"
        residues.append((chain, int(match.group(2))))
    return residues


def atom_coord(line):
    return (float(line[30:38]), float(line[38:46]), float(line[46:54]))


def box_from_residues(pdb_file, residues, size):
    coords = []
    wanted = set(residues)
    with open(pdb_file) as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            chain = line[21].strip() or "A"
            try:
                res_id = int(line[22:26])
            except ValueError:
                continue
            if (chain, res_id) in wanted:
                coords.append(atom_coord(line))
    if not coords:
        raise ValueError("No valid pocket residues found in receptor/reference PDB.")
    center = [sum(axis) / len(coords) for axis in zip(*coords)]
    return {
        "center_x": center[0],
        "center_y": center[1],
        "center_z": center[2],
        "size_x": size,
        "size_y": size,
        "size_z": size,
    }


def extract_ligand(reference_pdb, ligand_resname, out_file):
    ligand_resname = ligand_resname.upper()
    kept = 0
    with open(reference_pdb) as src, open(out_file, "w") as dst:
        for line in src:
            if not line.startswith("HETATM"):
                continue
            resname = line[17:20].strip().upper()
            if resname in WATER_NAMES:
                continue
            if ligand_resname and resname != ligand_resname:
                continue
            dst.write(line)
            kept += 1
        dst.write("END\n")
    if kept == 0:
        raise ValueError(f"No ligand atoms found for resname '{ligand_resname}' in {reference_pdb}")
    return out_file


def smiles_to_sdf(candidate_id, smiles, out_sdf):
    """Convierte un SMILES a SDF 3D con obabel (genera coords + protona a pH 7.4)."""
    cmd = [
        "obabel",
        f"-:{smiles}",
        "-osdf",
        "-O", out_sdf,
        "--gen3d",
        "-p", "7.4",
    ]
    proc = subprocess.run(cmd, text=True, capture_output=True)
    if proc.returncode != 0 or not Path(out_sdf).exists() or Path(out_sdf).stat().st_size == 0:
        raise RuntimeError(
            f"obabel falló para {candidate_id} ({smiles}):\n{proc.stdout}\n{proc.stderr}"
        )
    return out_sdf


def read_candidates(candidates_csv):
    """Lee un CSV candidate_id,smiles_ligand y genera un SDF 3D por ligando."""
    rows = []
    with open(candidates_csv, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            candidate_id = row.get("candidate_id") or row.get("id")
            smiles = row.get("smiles_ligand") or row.get("smiles")
            if not candidate_id or not smiles:
                print(f"WARNING: fila sin candidate_id o smiles, saltando: {row}")
                continue
            sdf = f"{candidate_id}_ligand.sdf"
            smiles_to_sdf(candidate_id, smiles.strip(), sdf)
            rows.append({"candidate_id": candidate_id, "ligand_sdf": sdf, "smiles": smiles.strip()})
    return rows


def parse_gnina_log(text):
    """Parsea la tabla de modos de GNINA. Devuelve scores del modo 1 (mejor pose)."""
    # Buscar líneas de la tabla: "    1       -2.85  ..."
    mode_pattern = re.compile(
        r"^\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)",
        re.MULTILINE
    )
    matches = mode_pattern.findall(text)
    if not matches:
        return {"cnn_score": None, "cnn_affinity": None, "affinity": None}

    # Modo 1 = mejor pose
    best = matches[0]
    return {
        "affinity":     float(best[1]),   # kcal/mol
        "cnn_score":    float(best[3]),   # CNN pose score
        "cnn_affinity": float(best[4]),   # CNN affinity
    }


def mean_std(xs):
    if not xs:
        return "", ""
    m = sum(xs) / len(xs)
    if len(xs) > 1:
        s = math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))
    else:
        s = 0.0
    return round(m, 4), round(s, 4)


def run_gnina(candidate_id, receptor, ligand, box, replicate, args, autobox_ref=None):
    out_sdf = f"{candidate_id}_rep{replicate}.sdf"
    log_file = f"{candidate_id}_rep{replicate}.gnina.log"

    cmd = ["gnina", "--receptor", receptor, "--ligand", ligand]

    if autobox_ref:
        # caja automatica alrededor del ligando de referencia + margen.
        # --autobox_add = buffer en Angstroms (default +4 en los seis lados).
        # --autobox_extend = booleano: expande la caja si el ligando no cabe
        #   rotando en su conformacion de entrada (clave con ligandos grandes).
        cmd += [
            "--autobox_ligand", autobox_ref,
            "--autobox_add", str(args.autobox_add),
            "--autobox_extend", "true",
        ]
    else:
        # fallback: caja explicita por centroide de residuos
        cmd += [
            "--center_x", str(box["center_x"]),
            "--center_y", str(box["center_y"]),
            "--center_z", str(box["center_z"]),
            "--size_x", str(box["size_x"]),
            "--size_y", str(box["size_y"]),
            "--size_z", str(box["size_z"]),
        ]

    cmd += [
        "--num_modes", str(args.gnina_num_modes),
        "--exhaustiveness", str(args.gnina_exhaustiveness),
        "--seed", str(args.gnina_seed + replicate),
        "--out", out_sdf,
    ]

    proc = subprocess.run(cmd, text=True, capture_output=True)
    log_text = proc.stdout + "\n" + proc.stderr
    with open(log_file, "w") as handle:
        handle.write(log_text)

    if proc.returncode != 0:
        print(f"WARNING: gnina returncode={proc.returncode} para {candidate_id} "
              f"rep{replicate}. Ver {log_file}", file=sys.stderr)

    parsed = parse_gnina_log(log_text)
    parsed.update({
        "candidate_id": candidate_id,
        "replicate": replicate,
        "returncode": proc.returncode,
        "docked_sdf": out_sdf,
        "log": log_file,
    })
    return parsed


def plot_scores(records, output_pdf="gnina_score_distributions.pdf"):
    """Pinta distribución de scores por candidato vs original y hace Mann-Whitney U."""

    metrics = {
        "affinity":     "Affinity (kcal/mol) — lower is better",
        "cnn_score":    "CNN Pose Score — higher is better",
        "cnn_affinity": "CNN Affinity — higher is better",
    }
    grouped = {}
    for rec in records:
        grouped.setdefault(rec["candidate_id"], []).append(rec)

    candidates = [c for c in grouped if c != "original"]
    from matplotlib.backends.backend_pdf import PdfPages

    if not candidates:
        print("WARNING: no hay candidatos para plotear, generando PDF vacío")
        with PdfPages(output_pdf) as pdf:
            fig, ax = plt.subplots()
            ax.axis("off")
            ax.text(0.5, 0.5, "No candidate ligands docked", ha="center")
            pdf.savefig(fig); plt.close(fig)
        return

    original_vals = {
        m: [r[m] for r in grouped.get("original", []) if r[m] is not None]
        for m in metrics
    }

    with PdfPages(output_pdf) as pdf:
        for metric, ylabel in metrics.items():
            fig, ax = plt.subplots(figsize=(max(8, len(candidates) * 2), 5))

            orig = original_vals[metric]
            orig_mean = sum(orig) / len(orig) if orig else None

            positions = []
            labels    = []
            stat_lines = []

            for i, cid in enumerate(candidates):
                vals = [r[metric] for r in grouped[cid] if r[metric] is not None]
                if not vals:
                    continue

                x = i + 1
                ax.scatter([x] * len(vals), vals, color="steelblue", zorder=3, s=60)
                m = sum(vals) / len(vals)
                ax.hlines(m, x - 0.25, x + 0.25, colors="steelblue", linewidth=2.5)

                if orig and len(vals) >= 1:
                    if len(vals) >= 2 and len(orig) >= 2:
                        stat, pval = stats.mannwhitneyu(vals, orig, alternative="two-sided")
                        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
                        stat_lines.append(f"{cid}: U={stat:.0f}, p={pval:.3f} {sig}")
                    else:
                        stat_lines.append(f"{cid}: n insuficiente para test")

                positions.append(x)
                labels.append(cid)

            if orig_mean is not None:
                ax.axhline(orig_mean, color="tomato", linestyle="--", linewidth=1.5,
                           label=f"Original mean ({orig_mean:.3f})")
                if len(orig) > 1:
                    orig_std = math.sqrt(sum((v - orig_mean)**2 for v in orig) / (len(orig)-1))
                    ax.axhspan(orig_mean - orig_std, orig_mean + orig_std,
                               alpha=0.15, color="tomato", label="Original ±1 SD")

            ax.set_xticks(positions)
            ax.set_xticklabels(labels, rotation=20, ha="right")
            ax.set_ylabel(ylabel)
            ax.set_title(f"{metric} — binders vs original")
            ax.legend(fontsize=8)

            if stat_lines:
                stats_text = "\n".join(stat_lines)
                fig.text(0.01, -0.05 * len(stat_lines),
                         "Mann-Whitney U (two-sided, vs original):\n" + stats_text,
                         fontsize=7, va="top",
                         bbox=dict(boxstyle="round", fc="lightyellow", ec="gray"))

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        # Página resumen con tabla
        fig, ax = plt.subplots(figsize=(10, max(3, len(candidates) * 0.5 + 2)))
        ax.axis("off")
        table_data = [["Candidate", "Affinity mean±SD", "CNN Score mean±SD",
                        "CNN Aff mean±SD", "p (affinity)", "p (cnn_aff)"]]

        def fmt(xs):
            if not xs:
                return "n/a"
            m = sum(xs)/len(xs)
            s = math.sqrt(sum((x-m)**2 for x in xs)/(len(xs)-1)) if len(xs) > 1 else 0
            return f"{m:.3f}±{s:.3f}"

        def pval_str(a, b):
            if len(a) >= 2 and len(b) >= 2:
                _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
                sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
                return f"{p:.3f} {sig}"
            return "n/a"

        for cid in candidates:
            vals_aff = [r["affinity"]     for r in grouped[cid] if r["affinity"] is not None]
            vals_cnn = [r["cnn_score"]    for r in grouped[cid] if r["cnn_score"] is not None]
            vals_ca  = [r["cnn_affinity"] for r in grouped[cid] if r["cnn_affinity"] is not None]
            table_data.append([
                cid,
                fmt(vals_aff),
                fmt(vals_cnn),
                fmt(vals_ca),
                pval_str(vals_aff, original_vals["affinity"]),
                pval_str(vals_ca,  original_vals["cnn_affinity"]),
            ])

        tbl = ax.table(cellText=table_data[1:], colLabels=table_data[0],
                       loc="center", cellLoc="center")
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(8)
        tbl.auto_set_column_width(col=list(range(len(table_data[0]))))
        ax.set_title("Summary table — Mann-Whitney U vs original", pad=12)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

    print(f"Saved {output_pdf}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--reference-pdb", required=True)
    parser.add_argument("--ligand-file")
    parser.add_argument("--ligand-resname", default="GOL")
    parser.add_argument("--pocket-residues", required=True)
    parser.add_argument("--box-size", type=float, default=20.0)
    parser.add_argument("--runs", type=int, default=3)
    parser.add_argument("--gnina-exhaustiveness", type=int, default=8)
    parser.add_argument("--gnina-num-modes", type=int, default=9)
    parser.add_argument("--gnina-seed", type=int, default=1000)
    parser.add_argument("--autobox-ligand", default=None,
                        help="PDB/SDF de referencia para autobox. Si se da, "
                             "la caja se construye alrededor de este ligando.")
    parser.add_argument("--autobox-add", type=float, default=4.0,
                        help="Buffer en Angstrom anadido a la autobox en los "
                             "seis lados (GNINA --autobox_add, default 4).")
    parser.add_argument("--no-autobox", action="store_true",
                        help="Desactiva autobox y fuerza la caja por residuos.")
    args = parser.parse_args()

    residues = parse_residues(args.pocket_residues)
    box = box_from_residues(args.reference_pdb, residues, args.box_size)


    receptor = "receptor_clean.pdb"
    clean_receptor(args.reference_pdb, receptor)

    # Cada candidato es un LIGANDO (SMILES -> SDF 3D)
    candidates = read_candidates(args.candidates)

    # "original": el ligando co-cristalizado extraido del reference (MK1/saquinavir)
    original_ligand = extract_ligand(args.reference_pdb, args.ligand_resname, "extracted_ligand.pdb")

    # ── Decidir autobox ───────────────────────────────────────────
    # Por defecto: autobox alrededor del ligando co-cristalizado extraido,
    # que garantiza que la caja contiene el sitio activo con margen
    # (resuelve el "ligand out of the box" con ligandos grandes/flexibles).
    autobox_ref = None
    if not args.no_autobox:
        autobox_ref = args.autobox_ligand or original_ligand
        if autobox_ref and not Path(autobox_ref).exists():
            print(f"WARNING: autobox ref '{autobox_ref}' no existe; "
                  f"usando caja por residuos.", file=sys.stderr)
            autobox_ref = None
    if autobox_ref:
        print(f"[gnina] usando autobox sobre '{autobox_ref}' "
              f"(add={args.autobox_add} A, extend=on)", file=sys.stderr)
    else:
        print(f"[gnina] usando caja por residuos centrada en "
              f"({box['center_x']:.1f}, {box['center_y']:.1f}, {box['center_z']:.1f}) "
              f"lado={args.box_size} A", file=sys.stderr)

    records = []
    # Ligandos nuevos
    for candidate in candidates:
        for rep in range(args.runs):
            records.append(
                run_gnina(candidate["candidate_id"], receptor,
                          candidate["ligand_sdf"], box, rep, args,
                          autobox_ref=autobox_ref)
            )
    # Ligando original como referencia de comparacion
    for rep in range(args.runs):
        records.append(
            run_gnina("original", receptor, original_ligand, box, rep, args,
                      autobox_ref=autobox_ref)
        )

    with open("gnina_scores.tsv", "w", newline="") as handle:
        fields = ["candidate_id", "replicate", "cnn_score", "cnn_affinity",
                  "affinity", "returncode", "docked_sdf", "log"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)

    grouped = {}
    for rec in records:
        grouped.setdefault(rec["candidate_id"], []).append(rec)

    with open("gnina_summary.tsv", "w", newline="") as handle:
        fields = [
            "candidate_id", "n",
            "cnn_score_mean", "cnn_score_std",
            "cnn_affinity_mean", "cnn_affinity_std",
            "affinity_mean", "affinity_std"
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for candidate_id, vals in grouped.items():
            row = {"candidate_id": candidate_id, "n": len(vals)}
            for metric in ["cnn_score", "cnn_affinity", "affinity"]:
                xs = [v[metric] for v in vals if v[metric] is not None]
                m, s = mean_std(xs)
                row[f"{metric}_mean"] = m
                row[f"{metric}_std"]  = s
            writer.writerow(row)

    plot_scores(records, output_pdf="gnina_score_distributions.pdf")


if __name__ == "__main__":
    main()
