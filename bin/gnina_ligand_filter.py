#!/usr/bin/env python3

import argparse
import csv
import os
import re
import shutil
import subprocess
from pathlib import Path
import math
import matplotlib
matplotlib.use('Agg')  # sin display
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats


WATER_NAMES = {"HOH", "WAT", "DOD"}


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


def read_candidates(candidates_csv):
    base_dir = Path(candidates_csv).resolve().parent
    rows = []
    with open(candidates_csv, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            candidate_id = row.get("candidate_id") or row.get("id") or Path(row.get("pdb", "")).stem
            pdb = row.get("pdb")
            if not pdb:
                continue
            pdb_path = Path(pdb)
            if not pdb_path.is_absolute():
                pdb_path = base_dir / pdb_path

            # ← NUEVO: copiar al workdir si no está ya aquí
            local_pdb = Path(pdb_path.name)
            if not local_pdb.exists():
                if pdb_path.exists():
                    shutil.copyfile(str(pdb_path), str(local_pdb))
                else:
                    print(f"WARNING: PDB not found: {pdb_path}, skipping")
                    continue

            rows.append({"candidate_id": candidate_id, "pdb": str(local_pdb), **row})
    return rows


import math

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



def run_gnina(candidate_id, receptor, ligand, box, replicate, args):
    out_sdf = f"{candidate_id}_rep{replicate}.sdf"
    log_file = f"{candidate_id}_rep{replicate}.gnina.log"
    cmd = [
        "gnina",
        "--receptor", receptor,
        "--ligand", ligand,
        "--center_x", str(box["center_x"]),
        "--center_y", str(box["center_y"]),
        "--center_z", str(box["center_z"]),
        "--size_x", str(box["size_x"]),
        "--size_y", str(box["size_y"]),
        "--size_z", str(box["size_z"]),
        "--num_modes", str(args.gnina_num_modes),
        "--exhaustiveness", str(args.gnina_exhaustiveness),
        "--seed", str(args.gnina_seed + replicate),
        "--out", out_sdf,
    ]
    proc = subprocess.run(cmd, text=True, capture_output=True)
    log_text = proc.stdout + "\n" + proc.stderr
    with open(log_file, "w") as handle:
        handle.write(log_text)
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
    # Separar original del resto
    grouped = {}
    for rec in records:
        grouped.setdefault(rec["candidate_id"], []).append(rec)

    original_vals = {
        m: [r[m] for r in grouped.get("original", []) if r[m] is not None]
        for m in metrics
    }
    candidates = [c for c in grouped if c != "original"]

    from matplotlib.backends.backend_pdf import PdfPages

    with PdfPages(output_pdf) as pdf:
        for metric, ylabel in metrics.items():
            fig, ax = plt.subplots(figsize=(max(8, len(candidates) * 2), 5))

            orig = original_vals[metric]
            orig_mean = sum(orig) / len(orig) if orig else None

            positions = []
            labels    = []
            colors    = []
            stat_lines = []

            for i, cid in enumerate(candidates):
                vals = [r[metric] for r in grouped[cid] if r[metric] is not None]
                if not vals:
                    continue

                x = i + 1
                # Puntos individuales
                ax.scatter([x] * len(vals), vals, color="steelblue", zorder=3, s=60)
                # Media
                m = sum(vals) / len(vals)
                ax.hlines(m, x - 0.25, x + 0.25, colors="steelblue", linewidth=2.5)

                # Mann-Whitney U vs original
                if orig and len(vals) >= 1:
                    if len(vals) >= 2 and len(orig) >= 2:
                        stat, pval = stats.mannwhitneyu(vals, orig, alternative="two-sided")
                        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
                        stat_lines.append(f"{cid}: U={stat:.0f}, p={pval:.3f} {sig}")
                    else:
                        stat_lines.append(f"{cid}: n insuficiente para test")

                positions.append(x)
                labels.append(cid)
                colors.append("steelblue")

            # Línea de la media del original
            if orig_mean is not None:
                ax.axhline(orig_mean, color="tomato", linestyle="--", linewidth=1.5,
                           label=f"Original mean ({orig_mean:.3f})")
                # Banda ±std del original
                if len(orig) > 1:
                    orig_std = math.sqrt(sum((v - orig_mean)**2 for v in orig) / (len(orig)-1))
                    ax.axhspan(orig_mean - orig_std, orig_mean + orig_std,
                               alpha=0.15, color="tomato", label="Original ±1 SD")

            ax.set_xticks(positions)
            ax.set_xticklabels(labels, rotation=20, ha="right")
            ax.set_ylabel(ylabel)
            ax.set_title(f"{metric} — binders vs original")
            ax.legend(fontsize=8)

            # Añadir stats como texto en la figura
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
        for cid in candidates:
            vals_aff = [r["affinity"]     for r in grouped[cid] if r["affinity"] is not None]
            vals_cnn = [r["cnn_score"]    for r in grouped[cid] if r["cnn_score"] is not None]
            vals_ca  = [r["cnn_affinity"] for r in grouped[cid] if r["cnn_affinity"] is not None]

            def fmt(xs):
                if not xs: return "n/a"
                m = sum(xs)/len(xs)
                s = math.sqrt(sum((x-m)**2 for x in xs)/(len(xs)-1)) if len(xs)>1 else 0
                return f"{m:.3f}±{s:.3f}"

            def pval_str(a, b):
                if len(a) >= 2 and len(b) >= 2:
                    _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
                    sig = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns"
                    return f"{p:.3f} {sig}"
                return "n/a"

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
    args = parser.parse_args()

    residues = parse_residues(args.pocket_residues)
    box = box_from_residues(args.reference_pdb, residues, args.box_size)

    ligand_file = args.ligand_file
    if ligand_file:
        staged_ligand = Path(ligand_file).name
        if Path(ligand_file).resolve() != Path(staged_ligand).resolve():
            shutil.copyfile(ligand_file, staged_ligand)
        ligand_file = staged_ligand
    else:
        ligand_file = extract_ligand(args.reference_pdb, args.ligand_resname, "extracted_ligand.pdb")

    candidates = read_candidates(args.candidates)
    candidates.append({"candidate_id": "original", "pdb": args.reference_pdb})

    records = []
    for candidate in candidates:
        for rep in range(args.runs):
            records.append(run_gnina(candidate["candidate_id"], candidate["pdb"], ligand_file, box, rep, args))

    with open("gnina_scores.tsv", "w", newline="") as handle:
        fields = ["candidate_id", "replicate", "cnn_score", "cnn_affinity", "affinity", "returncode", "docked_sdf", "log"]
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
