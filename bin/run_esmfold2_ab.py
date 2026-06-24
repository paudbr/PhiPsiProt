#!/usr/bin/env python3.12
"""
run_esmfold2_ab.py  —  batch mode
──────────────────────────────────
Carga ESMFold2-Fast UNA SOLA VEZ y procesa todos los diseños en un loop.
Escribe un TSV con métricas por diseño + un CIF por diseño.

Usage:
    python3.12 run_esmfold2_ab.py \
        --designs-dir  extracted_pdbs/ \
        --target       antigen_truncated.pdb \
        --framework    h-NbBCII10.pdb \
        --antigen-chain B \
        --ab-type      VHH \
        --out-tsv      esm2_all_metrics.tsv \
        --out-cif-dir  esm2_cifs/

Metrics per design (TSV columns):
    id, iptm_ht, iptm_global, pae_global, pae_interface,
    plddt_global, plddt_cdr, rmsd_global, rmsd_cdr
"""

import argparse
import sys
import os
import glob
import numpy as np
import torch


# ── AA lookup ─────────────────────────────────────────────────────────────────
AA3TO1 = dict(
    ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q",
    GLU="E", GLY="G", HIS="H", ILE="I", LEU="L", LYS="K",
    MET="M", PHE="F", PRO="P", SER="S", THR="T", TRP="W",
    TYR="Y", VAL="V"
)


# ── PDB helpers ───────────────────────────────────────────────────────────────
def read_ca(pdb_path, chain_id=None):
    ca = {}
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if line[13:15].strip() != "CA":
                continue
            ch  = line[21]
            if chain_id and ch != chain_id:
                continue
            rno = int(line[22:26])
            aa  = AA3TO1.get(line[17:20].strip(), "X")
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            ca[(ch, rno)] = (aa, np.array([x, y, z]))
    return ca


def extract_seq(pdb_path, chain_id):
    ca = read_ca(pdb_path, chain_id)
    return "".join(v[0] for k, v in sorted(ca.items()))


def extract_coords(pdb_path, chain_id):
    ca = read_ca(pdb_path, chain_id)
    if not ca:
        return np.array([]).reshape(0, 3)
    return np.array([v[1] for k, v in sorted(ca.items())])


def read_cdr_resnums(framework_pdb):
    resnums = []
    with open(framework_pdb) as fh:
        for line in fh:
            if "PDBinfo-LABEL" in line:
                parts = line.split()
                try:
                    resnums.append(int(parts[-2]))
                except (ValueError, IndexError):
                    pass
    return set(resnums)


# ── Geometry helpers ──────────────────────────────────────────────────────────
def superpose_rmsd(coords1, coords2):
    """Kabsch RMSD. Returns 999.0 on empty or mismatched input."""
    if len(coords1) == 0 or len(coords1) != len(coords2):
        return 999.0
    c1 = coords1 - coords1.mean(axis=0)
    c2 = coords2 - coords2.mean(axis=0)
    H  = c2.T @ c1
    U, S, Vt = np.linalg.svd(H)
    d  = np.linalg.det(Vt.T @ U.T)
    D  = np.diag([1.0, 1.0, d])
    R  = Vt.T @ D @ U.T
    c2r = c2 @ R.T
    return float(np.sqrt(((c1 - c2r) ** 2).sum(axis=1).mean()))


def parse_ca_from_cif(cif_str, chain_id):
    """Extract Cα coords from mmCIF string for a given label_asym_id chain."""
    coords = []
    for line in cif_str.splitlines():
        if not line.startswith("ATOM"):
            continue
        parts = line.split()
        if len(parts) < 19:
            continue
        if parts[3] != "CA" or parts[6] != chain_id:
            continue
        try:
            coords.append(np.array([float(parts[15]),
                                     float(parts[16]),
                                     float(parts[17])]))
        except (ValueError, IndexError):
            continue
    return np.array(coords) if coords else np.array([]).reshape(0, 3)


# ── Model loading ─────────────────────────────────────────────────────────────
def load_model(hf_home):
    from transformers.models.esmfold2.modeling_esmfold2 import (
        ESMFold2Model, ESMFold2Config
    )

    # Resolve ESMFold2-Fast snapshot
    snap_base  = os.path.join(hf_home, "models--biohub--ESMFold2-Fast", "snapshots")
    snaps      = sorted(glob.glob(os.path.join(snap_base, "*")))
    if not snaps:
        print(f"ERROR: no ESMFold2-Fast snapshots in {snap_base}", file=sys.stderr)
        sys.exit(1)
    snap_path  = snaps[-1]
    print(f"[ESMFold2] ESMFold2-Fast snapshot: {snap_path}", file=sys.stderr)

    # Load config, patch esmc_id to local ESMC-6B snapshot
    config     = ESMFold2Config.from_pretrained(snap_path, local_files_only=True)
    esmc_name  = config.esmc_id.split("/")[-1]          # e.g. "ESMC-6B"
    esmc_base  = os.path.join(hf_home, f"models--biohub--{esmc_name}", "snapshots")
    esmc_snaps = sorted(glob.glob(os.path.join(esmc_base, "*")))
    if not esmc_snaps:
        print(f"ERROR: no {esmc_name} snapshots in {esmc_base}", file=sys.stderr)
        sys.exit(1)
    esmc_snap      = esmc_snaps[-1]
    config.esmc_id = esmc_snap
    print(f"[ESMFold2] ESMC backbone: {esmc_snap}", file=sys.stderr)

    print("[ESMFold2] Loading model weights...", file=sys.stderr)
    model = ESMFold2Model.from_pretrained(
        snap_path,
        config=config,
        local_files_only=True
    ).cuda().eval()
    print("[ESMFold2] Model loaded OK", file=sys.stderr)
    return model


def load_builder(hf_home):
    from esm.models.esmfold2 import ESMFold2InputBuilder

    # Resolve ccd.pkl directory
    ccd_snaps = sorted(glob.glob(
        os.path.join(hf_home, "models--biohub--ESMFold2", "snapshots", "*", "ccd.pkl")
    ))
    ccd_dir = os.path.dirname(ccd_snaps[-1]) if ccd_snaps else None
    print(f"[ESMFold2] CCD dir: {ccd_dir}", file=sys.stderr)

    builder = ESMFold2InputBuilder(ccd_cache=ccd_dir) if ccd_dir else ESMFold2InputBuilder()
    return builder


# ── Per-design inference ──────────────────────────────────────────────────────
def run_one(model, builder, design_pdb, antigen_seq, antigen_chain,
            ab_type, cdr_set, num_loops, num_steps, seed):
    from esm.models.esmfold2 import ProteinInput, StructurePredictionInput

    sample_id   = os.path.splitext(os.path.basename(design_pdb))[0]
    heavy_seq   = extract_seq(design_pdb, "H")
    if not heavy_seq:
        print(f"  WARNING: no chain H in {design_pdb}, skipping", file=sys.stderr)
        return None

    L_H = len(heavy_seq)
    sequences = [ProteinInput(id="H", sequence=heavy_seq)]
    if ab_type == "scFv":
        light_seq = extract_seq(design_pdb, "L")
        if light_seq:
            sequences.append(ProteinInput(id="L", sequence=light_seq))
            L_H += len(light_seq)
    sequences.append(ProteinInput(id="T", sequence=antigen_seq))

    spi = StructurePredictionInput(sequences=sequences)

    with torch.no_grad():
        result = builder.fold(
            model, spi,
            num_loops=num_loops,
            num_sampling_steps=num_steps,
            num_diffusion_samples=1,
            seed=seed
        )

    # ── Metrics ──────────────────────────────────────────────────────────
    iptm_global = float(result.iptm)
    pm          = result.pair_chains_iptm.cpu().numpy()
    iptm_ht     = float(pm[0, -1])   # H↔T

    plddt        = result.plddt.cpu().numpy()
    plddt_global = float(plddt.mean())
    plddt_h      = plddt[:L_H]

    # Map CDR resnums → 0-based indices
    ca_design    = read_ca(design_pdb, "H")
    resno_list   = sorted(k[1] for k in ca_design.keys())
    resno_to_idx = {rno: i for i, rno in enumerate(resno_list)}
    cdr_indices  = [resno_to_idx[r] for r in sorted(cdr_set) if r in resno_to_idx]

    if cdr_indices:
        plddt_cdr = float(plddt_h[np.array(cdr_indices)].mean())
    else:
        plddt_cdr = plddt_global

    pae           = result.pae.cpu().numpy()
    pae_global    = float(pae.mean())
    pae_interface = float(pae[:L_H, L_H:].mean())

    # ── CIF + RMSD ───────────────────────────────────────────────────────
    cif_str         = result.complex.to_mmcif()
    pred_coords_H   = parse_ca_from_cif(cif_str, "H")
    design_coords_H = extract_coords(design_pdb, "H")
    min_len         = min(len(design_coords_H), len(pred_coords_H))

    rmsd_global = superpose_rmsd(
        design_coords_H[:min_len], pred_coords_H[:min_len]
    )
    if cdr_indices:
        valid_cdr = [i for i in cdr_indices if i < min_len]
        rmsd_cdr  = superpose_rmsd(
            design_coords_H[valid_cdr], pred_coords_H[valid_cdr]
        ) if valid_cdr else 999.0
    else:
        rmsd_cdr = 999.0

    print(f"  {sample_id}: iptm_ht={iptm_ht:.4f} pae_int={pae_interface:.2f} "
          f"plddt_cdr={plddt_cdr:.4f} rmsd_cdr={rmsd_cdr:.3f}", file=sys.stderr)

    return {
        "id":            sample_id,
        "iptm_ht":       iptm_ht,
        "iptm_global":   iptm_global,
        "pae_global":    pae_global,
        "pae_interface": pae_interface,
        "plddt_global":  plddt_global,
        "plddt_cdr":     plddt_cdr,
        "rmsd_global":   rmsd_global,
        "rmsd_cdr":      rmsd_cdr,
        "cif_str":       cif_str,
    }


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="ESMFold2 antibody batch metrics")
    parser.add_argument("--designs-dir",  required=True,
                        help="Directory containing design PDB files (chain H[L])")
    parser.add_argument("--target",       required=True,
                        help="Antigen PDB (truncated)")
    parser.add_argument("--framework",    required=True,
                        help="Framework HLT PDB (CDR REMARKs)")
    parser.add_argument("--antigen-chain", default="B",
                        help="Chain ID in target PDB (default: B)")
    parser.add_argument("--ab-type",      default="VHH", choices=["VHH", "scFv"])
    parser.add_argument("--out-tsv",      required=True,
                        help="Output TSV with all metrics")
    parser.add_argument("--out-cif-dir",  required=True,
                        help="Output directory for per-design CIF files")
    parser.add_argument("--num-loops",    type=int, default=1,
                        help="ESMFold2 diffusion loops (default 1 for speed)")
    parser.add_argument("--num-steps",    type=int, default=10,
                        help="ESMFold2 sampling steps (default 10 for speed)")
    parser.add_argument("--seed",         type=int, default=42)
    args = parser.parse_args()

    os.makedirs(args.out_cif_dir, exist_ok=True)

    # ── HuggingFace cache dir ─────────────────────────────────────────────
    hf_home = os.environ.get("HF_HOME", "/mnt/alphafold3_db/ESMFold2/hub")
    print(f"[ESMFold2] HF_HOME: {hf_home}", file=sys.stderr)

    # ── Read CDR positions ────────────────────────────────────────────────
    cdr_set = read_cdr_resnums(args.framework)
    print(f"[ESMFold2] CDR residues: {sorted(cdr_set)}", file=sys.stderr)

    # ── Read antigen sequence ─────────────────────────────────────────────
    antigen_seq = extract_seq(args.target, args.antigen_chain)
    if not antigen_seq:
        print(f"ERROR: no chain {args.antigen_chain} in {args.target}", file=sys.stderr)
        sys.exit(1)
    print(f"[ESMFold2] Antigen: {len(antigen_seq)}aa (chain {args.antigen_chain})",
          file=sys.stderr)

    # ── Collect design PDBs ───────────────────────────────────────────────
    pdbs = sorted(glob.glob(os.path.join(args.designs_dir, "*.pdb")))
    if not pdbs:
        print(f"ERROR: no PDB files in {args.designs_dir}", file=sys.stderr)
        sys.exit(1)
    print(f"[ESMFold2] {len(pdbs)} designs to process", file=sys.stderr)

    # ── Load model ONCE ───────────────────────────────────────────────────
    model   = load_model(hf_home)
    builder = load_builder(hf_home)

    # ── Process all designs ───────────────────────────────────────────────
    header = ("id\tiptm_ht\tiptm_global\tpae_global\tpae_interface\t"
              "plddt_global\tplddt_cdr\trmsd_global\trmsd_cdr\n")

    with open(args.out_tsv, "w") as tsv_fh:
        tsv_fh.write(header)

        for i, pdb_path in enumerate(pdbs):
            sample_id = os.path.splitext(os.path.basename(pdb_path))[0]
            print(f"[ESMFold2] [{i+1}/{len(pdbs)}] {sample_id}", file=sys.stderr)

            try:
                metrics = run_one(
                    model, builder,
                    design_pdb    = pdb_path,
                    antigen_seq   = antigen_seq,
                    antigen_chain = args.antigen_chain,
                    ab_type       = args.ab_type,
                    cdr_set       = cdr_set,
                    num_loops     = args.num_loops,
                    num_steps     = args.num_steps,
                    seed          = args.seed,
                )
            except torch.cuda.OutOfMemoryError:
                print(f"  WARNING: OOM for {sample_id}, skipping", file=sys.stderr)
                torch.cuda.empty_cache()
                continue
            except Exception as e:
                print(f"  WARNING: error for {sample_id}: {e}, skipping", file=sys.stderr)
                continue

            if metrics is None:
                continue

            # Write CIF
            cif_path = os.path.join(args.out_cif_dir, f"{sample_id}_esm2.cif")
            with open(cif_path, "w") as fh:
                fh.write(metrics["cif_str"])

            # Write TSV row
            tsv_fh.write(
                f"{metrics['id']}\t"
                f"{metrics['iptm_ht']:.4f}\t"
                f"{metrics['iptm_global']:.4f}\t"
                f"{metrics['pae_global']:.3f}\t"
                f"{metrics['pae_interface']:.3f}\t"
                f"{metrics['plddt_global']:.4f}\t"
                f"{metrics['plddt_cdr']:.4f}\t"
                f"{metrics['rmsd_global']:.3f}\t"
                f"{metrics['rmsd_cdr']:.3f}\n"
            )
            tsv_fh.flush()   # escribe línea a línea por si hay crash

    print(f"[ESMFold2] Done. Written {args.out_tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()