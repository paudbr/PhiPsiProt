#!/usr/bin/env python3
import os, glob, re
import numpy as np

TAB = chr(9)
NL  = chr(10)

# ── sample id desde variable de entorno ──────────────────────────────────────
sample = os.environ.get("SAMPLE_ID", "")
if not sample:
    raise RuntimeError("SAMPLE_ID env var not set")

print(f"sample={sample}, cwd={os.getcwd()}")
print(f"raw/ contents: {os.listdir('raw') if os.path.isdir('raw') else 'NO raw/ DIR'}")

# ── mapeo ranked_N → path cif ────────────────────────────────────────────────
ranked_cifs = glob.glob("raw/model_idx_*_ranked_*.cif")
print(f"ranked_cifs found: {ranked_cifs}")

ranked_map = {}
for p in ranked_cifs:
    m = re.search(r'ranked_(\d+)\.cif$', p)
    if m:
        ranked_map[int(m.group(1))] = p

# ── parsear CIF — columna B_iso_or_equiv para CA atoms ───────────────────────
def parse_cif_plddt(path):
    cols    = []
    col_idx = {}
    rows    = []
    in_atom = False

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line == 'loop_':
                cols    = []
                col_idx = {}
                in_atom = False
                continue
            if line.startswith('_atom_site.'):
                cols.append(line)
                in_atom = True
                continue
            if in_atom and cols and line.startswith(('ATOM', 'HETATM')):
                if not col_idx:
                    col_idx = {c: i for i, c in enumerate(cols)}
                parts = line.split()
                if len(parts) < len(cols):
                    continue
                try:
                    atom = parts[col_idx['_atom_site.label_atom_id']]
                except KeyError:
                    continue
                if atom != 'CA':
                    continue
                try:
                    chain   = parts[col_idx['_atom_site.auth_asym_id']]
                    res_num = int(parts[col_idx['_atom_site.auth_seq_id']])
                    bfac    = float(parts[col_idx['_atom_site.B_iso_or_equiv']])
                    rows.append((chain, res_num, bfac))
                except (KeyError, ValueError):
                    continue
            elif in_atom and cols and line.startswith('#'):
                in_atom = False
                col_idx = {}
                cols    = []
    return rows

# ── pLDDT TSV ─────────────────────────────────────────────────────────────────
plddt_rows = []
for rank_n in sorted(ranked_map):
    residues = parse_cif_plddt(ranked_map[rank_n])
    print(f"  ranked_{rank_n}: {len(residues)} CA atoms")
    for chain, res, plddt in residues:
        plddt_rows.append((f"ranked_{rank_n}", chain, res, plddt))

with open(f"{sample}_chai1_allmodels_plddt.tsv", "w") as fh:
    fh.write("model" + TAB + "chain" + TAB + "residue" + TAB + "plddt" + NL)
    for model, chain, res, plddt in plddt_rows:
        fh.write(model + TAB + chain + TAB + str(res) + TAB + f"{plddt:.4f}" + NL)

# ── scores NPZ ───────────────────────────────────────────────────────────────
idx_to_rank = {}
for rank_n, path in ranked_map.items():
    m = re.search(r'model_idx_(\d+)_ranked', path)
    if m:
        idx_to_rank[int(m.group(1))] = rank_n

score_rows = []
for npz_path in sorted(glob.glob("raw/scores.model_idx_*.npz")):
    m = re.search(r'model_idx_(\d+)\.npz$', npz_path)
    if not m:
        continue
    idx   = int(m.group(1))
    rank_n = idx_to_rank.get(idx)
    if rank_n is None:
        continue
    d = np.load(npz_path)
    score_rows.append({
        "model":           f"ranked_{rank_n}",
        "ptm":             float(d["ptm"].flat[0])             if "ptm"             in d else None,
        "iptm":            float(d["iptm"].flat[0])            if "iptm"            in d else None,
        "aggregate_score": float(d["aggregate_score"].flat[0]) if "aggregate_score" in d else None,
    })

score_rows.sort(key=lambda r: int(re.search(r'\d+$', r["model"]).group()))

with open(f"{sample}_chai1_scores.tsv", "w") as fh:
    fh.write("model" + TAB + "ptm" + TAB + "iptm" + TAB + "aggregate_score" + NL)
    for r in score_rows:
        ptm  = f"{r['ptm']:.4f}"             if r["ptm"]             is not None else "NA"
        iptm = f"{r['iptm']:.4f}"            if r["iptm"]            is not None else "NA"
        agg  = f"{r['aggregate_score']:.4f}" if r["aggregate_score"] is not None else "NA"
        fh.write(r["model"] + TAB + ptm + TAB + iptm + TAB + agg + NL)

print(f"pLDDT rows: {len(plddt_rows)}, score rows: {len(score_rows)}")