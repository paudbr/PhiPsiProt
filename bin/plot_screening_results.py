#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def read_rows(path):
    rows = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            ddg = row.get("ddg", "NA")
            if ddg == "NA":
                continue
            try:
                row["ddg_float"] = float(ddg)
            except ValueError:
                continue
            rows.append(row)
    return rows


parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--outdir", required=True)
args = parser.parse_args()

outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

rows = read_rows(args.input)
rows = sorted(rows, key=lambda r: r["ddg_float"])


if not rows:
    plt.figure(figsize=(8, 4))
    plt.text(0.5, 0.5, "No numeric ddG values found after filtering",
             ha="center", va="center")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(outdir / "ddg_ranking_barplot.png", dpi=300)
    plt.savefig(outdir / "top10_stabilizing_mutations.png", dpi=300)
    plt.savefig(outdir / "ddg_heatmap.png", dpi=300)
    plt.close()
    raise SystemExit(0)

mutations = [r["mutation"] for r in rows]
ddgs = [r["ddg_float"] for r in rows]

plt.figure(figsize=(max(8, len(rows) * 0.35), 6))
plt.bar(mutations, ddgs)
plt.axhline(0, linestyle="--", linewidth=1)
plt.xticks(rotation=90)
plt.ylabel("ddG mutant - WT")
plt.xlabel("Mutation")
plt.title("PyRosetta ddG ranking")
plt.tight_layout()
plt.savefig(outdir / "ddg_ranking_barplot.png", dpi=300)
plt.close()

top = rows[:10]
plt.figure(figsize=(8, 5))
plt.bar([r["mutation"] for r in top], [r["ddg_float"] for r in top])
plt.axhline(0, linestyle="--", linewidth=1)
plt.xticks(rotation=45, ha="right")
plt.ylabel("ddG mutant - WT")
plt.xlabel("Mutation")
plt.title("Top stabilizing mutations")
plt.tight_layout()
plt.savefig(outdir / "top10_stabilizing_mutations.png", dpi=300)
plt.close()

positions = sorted(set(r["position"] for r in rows), key=lambda x: int(x) if x.isdigit() else x)
mutants = list("ACDEFGHIKLMNPQRSTVWY")

matrix = []
for pos in positions:
    row_values = []
    for aa in mutants:
        vals = [r["ddg_float"] for r in rows if r["position"] == pos and r["mutant"] == aa]
        row_values.append(vals[0] if vals else float("nan"))
    matrix.append(row_values)

plt.figure(figsize=(12, max(4, len(positions) * 0.6)))
plt.imshow(matrix, aspect="auto")
plt.colorbar(label="ddG")
plt.xticks(range(len(mutants)), mutants)
plt.yticks(range(len(positions)), positions)
plt.xlabel("Mutant amino acid")
plt.ylabel("Position")
plt.title("ddG heatmap")
plt.tight_layout()
plt.savefig(outdir / "ddg_heatmap.png", dpi=300)
plt.close()
