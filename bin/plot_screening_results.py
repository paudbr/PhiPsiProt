#!/usr/bin/env python3

"""
Generate plots for PhiPsiProt Mode 1 screening results.

Purpose
-------
Create visual summaries from ranked screening candidates.

Generated plots
---------------
ddg_ranking_barplot.png:
    Barplot of all candidate ddG values ordered by final_score.

top10_stabilizing_mutations.png:
    Horizontal barplot of the top 10 most stabilizing mutations.

ddg_heatmap.png:
    Heatmap of ddG values by residue position and mutant amino acid.

Notes
-----
- The plots use PyRosetta ddG values.
- Some ddG values are clipped only for visualization readability.
- The underlying CSV values are not modified.
"""

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


def to_float(value):
    try:
        if value in ["", "NA", None]:
            return None
        return float(value)
    except Exception:
        return None


def read_rows(path):
    rows = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            ddg = to_float(row.get("ddg"))
            final_score = to_float(row.get("final_score"))

            if ddg is None:
                continue

            row["ddg_float"] = ddg
            row["final_score_float"] = final_score if final_score is not None else ddg
            rows.append(row)

    return rows


def clip_ddg(value, lower=-25, upper=25):
    return max(lower, min(upper, value))


parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--outdir", required=True)
args = parser.parse_args()

outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

rows = read_rows(args.input)
rows = sorted(rows, key=lambda r: r["final_score_float"])

if not rows:
    for name in [
        "ddg_ranking_barplot.png",
        "top10_stabilizing_mutations.png",
        "ddg_heatmap.png",
    ]:
        plt.figure(figsize=(8, 4))
        plt.text(0.5, 0.5, "No numeric ddG values found", ha="center", va="center")
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(outdir / name, dpi=300)
        plt.close()
    raise SystemExit(0)


# Ranking plot: show all candidates, clipped for readability
# Ranking plot: show only top candidates to keep the figure readable
ranking_rows = rows[:100]

mutations = [r["mutation"] for r in ranking_rows]
ddgs = [clip_ddg(r["ddg_float"]) for r in ranking_rows]

plt.figure(figsize=(max(9, len(ranking_rows) * 0.25), 6))
plt.bar(mutations, ddgs)
plt.axhline(0, linestyle="--", linewidth=1)
plt.xticks(rotation=90)
plt.ylabel("ddG mutant - WT (clipped)")
plt.xlabel("Mutation")
plt.title("PyRosetta ddG ranking - top 100 candidates")
plt.tight_layout()
plt.savefig(outdir / "ddg_ranking_barplot.png", dpi=300)
plt.close()


# Top 10 stabilizing mutations
top = rows[:10]

plt.figure(figsize=(8, 5))
plt.barh(
    [r["mutation"] for r in reversed(top)],
    [r["ddg_float"] for r in reversed(top)],
)
plt.axvline(0, linestyle="--", linewidth=1)
plt.xlabel("ddG mutant - WT")
plt.ylabel("Mutation")
plt.title("Top stabilizing mutations")
plt.tight_layout()
plt.savefig(outdir / "top10_stabilizing_mutations.png", dpi=300)
plt.close()


# Heatmap: focus on interpretable stabilizing/desestabilizing range
positions = sorted(
    set(r["position"] for r in rows),
    key=lambda x: int(x) if str(x).isdigit() else x,
)

mutants = list("ACDEFGHIKLMNPQRSTVWY")

matrix = []
for pos in positions:
    row_values = []
    for aa in mutants:
        vals = [
            r["ddg_float"]
            for r in rows
            if r["position"] == pos and r["mutant"] == aa
        ]
        row_values.append(clip_ddg(vals[0]) if vals else float("nan"))
    matrix.append(row_values)

plt.figure(figsize=(12, max(4, len(positions) * 0.7)))
norm = TwoSlopeNorm(vmin=-25, vcenter=0, vmax=25)
plt.imshow(matrix, aspect="auto", cmap="coolwarm", norm=norm)
plt.colorbar(label="ddG mutant - WT")
plt.xticks(range(len(mutants)), mutants)
plt.yticks(range(len(positions)), positions)
plt.xlabel("Mutant amino acid")
plt.ylabel("Position")
plt.title("ddG heatmap")
plt.tight_layout()
plt.savefig(outdir / "ddg_heatmap.png", dpi=300)
plt.close()