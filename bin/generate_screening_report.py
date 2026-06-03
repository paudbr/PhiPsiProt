#!/usr/bin/env python3

import argparse
import csv
import base64
from pathlib import Path


def img_to_base64(path):
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode("utf-8")


parser = argparse.ArgumentParser()
parser.add_argument("--results_csv", required=True)
parser.add_argument("--plots_dir", required=True)
parser.add_argument("--outdir", required=True)
args = parser.parse_args()

results_csv = Path(args.results_csv)
plots_dir = Path(args.plots_dir)
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

rows = []
with open(results_csv, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            row["ddg_float"] = float(row["ddg"])
            rows.append(row)
        except Exception:
            pass

rows = sorted(rows, key=lambda r: r["ddg_float"])
top10 = rows[:10]

best = rows[0]["ddg_float"] if rows else "NA"
worst = rows[-1]["ddg_float"] if rows else "NA"
mean = sum(r["ddg_float"] for r in rows) / len(rows) if rows else "NA"

plot_files = [
    ("ddg_ranking_barplot.png", "PyRosetta ΔΔG ranking"),
    ("top10_stabilizing_mutations.png", "Top 10 stabilizing mutations"),
    ("ddg_heatmap.png", "ΔΔG heatmap"),
]

plots_html = ""
for filename, title in plot_files:
    path = plots_dir / filename
    if path.exists():
        b64 = img_to_base64(path)
        plots_html += f"""
        <section class="card">
            <h2>{title}</h2>
            <img src="data:image/png;base64,{b64}" />
        </section>
        """

table_rows = ""
for r in top10:
    table_rows += f"""
    <tr>
        <td>{r.get("candidate_id", "")}</td>
        <td>{r.get("mutation", "")}</td>
        <td>{r.get("wildtype", "")}</td>
        <td>{r.get("mutant", "")}</td>
        <td>{r.get("position", "")}</td>
        <td>{r.get("ddg", "")}</td>
        <td>{r.get("wt_score", "")}</td>
        <td>{r.get("mutant_score", "")}</td>
    </tr>
    """

html = f"""
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>PhiPsiProt Screening Report</title>
<style>
body {{
    font-family: Arial, sans-serif;
    margin: 0;
    background: #f5f7fb;
    color: #222;
}}
header {{
    background: linear-gradient(90deg, #1f2937, #334155);
    color: white;
    padding: 30px 45px;
}}
header h1 {{
    margin: 0;
    font-size: 34px;
}}
header p {{
    margin-top: 8px;
    color: #d1d5db;
}}
.container {{
    padding: 30px 45px;
}}
.summary {{
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 18px;
    margin-bottom: 28px;
}}
.metric {{
    background: white;
    border-radius: 14px;
    padding: 22px;
    box-shadow: 0 2px 10px rgba(0,0,0,0.08);
}}
.metric .label {{
    color: #64748b;
    font-size: 14px;
}}
.metric .value {{
    font-size: 28px;
    font-weight: bold;
    margin-top: 8px;
}}
.card {{
    background: white;
    border-radius: 14px;
    padding: 24px;
    margin-bottom: 28px;
    box-shadow: 0 2px 10px rgba(0,0,0,0.08);
}}
.card h2 {{
    margin-top: 0;
}}
img {{
    max-width: 100%;
    border: 1px solid #e5e7eb;
    border-radius: 10px;
}}
table {{
    width: 100%;
    border-collapse: collapse;
    margin-top: 12px;
}}
th {{
    background: #e5e7eb;
    text-align: left;
    padding: 10px;
}}
td {{
    border-bottom: 1px solid #e5e7eb;
    padding: 9px;
}}
.badge {{
    display: inline-block;
    background: #dcfce7;
    color: #166534;
    padding: 6px 10px;
    border-radius: 999px;
    font-weight: bold;
}}
footer {{
    color: #64748b;
    padding: 25px 45px;
    font-size: 13px;
}}
</style>
</head>
<body>

<header>
    <h1>PhiPsiProt Mutational Screening Report</h1>
    <p>PyRosetta ΔΔG screening and candidate ranking</p>
</header>

<div class="container">

    <div class="summary">
        <div class="metric">
            <div class="label">Candidates after filtering</div>
            <div class="value">{len(rows)}</div>
        </div>
        <div class="metric">
            <div class="label">Best ΔΔG</div>
            <div class="value">{best:.3f}</div>
        </div>
        <div class="metric">
            <div class="label">Worst ΔΔG</div>
            <div class="value">{worst:.3f}</div>
        </div>
        <div class="metric">
            <div class="label">Mean ΔΔG</div>
            <div class="value">{mean:.3f}</div>
        </div>
    </div>

    <section class="card">
        <h2>Top stabilizing mutations</h2>
        <p><span class="badge">Lower ΔΔG = more stabilizing</span></p>
        <table>
            <thead>
                <tr>
                    <th>Candidate</th>
                    <th>Mutation</th>
                    <th>WT</th>
                    <th>Mutant</th>
                    <th>Position</th>
                    <th>ΔΔG</th>
                    <th>WT score</th>
                    <th>Mutant score</th>
                </tr>
            </thead>
            <tbody>
                {table_rows}
            </tbody>
        </table>
    </section>

    {plots_html}

</div>

<footer>
    Generated by PhiPsiProt · Mode 1 Screening
</footer>

</body>
</html>
"""

with open(outdir / "screening_report.html", "w") as f:
    f.write(html)
