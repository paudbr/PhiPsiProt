#!/usr/bin/env python3

import argparse
import csv
import base64
from pathlib import Path


def img_to_base64(path):
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode("utf-8")


def to_float(value):
    try:
        if value in ["", "NA", None]:
            return None
        return float(value)
    except Exception:
        return None


def fmt_metric(value):
    if isinstance(value, (int, float)):
        return f"{value:.3f}"
    return str(value)

def css_delta(value, good_if_positive=True):
    v = to_float(value)
    if v is None:
        return ""
    if good_if_positive:
        return "good" if v > 0 else "bad" if v < 0 else ""
    return "good" if v < 0 else "bad" if v > 0 else ""

def css_ddg(value):
    v = to_float(value)
    if v is None:
        return ""
    if v < 0:
        return "good"
    if v > 0:
        return "bad"
    return ""


def css_final_score(value):
    v = to_float(value)
    if v is None:
        return ""
    if v < 0:
        return "good"
    if v > 0:
        return "bad"
    return ""

def top_candidate_summary(rows):
    if not rows:
        return ""

    r = rows[0]

    mutation = r.get("mutation", "NA")
    ddg = r.get("ddg", "NA")
    final_score = r.get("final_score", "NA")
    dsol = r.get("delta_solubility_score", "NA")
    dinst = r.get("delta_instability_index", "NA")

    return f"""
    <section class="card">
        <h2>Top candidate summary</h2>
        <p class="note">
            <b>{mutation}</b> is the highest ranked candidate.
            It has ΔΔG = <b>{ddg}</b> and final score = <b>{final_score}</b>.
            Relative to WT, ΔSolubility = <b>{dsol}</b> and
            ΔInstability = <b>{dinst}</b>.
        </p>
    </section>
    """

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
        ddg = to_float(row.get("ddg"))
        final_score = to_float(row.get("final_score"))

        if ddg is None:
            continue

        row["ddg_float"] = ddg
        row["final_score_float"] = final_score if final_score is not None else ddg
        rows.append(row)

# Final ranking uses the combined score when available
rows = sorted(rows, key=lambda r: r["final_score_float"])
top10 = rows[:10]
display_rows = rows


best_ddg = min((r["ddg_float"] for r in rows), default="NA")
worst_ddg = max((r["ddg_float"] for r in rows), default="NA")
mean_ddg = (
    sum(r["ddg_float"] for r in rows) / len(rows)
    if rows else "NA"
)
best_final_score = (
    rows[0]["final_score_float"]
    if rows else "NA"
)

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
for r in display_rows:
    table_rows += f"""
    <tr>
        <td>{r.get("candidate_id", "")}</td>
        <td>{r.get("mutation", "")}</td>
        <td>{r.get("position", "")}</td>
        <td>{r.get("wildtype", "")}</td>
        <td>{r.get("mutant", "")}</td>
        <td class="{css_ddg(r.get("ddg"))}">{r.get("ddg", "")}</td>
        <td class="{css_final_score(r.get("final_score"))}">{r.get("final_score", "")}</td>
        <td class="{css_delta(r.get("delta_solubility_score"), True)}">{r.get("delta_solubility_score", "")}</td>
        <td class="{css_delta(r.get("delta_instability_index"), False)}">{r.get("delta_instability_index", "")}</td>
        <td class="{css_delta(r.get("delta_gravy"), False)}">{r.get("delta_gravy", "")}</td>
        <td>{r.get("isoelectric_point", "")}</td>
        <td>{r.get("molecular_weight", "")}</td>
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
    font-size: 14px;
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

.table-container {{
    overflow-x: auto;
    overflow-y: auto;
    max-height: 700px;
    border-radius: 10px;
    border: 1px solid #e5e7eb;
}}

thead th {{
    position: sticky;
    top: 0;
    background: #e5e7eb;
    z-index: 2;
}}
.badge {{
    display: inline-block;
    background: #dcfce7;
    color: #166534;
    padding: 6px 10px;
    border-radius: 999px;
    font-weight: bold;
}}

.good {{
    color: #15803d;
    font-weight: bold;
}}

.bad {{
    color: #dc2626;
    font-weight: bold;
}}
.note {{
    color: #475569;
    line-height: 1.5;
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
    <p>PyRosetta ΔΔG screening with biophysical candidate prioritization</p>
</header>

<div class="container">

    <div class="summary">
        <div class="metric">
            <div class="label">Candidates ranked</div>
            <div class="value">{len(rows)}</div>
        </div>
        <div class="metric">
            <div class="label">Best ΔΔG</div>
            <div class="value">{fmt_metric(best_ddg)}</div>
        </div>
        <div class="metric">
            <div class="label">Mean ΔΔG</div>
            <div class="value">{fmt_metric(mean_ddg)}</div>
        </div>
        <div class="metric">
            <div class="label">Best final score</div>
            <div class="value">{fmt_metric(best_final_score)}</div>
        </div>
    </div>
    {top_candidate_summary(rows)}
    <section class="card">
        <h2>Top ranked mutations</h2>
        <p><span class="badge">Lower final score = better combined candidate</span></p>
        <p class="note">
            The final score combines PyRosetta ΔΔG with biophysical penalties,
            including instability and solubility. ΔΔG remains the main driver,
            while biophysical properties help prioritize more developable candidates.
        </p>

        <div class="table-container">
        <table>
            <thead>
                <tr>
                    <th>Candidate</th>
                    <th>Mutation</th>
                    <th>Position</th>
                    <th>WT</th>
                    <th>Mutant</th>
                    <th>ΔΔG</th>
                    <th>Final score</th>
                    <th>ΔSolubility</th>
                    <th>ΔInstability</th>
                    <th>ΔGRAVY</th>
                    <th>pI</th>
                    <th>MW</th>
                </tr>
            </thead>
            <tbody>
                {table_rows}
            </tbody>
        </table>
        </div>
    </section>

    <section class="card">
        <h2>Biophysical descriptors</h2>
        <p class="note">
            <b>Solubility score:</b> simple sequence-based proxy where higher values are better.<br>
            <b>Instability index:</b> values below 40 are generally considered more stable.<br>
            <b>GRAVY:</b> negative values indicate more hydrophilic sequences; positive values indicate more hydrophobic sequences.<br>
            <b>pI:</b> predicted isoelectric point of the mutated sequence.<br>
            <b>MW:</b> molecular weight of the mutated sequence.
        </p>
    </section>

    <section class="card">
        <h2>How to interpret this table</h2>
        <p class="note">
            Green values indicate favorable changes. Red values indicate potentially unfavorable changes.
            Negative ΔΔG and lower final score are preferred. Positive ΔSolubility is preferred.
            Negative ΔInstability and negative ΔGRAVY are generally preferred.
        </p>
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