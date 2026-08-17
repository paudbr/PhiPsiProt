#!/usr/bin/env python3

"""
Generate the final HTML report for PhiPsiProt Mode 1 screening.

Purpose
-------
Create a self-contained HTML report from final_screening_results.csv,
screening plots and optional PyMOL structural images.

Current ranking scheme
----------------------
Candidates are ranked by PyRosetta ΔΔG.

final_score = ddg

Biophysical descriptors are reported as annotations only. They help interpret
candidate developability but do not modify the ranking score.
"""

import argparse
import base64
import csv
from pathlib import Path


def img_to_base64(path):
    """Encode image as base64 for embedding in HTML."""
    with open(path, "rb") as handle:
        return base64.b64encode(handle.read()).decode("utf-8")


def to_float(value):
    """Convert value to float, returning None for missing values."""
    try:
        if value in ["", "NA", None]:
            return None
        return float(value)
    except Exception:
        return None


def fmt_metric(value):
    """Format numeric metrics for display."""
    if isinstance(value, (int, float)):
        return f"{value:.3f}"
    return str(value)


def css_delta(value, good_if_positive=True):
    """Return CSS class for delta-style metrics."""
    v = to_float(value)
    if v is None:
        return ""

    if good_if_positive:
        return "good" if v > 0 else "bad" if v < 0 else ""

    return "good" if v < 0 else "bad" if v > 0 else ""


def css_ddg(value):
    """Return CSS class for ddG-like metrics."""
    v = to_float(value)
    if v is None:
        return ""

    if v < 0:
        return "good"
    if v > 0:
        return "bad"

    return ""


def top_candidate_summary(rows):
    """Create HTML summary for the top-ranked candidate."""
    if not rows:
        return ""

    top = rows[0]

    mutation = top.get("mutation", "NA")
    ddg = top.get("ddg", "NA")
    final_score = top.get("final_score", "NA")
    dsol = top.get("delta_solubility_score", "NA")
    dinst = top.get("delta_instability_index", "NA")

    return f"""
    <section class="card">
        <h2>Top candidate summary</h2>
        <p class="note">
            <b>{mutation}</b> is the highest ranked candidate.
            It has ΔΔG = <b>{ddg}</b>. In the current ranking scheme,
            final score = ΔΔG = <b>{final_score}</b>.
            Relative to WT, ΔSolubility = <b>{dsol}</b> and
            ΔInstability = <b>{dinst}</b>.
        </p>
    </section>
    """


def target_summary(rows):
    """Create HTML summary of the screened protein and residue selection."""
    if not rows:
        return ""

    first = rows[0]

    input_pdb = first.get("input_pdb", "NA")
    chain = first.get("chain", "NA")
    selection_mode = first.get("selection_mode", "NA")
    selected_positions = first.get("selected_positions", "NA")
    ligand_resname = first.get("ligand_resname", "")
    interface_chain = first.get("interface_chain", "")
    distance_cutoff = first.get("distance_cutoff", "NA")

    extra_rows = ""

    if ligand_resname:
        extra_rows += f"""
        <tr>
            <th>Ligand</th>
            <td>{ligand_resname}</td>
        </tr>
        """

    if interface_chain:
        extra_rows += f"""
        <tr>
            <th>Interface chain</th>
            <td>{interface_chain}</td>
        </tr>
        """

    return f"""
    <section class="card">
        <h2>Target protein and residue selection</h2>
        <table class="metadata-table">
            <tr>
                <th>Input structure</th>
                <td>{input_pdb}</td>
            </tr>
            <tr>
                <th>Target chain</th>
                <td>{chain}</td>
            </tr>
            <tr>
                <th>Selection mode</th>
                <td>{selection_mode}</td>
            </tr>
            <tr>
                <th>Selected positions</th>
                <td>{selected_positions}</td>
            </tr>
            {extra_rows}
            <tr>
                <th>Distance cutoff</th>
                <td>{distance_cutoff}</td>
            </tr>
        </table>
    </section>
    """


def build_pymol_html(pymol_dir, rows):
    """Create HTML block for optional PyMOL structural images."""
    if not pymol_dir:
        return ""

    target_img = pymol_dir / "target_structure.png"
    mutation_img = pymol_dir / "top_candidate_mutation.png"

    target_html = ""
    mutation_html = ""

    if target_img.exists():
        target_b64 = img_to_base64(target_img)
        target_html = f"""
        <div>
            <h3>Input structure</h3>
            <img src="data:image/png;base64,{target_b64}" />
        </div>
        """

    if mutation_img.exists():
        mutation_b64 = img_to_base64(mutation_img)
        mutation_label = rows[0].get("mutation", "Top candidate") if rows else "Top candidate"
        mutation_html = f"""
        <div>
            <h3>Top candidate: {mutation_label}</h3>
            <img src="data:image/png;base64,{mutation_b64}" />
        </div>
        """

    if not target_html and not mutation_html:
        return ""

    return f"""
    <section class="card">
        <h2>Structural visualization</h2>
        <p class="note">
            PyMOL renderings of the input protein structure and the top-ranked
            mutation site highlighted on the structure.
        </p>
        <div class="structure-grid">
            {target_html}
            {mutation_html}
        </div>
    </section>
    """


def build_plots_html(plots_dir):
    """Create HTML blocks for available screening plots."""
    plot_files = [
        ("ddg_ranking_barplot.png", "PyRosetta ΔΔG ranking"),
        ("top10_stabilizing_mutations.png", "Top 10 stabilizing mutations"),
        ("ddg_heatmap.png", "ΔΔG heatmap"),
    ]

    plots_html = ""

    for filename, title in plot_files:
        path = plots_dir / filename
        if path.exists():
            image_b64 = img_to_base64(path)
            plots_html += f"""
            <section class="card">
                <h2>{title}</h2>
                <img src="data:image/png;base64,{image_b64}" />
            </section>
            """

    return plots_html


def build_table_rows(rows):
    """Create HTML table rows for ranked candidates."""
    table_rows = ""

    for row in rows:
        table_rows += f"""
        <tr>
            <td>{row.get("candidate_id", "")}</td>
            <td>{row.get("mutation", "")}</td>
            <td>{row.get("position", "")}</td>
            <td>{row.get("wildtype", "")}</td>
            <td>{row.get("mutant", "")}</td>
            <td class="{css_ddg(row.get("ddg"))}">{row.get("ddg", "")}</td>
            <td class="{css_ddg(row.get("final_score"))}">{row.get("final_score", "")}</td>
            <td>{row.get("ddg_pass", "")}</td>
            <td class="{css_delta(row.get("delta_solubility_score"), True)}">{row.get("delta_solubility_score", "")}</td>
            <td class="{css_delta(row.get("delta_instability_index"), False)}">{row.get("delta_instability_index", "")}</td>
            <td class="{css_delta(row.get("delta_gravy"), False)}">{row.get("delta_gravy", "")}</td>
            <td>{row.get("isoelectric_point", "")}</td>
            <td>{row.get("molecular_weight", "")}</td>
        </tr>
        """

    return table_rows


def read_results(results_csv):
    """Read final screening results and keep candidates with numeric ddG."""
    rows = []

    with open(results_csv, newline="") as handle:
        reader = csv.DictReader(handle)

        for row in reader:
            ddg = to_float(row.get("ddg"))
            final_score = to_float(row.get("final_score"))

            if ddg is None:
                continue

            row["ddg_float"] = ddg
            row["final_score_float"] = final_score if final_score is not None else ddg
            rows.append(row)

    return sorted(rows, key=lambda row: row["final_score_float"])


def main():
    """Generate screening HTML report."""
    parser = argparse.ArgumentParser(
        description="Generate PhiPsiProt screening HTML report."
    )
    parser.add_argument("--results_csv", required=True)
    parser.add_argument("--plots_dir", required=True)
    parser.add_argument("--pymol_dir", required=False, default="")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    results_csv = Path(args.results_csv)
    plots_dir = Path(args.plots_dir)
    pymol_dir = Path(args.pymol_dir) if args.pymol_dir else None
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    rows = read_results(results_csv)

    best_ddg = min((row["ddg_float"] for row in rows), default="NA")
    mean_ddg = (
        sum(row["ddg_float"] for row in rows) / len(rows)
        if rows else "NA"
    )
    best_final_score = (
        rows[0]["final_score_float"]
        if rows else "NA"
    )

    table_rows = build_table_rows(rows)
    plots_html = build_plots_html(plots_dir)
    pymol_html = build_pymol_html(pymol_dir, rows)
    target_html = target_summary(rows)

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
.metadata-table th {{
    width: 220px;
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
.structure-grid {{
    display: grid;
    grid-template-columns: repeat(2, minmax(300px, 1fr));
    gap: 24px;
}}
.structure-grid h3 {{
    margin-top: 0;
    color: #334155;
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
    <p>PyRosetta ΔΔG screening with biophysical candidate annotation</p>
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

    {target_html}

    {top_candidate_summary(rows)}

    {pymol_html}

    <section class="card">
        <h2>Top ranked mutations</h2>
        <p><span class="badge">Lower final score = better PyRosetta ΔΔG</span></p>
        <p class="note">
            Candidates are ranked by PyRosetta ΔΔG. The final score is equal to
            the predicted ΔΔG. Biophysical descriptors are reported as annotations
            to help interpret candidate developability, but they are not used to
            modify the ranking.
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
                    <th>ddG pass</th>
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
            Green values indicate favorable changes. Red values indicate potentially
            unfavorable changes. Negative ΔΔG and lower final score are preferred.
            Positive ΔSolubility, negative ΔInstability and negative ΔGRAVY may
            indicate more favorable sequence-level properties, but these descriptors
            are reported separately from the ranking score.
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

    with open(outdir / "screening_report.html", "w") as handle:
        handle.write(html)


if __name__ == "__main__":
    main()