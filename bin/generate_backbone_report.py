#!/usr/bin/env python3

import argparse
import base64
import csv
import json
import re
from pathlib import Path
from html import escape
from collections import defaultdict


def read_rows(csv_path):
    with open(csv_path, newline="") as handle:
        return list(csv.DictReader(handle))


def safe(value):
    if value is None or value == "":
        return "NA"
    return escape(str(value))


def to_float(value):
    try:
        return float(value)
    except Exception:
        return None


def mean_numeric(rows, column):
    values = [to_float(r.get(column)) for r in rows]
    values = [v for v in values if v is not None]
    if not values:
        return "NA"
    return round(sum(values) / len(values), 4)


def unique_count(rows, column):
    return len({r.get(column, "") for r in rows if r.get(column, "")})


def best_rows(rows, score_column="mpnn_score", n=10):
    def score(row):
        value = to_float(row.get(score_column))
        return value if value is not None else 999999
    return sorted(rows, key=score)[:n]


def best_row(rows):
    top = best_rows(rows, "mpnn_score", 1)
    return top[0] if top else {}


def color_class(column, value):
    v = to_float(value)
    if v is None:
        return ""
    if column == "mpnn_score":
        return "good" if v <= 0.9 else "warn" if v <= 1.0 else "bad"
    if column == "instability_index":
        return "good" if v < 40 else "bad"
    if column == "gravy":
        return "good" if v <= 0 else "warn"
    return ""


def make_table(rows, columns):
    html = ["<table>", "<thead><tr>"]
    for col in columns:
        html.append(f"<th>{safe(col)}</th>")
    html.append("</tr></thead><tbody>")
    for row in rows:
        html.append("<tr>")
        for col in columns:
            value = row.get(col, "NA")
            klass = color_class(col, value)
            if col in {"final_sequence", "biophysical_sequence", "designed_sequence"}:
                value = f"<code>{safe(value)}</code>"
            else:
                value = safe(value)
            html.append(f'<td class="{klass}">{value}</td>')
        html.append("</tr>")
    html.append("</tbody></table>")
    return "\n".join(html)


def mini_bar_plot(rows, column, label, n=15):
    values = []
    for row in rows[:n]:
        value = to_float(row.get(column))
        if value is not None:
            values.append((row.get("design_id", row.get("candidate_id", "NA")), value))
    if not values:
        return "<p>No numeric values available.</p>"
    max_value = max(v for _, v in values) or 1
    html = [f"<h3>{safe(label)}</h3>", '<div class="barplot">']
    for name, value in values:
        width = value / max_value * 100
        klass = "bar goodbar" if value <= 0.9 else "bar warnbar"
        html.append(f"""
            <div class="barrow">
                <div class="barlabel">{safe(name)}</div>
                <div class="barwrap"><div class="{klass}" style="width:{width:.1f}%"></div></div>
                <div class="barvalue">{value:.4f}</div>
            </div>
        """)
    html.append("</div>")
    return "\n".join(html)


def backbone_summary_table(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row.get("candidate_id", "NA")].append(row)
    summary = []
    for candidate_id, group in grouped.items():
        scores = [to_float(r.get("mpnn_score")) for r in group]
        scores = [s for s in scores if s is not None]
        best = min(group, key=lambda r: to_float(r.get("mpnn_score")) or 999999)
        summary.append({
            "candidate_id": candidate_id,
            "num_sequences": len(group),
            "best_design_id": best.get("design_id", "NA"),
            "best_mpnn_score": min(scores) if scores else "NA",
            "mean_mpnn_score": round(sum(scores) / len(scores), 4) if scores else "NA",
            "pdb": group[0].get("pdb", "NA"),
        })
    return summary


def b64_file(path):
    with open(path, "rb") as handle:
        return base64.b64encode(handle.read()).decode()


def collect_structures(structure_dir):
    structure_dir = Path(structure_dir)
    structures = []
    input_pdb = structure_dir / "input.pdb"
    if input_pdb.exists():
        structures.append({"id": "input", "label": "Input/reference structure", "filename": input_pdb.name, "b64": b64_file(input_pdb), "type": "input"})
    for pdb in sorted(structure_dir.glob("backbone_design_*.pdb")):
        structures.append({"id": pdb.stem, "label": pdb.stem, "filename": pdb.name, "b64": b64_file(pdb), "type": "design"})
    return structures


def parse_generated_region_from_contig(contig):
    if not contig or contig in {"NA", "null"}:
        return []
    contig = contig.strip().replace("[", "").replace("]", "").replace("'", "").replace('"', "")
    parts = contig.split("/")
    current_position = 1
    generated_regions = []
    for part in parts:
        part = part.strip()
        if not part:
            continue
        motif = re.match(r"^[A-Za-z](\d+)-(\d+)$", part)
        if motif:
            start = int(motif.group(1)); end = int(motif.group(2))
            current_position += abs(end - start) + 1
            continue
        gen = re.match(r"^(\d+)-(\d+)$", part)
        if gen:
            length = int(gen.group(2))
            start = current_position
            end = current_position + length - 1
            generated_regions.append([start, end])
            current_position += length
    return generated_regions


def parse_redesign_region(redesign_region):
    if not redesign_region or redesign_region in {"NA", "null"}:
        return []
    match = re.search(r"(\d+)\s*-\s*(\d+)", redesign_region)
    if not match:
        return []
    return [[int(match.group(1)), int(match.group(2))]]


def formatted_regions(regions, prefix=""):
    if not regions:
        return "NA"
    return ", ".join(f"{prefix}{start}-{end}" for start, end in regions)


def get_display_redesign_region(redesign_region, chain):
    regions = parse_redesign_region(redesign_region)
    if not regions:
        return "NA"
    return formatted_regions(regions, prefix=chain)


def sequence_region_summary(rows, redesign_region, chain):
    if not rows:
        return "<p>No candidate rows available.</p>"
    top = best_row(rows)
    region_label = get_display_redesign_region(redesign_region, chain)
    data = [{
        "Region requested by user": region_label,
        "Best design": top.get("design_id", "NA"),
        "Best backbone": top.get("candidate_id", "NA"),
        "ProteinMPNN score": top.get("mpnn_score", "NA"),
        "Sequence recovery": top.get("mpnn_seq_recovery", "NA"),
    }]
    columns = ["Region requested by user", "Best design", "Best backbone", "ProteinMPNN score", "Sequence recovery"]
    return make_table(data, columns)


def generate_report(rows, output_html, structure_dir=None, contig="NA", design_strategy="NA", redesign_region="NA", chain="A"):
    total_designs = len(rows)
    num_backbones = unique_count(rows, "candidate_id")
    num_sequences = unique_count(rows, "design_id")
    mean_mpnn = mean_numeric(rows, "mpnn_score")

    stable_count = sum(1 for r in rows if to_float(r.get("instability_index")) is not None and float(r["instability_index"]) < 40)
    hydrophilic_count = sum(1 for r in rows if to_float(r.get("gravy")) is not None and float(r["gravy"]) < 0)
    basic_count = sum(1 for r in rows if to_float(r.get("isoelectric_point")) is not None and float(r["isoelectric_point"]) > 7)

    top = best_row(rows)
    most_stable = min(rows, key=lambda r: to_float(r.get("instability_index")) or 999999) if rows else {}
    lowest_gravy = min(rows, key=lambda r: to_float(r.get("gravy")) or 999999) if rows else {}

    top_rows = best_rows(rows, "mpnn_score", 10)
    backbone_rows = backbone_summary_table(rows)

    top_columns = ["design_id", "candidate_id", "mpnn_score", "mpnn_seq_recovery", "gravy", "instability_index", "isoelectric_point", "final_sequence"]
    top_columns = [c for c in top_columns if rows and c in rows[0]]
    backbone_columns = ["candidate_id", "num_sequences", "best_design_id", "best_mpnn_score", "mean_mpnn_score", "pdb"]

    structures = collect_structures(structure_dir) if structure_dir else []
    generated_regions = parse_generated_region_from_contig(contig)
    input_redesign_regions = parse_redesign_region(redesign_region)

    structures_json = json.dumps(structures)
    generated_regions_json = json.dumps(generated_regions)
    input_regions_json = json.dumps(input_redesign_regions)

    display_redesign_region = get_display_redesign_region(redesign_region, chain)
    display_generated_region = formatted_regions(generated_regions)

    html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>PhiPsiProt - Mode 2 Backbone Stabilization Report</title>
<script src="https://cdn.rawgit.com/arose/ngl/v2.0.0-dev.37/dist/ngl.js"></script>
<style>
body {{ font-family: Arial, sans-serif; margin: 0; background: #f6f7fb; color: #222; }}
.page {{ max-width: 1300px; margin: 0 auto; padding: 34px; }}
.hero {{ background: linear-gradient(135deg, #102a43, #315f86); color: white; padding: 32px; border-radius: 18px; margin-bottom: 24px; box-shadow: 0 4px 14px rgba(0,0,0,0.16); }}
.hero h1 {{ margin: 0; font-size: 34px; }}
.hero p {{ margin-top: 8px; font-size: 15px; opacity: 0.92; }}
.card {{ background: white; padding: 22px; margin-bottom: 24px; border-radius: 14px; box-shadow: 0 2px 10px rgba(0,0,0,0.08); }}
h2 {{ color: #102a43; margin-top: 0; }}
h3 {{ color: #1f3b57; }}
.grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr)); gap: 14px; }}
.metric {{ background: #eef3f8; padding: 16px; border-radius: 12px; }}
.metric .value {{ font-size: 24px; font-weight: bold; color: #082032; word-break: break-word; }}
.metric .label {{ font-size: 13px; color: #555; }}
.summarybox {{ border-left: 5px solid #315f86; background: #f1f6fb; padding: 16px; border-radius: 10px; }}
table {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
th {{ background: #173b57; color: white; padding: 9px; text-align: left; }}
td {{ border-bottom: 1px solid #ddd; padding: 8px; vertical-align: top; }}
code {{ font-family: monospace; font-size: 12px; word-break: break-all; }}
.good {{ color: #11823b; font-weight: bold; }}
.warn {{ color: #c47f00; font-weight: bold; }}
.bad {{ color: #b00020; font-weight: bold; }}
.barplot {{ margin-top: 10px; }}
.barrow {{ display: grid; grid-template-columns: 230px 1fr 80px; gap: 10px; align-items: center; margin: 6px 0; }}
.barlabel {{ font-size: 12px; overflow: hidden; white-space: nowrap; text-overflow: ellipsis; }}
.barwrap {{ background: #edf0f5; height: 14px; border-radius: 8px; overflow: hidden; }}
.bar {{ height: 14px; }}
.goodbar {{ background: #2f80ed; }}
.warnbar {{ background: #f2994a; }}
.barvalue {{ font-size: 12px; }}
.viewer-layout {{ display: grid; grid-template-columns: 300px 1fr; gap: 18px; }}
.viewer-controls {{ background: #f1f6fb; border-radius: 12px; padding: 14px; }}
.viewer-controls label {{ display: block; font-size: 13px; font-weight: bold; color: #102a43; margin-top: 12px; margin-bottom: 6px; }}
.viewer-select {{ width: 100%; padding: 10px; border: 1px solid #cbd5e1; background: white; border-radius: 8px; font-size: 14px; }}
.viewer-button {{ width: 100%; display: block; margin-top: 10px; padding: 10px; border: 1px solid #315f86; background: white; color: #102a43; border-radius: 8px; cursor: pointer; text-align: center; }}
.viewer-button:hover {{ background: #eaf2fb; }}
.viewer-row {{ margin-bottom: 10px; }}
#ngl-viewer {{ width: 100%; height: 700px; border-radius: 14px; background: #020617; overflow: hidden; }}
.legend {{ margin-top: 14px; font-size: 13px; color: #475569; }}
.legend span {{ display: inline-block; margin-bottom: 6px; }}
.dot {{ width: 12px; height: 12px; border-radius: 999px; display: inline-block; margin-right: 4px; }}
.viewer-note {{ margin-top: 12px; font-size: 12px; color: #475569; line-height: 1.4; }}
.footer {{ font-size: 12px; color: #666; text-align: center; }}
</style>
</head>
<body>
<div class="page">
<div class="hero"><h1>PhiPsiProt Backbone Stabilization Report</h1><p>Mode 2 · RFdiffusion backbone redesign · ProteinMPNN sequence design</p></div>

<div class="card"><h2>Design strategy</h2><table>
<tr><th>Step</th><th>Description</th></tr>
<tr><td>Design strategy</td><td><code>{safe(design_strategy)}</code></td></tr>
<tr><td>Requested redesigned region</td><td><code>{safe(display_redesign_region)}</code></td></tr>
<tr><td>Internal RFdiffusion contig</td><td><code>{safe(contig)}</code></td></tr>
<tr><td>Generated region in RFdiffusion output numbering</td><td><code>{safe(display_generated_region)}</code></td></tr>
<tr><td>Backbone generation</td><td>RFdiffusion backbone redesign</td></tr>
<tr><td>Sequence design</td><td>ProteinMPNN sequences generated for each designed backbone</td></tr>
<tr><td>Biophysical annotation</td><td>Sequence-derived descriptors calculated for the final designed sequences</td></tr>
<tr><td>Total backbone candidates</td><td>{num_backbones}</td></tr>
<tr><td>Total ProteinMPNN sequences</td><td>{num_sequences}</td></tr>
<tr><td>Total final candidates</td><td>{total_designs}</td></tr>
</table></div>

<div class="card"><h2>Summary metrics</h2><div class="grid">
<div class="metric"><div class="value">{num_backbones}</div><div class="label">Backbone candidates</div></div>
<div class="metric"><div class="value">{num_sequences}</div><div class="label">ProteinMPNN sequences</div></div>
<div class="metric"><div class="value">{total_designs}</div><div class="label">Final candidates</div></div>
<div class="metric"><div class="value">{mean_mpnn}</div><div class="label">Mean ProteinMPNN score</div></div>
<div class="metric"><div class="value">{stable_count}/{total_designs}</div><div class="label">Instability Index &lt; 40</div></div>
<div class="metric"><div class="value">{safe(display_redesign_region)}</div><div class="label">Requested redesigned region</div></div>
</div></div>

<div class="card"><h2>Top candidate summary</h2><div class="summarybox">
<p><strong>{safe(top.get('design_id'))}</strong> is the best-ranked candidate according to the ProteinMPNN score (<strong>{safe(top.get('mpnn_score'))}</strong>).</p>
<p>Sequence recovery: <strong>{safe(top.get('mpnn_seq_recovery'))}</strong></p>
<p>Candidate with lowest Instability Index: <strong>{safe(most_stable.get('design_id'))}</strong> (Instability Index = {safe(most_stable.get('instability_index'))})</p>
<p>Candidate with lowest GRAVY value: <strong>{safe(lowest_gravy.get('design_id'))}</strong> (GRAVY = {safe(lowest_gravy.get('gravy'))})</p>
</div></div>

<div class="card"><h2>Designed region summary</h2><p>This table links the user-facing redesign request to the best-ranked ProteinMPNN sequence generated for the redesigned RFdiffusion backbone. ProteinMPNN designs sequences for RFdiffusion backbones; it does not generate new backbone coordinates.</p>{sequence_region_summary(rows, redesign_region, chain)}</div>
<div class="card"><h2>Backbone candidate summary</h2>{make_table(backbone_rows, backbone_columns)}</div>
<div class="card"><h2>Top candidates by ProteinMPNN score</h2><p>Lower ProteinMPNN scores are generally interpreted as better sequence compatibility with the designed backbone. Biophysical descriptors are shown as annotations and are not used here to modify the ProteinMPNN ranking.</p>{make_table(top_rows, top_columns)}</div>
<div class="card"><h2>ProteinMPNN score overview</h2>{mini_bar_plot(top_rows, 'mpnn_score', 'Top candidates')}</div>

<div class="card"><h2>Biophysical interpretation</h2><ul>
<li><strong>Instability Index:</strong> {stable_count}/{total_designs} candidates have Instability Index &lt; 40.</li>
<li><strong>GRAVY:</strong> {hydrophilic_count}/{total_designs} candidates have negative GRAVY values.</li>
<li><strong>Isoelectric point:</strong> {basic_count}/{total_designs} candidates have pI &gt; 7.</li>
</ul><p>Instability Index is interpreted using the classical threshold proposed by Guruprasad et al. (1990), where values below 40 are generally associated with stable proteins. GRAVY values are interpreted according to the hydropathy scale of Kyte & Doolittle (1982). These annotations are exploratory sequence-level descriptors and should not replace downstream structural validation.</p></div>

<div class="card"><h2>Structural redesign overview</h2><p>Interactive NGL viewer showing the input structure and RFdiffusion-designed backbone candidates. For the input/reference structure, the requested redesign region is highlighted. For RFdiffusion outputs, the generated region inferred from the contig is highlighted in orange.</p>
<div class="viewer-layout"><div class="viewer-controls">
<div class="viewer-row"><label for="structure-select">Structure</label><select id="structure-select" class="viewer-select" onchange="loadSelectedStructure()"></select></div>
<div class="viewer-row"><label for="representation-select">Representation</label><select id="representation-select" class="viewer-select" onchange="updateRepresentation()"><option value="cartoon">Cartoon</option><option value="backbone">Backbone</option><option value="ball+stick">Ball + stick</option><option value="licorice">Licorice</option><option value="surface">Surface</option></select></div>
<div class="viewer-row"><label for="color-select">Color mode</label><select id="color-select" class="viewer-select" onchange="updateRepresentation()"><option value="redesign">Highlight redesigned region</option><option value="default">Default</option><option value="chain">By chain</option><option value="residueindex">Residue index</option></select></div>
<button class="viewer-button" onclick="downloadCurrentStructure()">Download selected PDB</button><button class="viewer-button" onclick="resetView()">Reset view</button>
<div class="legend"><span><i class="dot" style="background:#9ca3af"></i>Input/reference structure</span><br><span><i class="dot" style="background:#60a5fa"></i>Designed backbone</span><br><span><i class="dot" style="background:#f97316"></i>Highlighted redesigned/generated region</span></div>
<div class="viewer-note">Requested region: <strong>{safe(display_redesign_region)}</strong><br>Generated output region: <strong>{safe(display_generated_region)}</strong></div>
</div><div id="ngl-viewer"></div></div></div>

<div class="card footer">Generated by PhiPsiProt · Mode 2 Backbone Stabilization</div>
</div>

<script>
const STRUCTURES = {structures_json};
const GENERATED_REGIONS = {generated_regions_json};
const INPUT_REDESIGN_REGIONS = {input_regions_json};
let stage = null;
let currentComponent = null;
let currentStructureIndex = 0;
function b64ToBlob(b64, mime) {{ const binary = atob(b64); const array = new Uint8Array(binary.length); for (let i = 0; i < binary.length; i++) {{ array[i] = binary.charCodeAt(i); }} return new Blob([array], {{type: mime}}); }}
function regionsToSelection(regions) {{ if (!regions || regions.length === 0) {{ return ""; }} return regions.map(r => r[0] + "-" + r[1]).join(" or "); }}
function getHighlightSelection(item) {{ if (!item) {{ return ""; }} if (item.type === "input") {{ return regionsToSelection(INPUT_REDESIGN_REGIONS); }} return regionsToSelection(GENERATED_REGIONS); }}
function getSelectedRepresentation() {{ const select = document.getElementById("representation-select"); return select ? select.value : "cartoon"; }}
function getSelectedColorMode() {{ const select = document.getElementById("color-select"); return select ? select.value : "redesign"; }}
function clearViewer() {{ if (stage) {{ stage.removeAllComponents(); }} currentComponent = null; }}
function addBaseRepresentation(component, item) {{ const representation = getSelectedRepresentation(); const colorMode = getSelectedColorMode(); let repParams = {{ color: item.type === "input" ? "#9ca3af" : "#60a5fa", opacity: 1.0 }}; if (colorMode === "chain") {{ repParams.color = "chainname"; }} if (colorMode === "residueindex") {{ repParams.color = "residueindex"; }} if (representation === "surface") {{ repParams.opacity = 0.62; }} component.addRepresentation(representation, repParams); }}
function addHighlightedRegionRepresentation(component, item) {{ const colorMode = getSelectedColorMode(); const selection = getHighlightSelection(item); if (!selection || colorMode !== "redesign") {{ return; }} component.addRepresentation("cartoon", {{ sele: selection, color: "#f97316", radiusScale: 1.6 }}); component.addRepresentation("ball+stick", {{ sele: selection, color: "#f97316" }}); component.addRepresentation("label", {{ sele: selection + " and .CA", labelType: "residue", color: "#111827", radius: 0.8, zOffset: 2.0 }}); }}
function renderCurrentStructure() {{ const item = STRUCTURES[currentStructureIndex]; if (!item || !stage) {{ return; }} clearViewer(); const blob = b64ToBlob(item.b64, "chemical/x-pdb"); stage.loadFile(blob, {{ext: "pdb"}}).then(component => {{ currentComponent = component; addBaseRepresentation(component, item); addHighlightedRegionRepresentation(component, item); component.autoView(); }}); }}
function loadStructure(index) {{ currentStructureIndex = parseInt(index); renderCurrentStructure(); }}
function loadSelectedStructure() {{ const select = document.getElementById("structure-select"); loadStructure(select.value); }}
function updateRepresentation() {{ renderCurrentStructure(); }}
function resetView() {{ if (currentComponent) {{ currentComponent.autoView(); }} }}
function downloadCurrentStructure() {{ const item = STRUCTURES[currentStructureIndex]; if (!item) {{ return; }} const blob = b64ToBlob(item.b64, "chemical/x-pdb"); const url = URL.createObjectURL(blob); const link = document.createElement("a"); link.href = url; link.download = item.filename || item.id + ".pdb"; document.body.appendChild(link); link.click(); document.body.removeChild(link); URL.revokeObjectURL(url); }}
function populateStructureSelect() {{ const select = document.getElementById("structure-select"); if (!select) {{ return; }} STRUCTURES.forEach((item, index) => {{ const option = document.createElement("option"); option.value = index; option.textContent = item.label; select.appendChild(option); }}); }}
function initViewer() {{ const viewer = document.getElementById("ngl-viewer"); if (!viewer || STRUCTURES.length === 0) {{ return; }} stage = new NGL.Stage("ngl-viewer", {{ backgroundColor: "white", quality: "high" }}); window.addEventListener("resize", function() {{ stage.handleResize(); }}); populateStructureSelect(); const colorSelect = document.getElementById("color-select"); if (colorSelect) {{ colorSelect.value = "redesign"; }} renderCurrentStructure(); }}
document.addEventListener("DOMContentLoaded", initViewer);
</script>
</body>
</html>
"""
    Path(output_html).write_text(html)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_csv", required=True)
    parser.add_argument("--output_html", required=True)
    parser.add_argument("--structure_dir", default=None)
    parser.add_argument("--contig", default="NA")
    parser.add_argument("--design_strategy", default="NA")
    parser.add_argument("--redesign_region", default="NA")
    parser.add_argument("--chain", default="A")
    args = parser.parse_args()
    rows = read_rows(args.input_csv)
    generate_report(rows=rows, output_html=args.output_html, structure_dir=args.structure_dir, contig=args.contig, design_strategy=args.design_strategy, redesign_region=args.redesign_region, chain=args.chain)
    print(f"Backbone report written to {args.output_html}")


if __name__ == "__main__":
    main()
