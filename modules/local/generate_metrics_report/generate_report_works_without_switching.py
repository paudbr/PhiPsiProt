#!/usr/bin/env python3
"""
PhiPsiProt – unified protein structure prediction metrics report.
Reads staged inputs:
  plddt/   *.tsv  — single-model:   chain  residue  plddt
                  — multi-model:    model  chain  residue  plddt   (*_allmodels_plddt.tsv)
  scores/  *.tsv  — simple:         metric  value
                  — ranked:         model  ptm  iptm  aggregate_score
  pdbs/    *.pdb / *.cif
Writes PhiPsiProt_metrics_report.html
"""
import os
import glob
import json
import base64
import re
from collections import defaultdict

# ── helpers ───────────────────────────────────────────────────────────────────

def read_tsv(path):
    rows = []
    with open(path) as fh:
        raw_header = fh.readline()
        if not raw_header:
            return rows
        header = [h.strip().lower() for h in raw_header.strip().split("\t")]
        for line in fh:
            line = line.strip()
            if not line:
                continue
            vals = line.split("\t")
            if len(vals) == len(header):
                rows.append(dict(zip(header, vals)))
    return rows


def b64_file(path):
    if not path or not os.path.exists(path):
        return ""
    with open(path, "rb") as fh:
        return base64.b64encode(fh.read()).decode()


TOOLS = {"af2", "af3", "esm", "chai1", "chai", "alphafold3", "esmfold"}

# Canonical names — collapse aliases to one key
TOOL_ALIASES = {
    "chai":        "chai1",
    "alphafold3":  "af3",
    "esmfold":     "esm",
}

def normalise_tool(t):
    return TOOL_ALIASES.get(t.lower(), t.lower())

def split_sample_tool(base):
    """Return (sample, tool) from a basename like Dsup_01_chai1 or Dsup_01."""
    parts = base.rsplit("_", 1)
    if len(parts) == 2 and parts[1].lower() in TOOLS:
        return parts[0], normalise_tool(parts[1])
    return base, "unknown"


def empty_entry(sample, tool):
    return {
        "models":          {},     # ranked_N → {plddt:[…], residues:[…]}
        "plddt":           [],     # ranked_0 / only model
        "residues":        [],
        "ptm":             None,
        "iptm":            None,
        "aggregate_score": None,
        "ranked_scores":   [],     # [{model, ptm, iptm, aggregate_score}, …]
        "pdb_b64":         "",
        "pdb_filename":    f"{sample}_{tool}.cif",
    }


# ── collect structure files ───────────────────────────────────────────────────

struct_files     = {}   # basename → path
struct_by_key    = {}   # (sample, tool) → path

for ext in ("pdb", "cif"):
    for p in glob.glob(f"pdbs/*.{ext}"):
        bn   = os.path.basename(p)
        base = bn.rsplit(".", 1)[0]
        struct_files[bn] = p
        s, t = split_sample_tool(base)
        for tool_key in [t, TOOL_ALIASES.get(t, t)]:
            struct_by_key[(s, tool_key)] = p

print(f"[DEBUG] pdbs found: {list(struct_files.keys())}")
print(f"[DEBUG] struct_by_key: {list(struct_by_key.keys())}")


# ── collect pLDDT TSVs ────────────────────────────────────────────────────────

data = {}   # data[sample][tool] = entry dict

for p in sorted(glob.glob("plddt/*.tsv")):
    fname = os.path.basename(p)

    is_multimodel = "_allmodels" in fname
    base = fname.replace("_allmodels_plddt.tsv", "").replace("_plddt.tsv", "")
    sample, tool = split_sample_tool(base)

    rows = read_tsv(p)
    if not rows:
        continue

    data.setdefault(sample, {})
    if tool not in data[sample]:
        data[sample][tool] = empty_entry(sample, tool)

    entry = data[sample][tool]

    resolved = None
    for ext in ("pdb", "cif"):
        key = f"{sample}_{tool}.{ext}"
        if key in struct_files:
            resolved = struct_files[key]
            break
    if not resolved:
        resolved = struct_by_key.get((sample, tool))
    if not resolved:
        for bn, p in struct_files.items():
            if sample in bn:
                resolved = p
                break
    if resolved:
        entry["pdb_b64"]      = b64_file(resolved)
        entry["pdb_filename"] = os.path.basename(resolved)
        print(f"[DEBUG] resolved struct for {sample}/{tool}: {resolved}")
    else:
        print(f"[WARN]  no struct found for {sample}/{tool}")

    if is_multimodel:
        by_model = defaultdict(list)
        for r in rows:
            by_model[r.get("model", "ranked_0")].append(r)

        def rank_key(name):
            m = re.search(r'\d+$', name)
            return int(m.group()) if m else 99

        for model_name in sorted(by_model, key=rank_key):
            mrows = by_model[model_name]
            entry["models"][model_name] = {
                "plddt":    [float(r["plddt"])  for r in mrows],
                "residues": [int(r["residue"])  for r in mrows],
            }

        best = entry["models"].get(
            "ranked_0",
            next(iter(entry["models"].values()), {})
        )
        entry["plddt"]    = best.get("plddt", [])
        entry["residues"] = best.get("residues", [])

    else:
        plddt    = [float(r["plddt"])  for r in rows]
        residues = [int(r["residue"])  for r in rows]
        entry["plddt"]    = plddt
        entry["residues"] = residues
        if "ranked_0" not in entry["models"]:
            entry["models"]["ranked_0"] = {"plddt": plddt, "residues": residues}


# ── collect score TSVs ────────────────────────────────────────────────────────

for p in sorted(glob.glob("scores/*.tsv")):
    fname = os.path.basename(p)
    base  = fname.replace("_scores.tsv", "")
    sample, tool = split_sample_tool(base)

    if sample not in data or tool not in data[sample]:
        data.setdefault(sample, {})
        if tool not in data[sample]:
            data[sample][tool] = empty_entry(sample, tool)

    entry = data[sample][tool]
    rows  = read_tsv(p)
    if not rows:
        continue

    cols = list(rows[0].keys())

    if "model" in cols:
        ranked = []
        for r in rows:
            rec = {"model": r["model"]}
            for k in ("ptm", "iptm", "aggregate_score", "plddt"):
                if k in r and r[k] not in ("", "NA", "None", "nan"):
                    try:
                        rec[k] = float(r[k])
                    except ValueError:
                        rec[k] = None
                else:
                    rec[k] = None

            if rec["plddt"] is None and rec["model"] in entry["models"]:
                plddt_list = entry["models"][rec["model"]]["plddt"]
                if plddt_list:
                    rec["plddt"] = sum(plddt_list) / len(plddt_list)

            ranked.append(rec)

        def rank_key(rec):
            m = re.search(r'\d+$', rec.get("model", ""))
            return int(m.group()) if m else 99

        ranked.sort(key=rank_key)
        entry["ranked_scores"] = ranked

        best = ranked[0] if ranked else {}
        entry["ptm"]             = best.get("ptm")
        entry["iptm"]            = best.get("iptm")
        entry["aggregate_score"] = best.get("aggregate_score")
        entry["mean_plddt"]      = best.get("plddt")

    elif "metric" in cols:
        scores = {}
        for r in rows:
            try:
                scores[r["metric"]] = float(r["value"])
            except (ValueError, KeyError):
                pass
        entry["ptm"]  = scores.get("ptm")
        entry["iptm"] = scores.get("iptm")

        if "plddt" in scores:
            entry["mean_plddt"] = scores.get("plddt")
        elif entry["plddt"]:
            entry["mean_plddt"] = sum(entry["plddt"]) / len(entry["plddt"])


print(f"[DEBUG] struct_by_key: {list(entry.keys())}")

# ── serialise ──────────────────────────────────────────────────────────────────

DATA_JSON = json.dumps(data)


# ── HTML ──────────────────────────────────────────────────────────────────────

HTML = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1.0"/>
  <title>PhiPsiProt Prediction Report</title>
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;600&family=Syne:wght@400;600;800&display=swap" rel="stylesheet">
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap-icons@1.11.3/font/bootstrap-icons.min.css"/>
  <script src="https://cdn.tailwindcss.com"></script>
  <script src="https://cdn.rawgit.com/arose/ngl/v2.0.0-dev.37/dist/ngl.js"></script>
  <script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
  <style>
    :root {
      --bg:       #f6f8fc;
      --surface:  #ffffff;
      --border:   #d9e2ef;

      --accent:   #2563eb;
      --accent2:  #7c3aed;

      --text:     #0f172a;
      --muted:    #475569;

      --ranked0:  #3b82f6;
      --ranked1:  #10b981;
      --ranked2:  #f59e0b;
      --ranked3:  #ef4444;
      --ranked4:  #8b5cf6;
    }

    * { box-sizing: border-box; margin: 0; padding: 0; }

    body {
      background:
        radial-gradient(circle at top left, #dbeafe, transparent 28%),
        radial-gradient(circle at top right, #ede9fe, transparent 24%),
        var(--bg);
      color: var(--text);
      font-family: 'Syne', sans-serif;
      min-height: 100vh;
      font-size: 18px;
      line-height: 1.6;
    }

    code, .mono { font-family: 'JetBrains Mono', monospace; }

    /* ── nav ── */
    .nav {
      display: flex; align-items: center; justify-content: space-between;
      padding: 1.5rem 3rem; border-bottom: 1px solid var(--border);
      background: rgba(255,255,255,.92); backdrop-filter: blur(16px);
      position: sticky; top: 0; z-index: 100;
    }
    .nav-brand { display: flex; align-items: center; gap: .75rem; }
    .nav-brand svg { width: 28px; height: 28px; }
    .nav-title {
      font-size: 2rem; font-weight: 800; letter-spacing: -.04em;
      background: linear-gradient(90deg, var(--accent), var(--accent2));
      -webkit-background-clip: text; -webkit-text-fill-color: transparent;
    }
    .nav-badge { font-family: 'JetBrains Mono', monospace; font-size: .65rem;
                 padding: .2rem .5rem; border-radius: 9999px;
                 border: 1px solid var(--border); color: var(--muted); }

    /* ── tool tabs ── */
    .tabs-wrap {
      display: flex; gap: 1rem; padding: 1.5rem 2rem 1rem;
      border-bottom: 1px solid var(--border); overflow-x: auto;
    }
    .tool-tab {
      font-family: 'JetBrains Mono', monospace; font-size: 1rem; font-weight: 700;
      letter-spacing: .08em; padding: 1.1rem 2.2rem; border-radius: 1rem;
      border: 1px solid var(--border); background: #ffffff; color: #0f172a;
      cursor: pointer; transition: all .25s ease; min-height: 64px;
      display: flex; align-items: center; justify-content: center;
    }
    .tool-tab:hover:not(.active) { background: #eef2ff; transform: translateY(-2px); }
    .tool-tab.active {
      background: linear-gradient(135deg, rgba(37,99,235,.14), rgba(124,58,237,.14));
      border-color: rgba(37,99,235,.35); color: #111827;
      box-shadow: 0 10px 30px rgba(37,99,235,.10), 0 4px 10px rgba(124,58,237,.08);
    }

    /* ── panels ── */
    .tool-panel        { display: none; padding: 1.5rem 2rem 3rem; }
    .tool-panel.active { display: block; }

    /* ── sample card ── */
    .sample-card {
      border: 1px solid var(--border); border-radius: 1.5rem;
      background: rgba(255,255,255,.92); backdrop-filter: blur(10px);
      overflow: hidden; margin-bottom: 3rem; transition: all .25s ease;
      box-shadow: 0 10px 30px rgba(15,23,42,.04);
    }
    .sample-card:hover { transform: translateY(-4px); box-shadow: 0 20px 40px rgba(15,23,42,.08); }

    .card-header {
      display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap;
      gap: .75rem; padding: 1rem 1.5rem; border-bottom: 1px solid var(--border);
      background: linear-gradient(90deg, rgba(15,39,68,.6), transparent);
    }
    .card-title { font-size: 1.5rem; font-weight: 800; color: #0f172a; }

    .sec-label { font-size: .85rem; letter-spacing: .12em; color: #64748b; }
    .card-tool  { font-family: 'JetBrains Mono', monospace; font-size: .7rem;
                  padding: .15rem .5rem; border-radius: 9999px;
                  background: rgba(56,189,248,.12); color: var(--accent);
                  border: 1px solid rgba(56,189,248,.25); }
    .metric-pill {
      font-family: 'JetBrains Mono', monospace; font-size: .7rem;
      padding: .2rem .6rem; border-radius: .25rem;
      background: rgba(255,255,255,.04); border: 1px solid var(--border);
      color: var(--muted);
    }
    .metric-pill strong { color: #e2e8f0; }
    .rep-btn, .model-btn, .dl-btn {
      font-size: .9rem !important; padding: .7rem 1rem !important; border-radius: .8rem !important;
    }
    .dl-btn:hover { background: rgba(255,255,255,.1); border-color: var(--accent); color: var(--accent); }

    /* ── two-col layout ── */
    .card-body { display: grid; grid-template-columns: 1fr 1fr; }
    @media (max-width: 900px) { .card-body { grid-template-columns: 1fr; } }
    .col-left  { padding: 1.25rem 1.5rem; border-right: 1px solid var(--border); display: flex; flex-direction: column; gap: 1rem; }
    .col-right { padding: 1.25rem 1.5rem; display: flex; flex-direction: column; gap: .75rem; }

    /* ── pLDDT legend ── */
    .legend { display: flex; flex-wrap: wrap; gap: 1rem; padding: 1rem 0; }
    .legend-item {
      display: flex; align-items: center; gap: .65rem;
      font-size: 1rem; font-weight: 600; color: #0f172a;
      background: rgba(255,255,255,.7); border: 1px solid #dbe4f0;
      padding: .7rem 1rem; border-radius: .9rem;
      box-shadow: 0 4px 10px rgba(15,23,42,.04);
    }
    .legend-dot { width: 18px; height: 18px; border-radius: .35rem; flex-shrink: 0; }

    /* ── model selector ── */
    .model-selector { display: flex; flex-wrap: wrap; gap: .4rem; }
    .model-btn {
      font-family: 'JetBrains Mono', monospace; font-size: .65rem;
      padding: .2rem .55rem; border-radius: .25rem;
      border: 1px solid var(--border); background: transparent;
      color: var(--muted); cursor: pointer; transition: all .15s;
    }
    .model-btn.active  { border-color: var(--accent); color: var(--accent); background: rgba(56,189,248,.08); }
    .model-btn:hover:not(.active) { border-color: var(--muted); color: var(--text); }

    /* ── ranked scores table ── */
    .ranked-wrap { overflow-x: auto; border-radius: .375rem; border: 1px solid var(--border); }
    .ranked-table { width: 100%; border-collapse: collapse; font-family: 'JetBrains Mono', monospace; font-size: .7rem; }
    .ranked-table th { padding: .35rem .75rem; text-align: left; color: var(--muted);
                       font-weight: 600; border-bottom: 1px solid var(--border);
                       background: rgba(255,255,255,.02); }
    .ranked-table td { padding: .3rem .75rem; border-bottom: 1px solid rgba(30,45,69,.5); }
    .ranked-table tr:last-child td { border-bottom: none; }
    .ranked-table tr.best td { color: #e2e8f0; }
    .ranked-table tr:not(.best) td { color: var(--muted); }
    .rank-dot { display: inline-block; width: 7px; height: 7px; border-radius: 50%; margin-right: 5px; }

    /* ── rep buttons ── */
    .rep-btns { display: flex; flex-wrap: wrap; gap: .35rem; }
    .rep-btn {
      font-size: .7rem; padding: .25rem .65rem; border-radius: 9999px;
      border: 1px solid var(--border); background: transparent;
      color: var(--muted); cursor: pointer; transition: all .15s;
    }
    .rep-btn.active  { background: rgba(56,189,248,.15); border-color: var(--accent); color: var(--accent); }
    .rep-btn:hover:not(.active) { border-color: var(--muted); color: var(--text); }

    /* ── NGL viewer ── */
    .ngl-stage { width: 100%; height: 520px; border-radius: 1rem; background: #020617; overflow: hidden; }
    .ngl-offscreen { position: absolute; visibility: hidden; top: -9999px; left: -9999px; width: 1px; height: 1px; }

    /* ── section label ── */
    .sec-label { font-size: .65rem; letter-spacing: .08em; text-transform: uppercase;
                 color: var(--muted); font-family: 'JetBrains Mono', monospace; }
  </style>
</head>
<body>

<nav class="nav">
  <div class="nav-brand">
    <svg viewBox="0 0 28 28" fill="none" xmlns="http://www.w3.org/2000/svg">
      <circle cx="14" cy="14" r="13" stroke="#38bdf8" stroke-width="1.5"/>
      <path d="M7 14 Q10 7 14 14 Q18 21 21 14" stroke="#818cf8" stroke-width="2" fill="none" stroke-linecap="round"/>
      <circle cx="14" cy="14" r="2.5" fill="#38bdf8"/>
    </svg>
    <span class="nav-title">PhiPsiProt</span>
    <span class="nav-badge mono">Prediction Report</span>
  </div>
</nav>

<div class="tabs-wrap" id="tool-tabs"></div>
<div id="tool-panels"></div>
<canvas id="ngl-off" class="ngl-offscreen"></canvas>

<script>
const DATA    = """ + DATA_JSON + """;
const TOOLS   = [...new Set(Object.values(DATA).flatMap(s => Object.keys(s)))].sort();
const SAMPLES = Object.keys(DATA).sort();

const RANK_COLORS = ['#38bdf8','#34d399','#f59e0b','#f87171','#a78bfa'];
const REPS        = ['cartoon','ball+stick','surface','licorice'];
const REP_LABELS  = {cartoon:'Cartoon','ball+stick':'Ball+Stick',surface:'Surface',licorice:'Ligands'};

const viewerState = {};
const offStage    = new NGL.Stage('ngl-off', {backgroundColor:'black', quality:'low'});

// ── utilities ─────────────────────────────────────────────────────────────────
function b64toBlob(b64, mime) {
  const bin = atob(b64), arr = new Uint8Array(bin.length);
  for (let i=0; i<bin.length; i++) arr[i] = bin.charCodeAt(i);
  return new Blob([arr], {type: mime});
}

function guessExt(filename) {
  return filename && filename.endsWith('.cif') ? 'cif' : 'pdb';
}

function downloadStruct(b64, filename) {
  if (!b64) return;
  const ext  = guessExt(filename);
  const mime = ext === 'cif' ? 'chemical/x-cif' : 'chemical/x-pdb';
  const url  = URL.createObjectURL(b64toBlob(b64, mime));
  const a    = Object.assign(document.createElement('a'), {href:url, download:filename});
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

// ── pLDDT color ───────────────────────────────────────────────────────────────
function plddtColor(v) {
  if (v >= 90) return '#2563eb';
  if (v >= 70) return '#06b6d4';
  if (v >= 50) return '#eab308';
  return '#f97316';
}

// ── Plotly traces — one per ranked model ─────────────────────────────────────
function buildTraces(models, activeModel) {
  return Object.entries(models).map(([name, m], i) => {
    const isActive = name === activeModel;
    return {
      x: m.residues, y: m.plddt,
      mode: 'lines', name: name,
      line: {color: RANK_COLORS[i % RANK_COLORS.length], width: isActive ? 2.5 : 1},
      opacity: isActive ? 1 : 0.35,
      hovertemplate: `<b>${name}</b><br>Residue %{x}<br>pLDDT: %{y:.1f}<extra></extra>`,
    };
  });
}

const PLOTLY_LAYOUT = {
  paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
  font: {color:'#64748b', size:10, family:"'JetBrains Mono', monospace"},
  margin: {t:10, b:40, l:50, r:10},
  xaxis: {title:'Residue', gridcolor:'#1e2d45', zerolinecolor:'#1e2d45', tickfont:{size:9}},
  yaxis: {title:'pLDDT', range:[0,100], gridcolor:'#1e2d45', zerolinecolor:'#1e2d45', tickfont:{size:9}},
  legend: {bgcolor:'rgba(0,0,0,0)', font:{size:9}, orientation:'h', y:-0.2},
  shapes: [
    {type:'rect', x0:0, x1:1, xref:'paper', y0:90, y1:100, fillcolor:'rgba(37,99,235,.06)', line:{width:0}},
    {type:'rect', x0:0, x1:1, xref:'paper', y0:70, y1:90,  fillcolor:'rgba(6,182,212,.06)',  line:{width:0}},
    {type:'rect', x0:0, x1:1, xref:'paper', y0:50, y1:70,  fillcolor:'rgba(234,179,8,.06)',  line:{width:0}},
    {type:'rect', x0:0, x1:1, xref:'paper', y0:0,  y1:50,  fillcolor:'rgba(249,115,22,.06)', line:{width:0}},
  ],
};

// ── ranked scores table ───────────────────────────────────────────────────────
function buildRankedTable(ranked) {
  if (!ranked || ranked.length === 0) return '';
  const metricCols = ['plddt','ptm','iptm','aggregate_score'].filter(k => ranked.some(r => r[k] != null));
  if (metricCols.length === 0) return '';

  const headers = ['<th>Model</th>', ...metricCols.map(c => `<th>${c.replace('_',' ').toUpperCase()}</th>`)].join('');
  const rows = ranked.map((r, i) => {
    const dot   = `<span class="rank-dot" style="background:${RANK_COLORS[i % RANK_COLORS.length]}"></span>`;
    const cells = metricCols.map(k => `<td>${r[k] != null ? r[k].toFixed(4) : '—'}</td>`).join('');
    return `<tr class="${i===0?'best':''}"><td>${dot}${r.model}</td>${cells}</tr>`;
  }).join('');

  return `<div>
    <p class="sec-label" style="margin-bottom:.5rem">Ranked model scores</p>
    <div class="ranked-wrap">
      <table class="ranked-table">
        <thead><tr>${headers}</tr></thead>
        <tbody>${rows}</tbody>
      </table>
    </div>
  </div>`;
}

// ── NGL: custom pLDDT colour scheme ──────────────────────────────────────────
// Registers a per-residue colour scheme using the 4-band legend palette.
// plddtArray: Float array indexed by residueIndex (0-based within the component).
function makePlddtScheme(plddtArray) {
  return NGL.ColormakerRegistry.addScheme(function() {
    this.atomColor = function(atom) {
      const val = plddtArray[atom.residueIndex];
      if (val === undefined || val === null) return 0xaaaaaa;
      if (val >= 90) return 0x2563eb;   // Very high  — blue
      if (val >= 70) return 0x06b6d4;   // Confident  — cyan
      if (val >= 50) return 0xeab308;   // Low        — yellow
      return 0xf97316;                   // Very low   — orange
    };
  });
}

// ── NGL helpers ───────────────────────────────────────────────────────────────
function clearReps(comp) { comp.removeAllRepresentations(); }

function addRep(comp, name, plddtArray) {
  // Use custom pLDDT scheme when data is available, otherwise fall back to bfactor
  const scheme = (plddtArray && plddtArray.length) ? makePlddtScheme(plddtArray) : 'bfactor';
  if (name==='cartoon')    comp.addRepresentation('cartoon',    {color: scheme});
  if (name==='ball+stick') comp.addRepresentation('ball+stick', {color: scheme});
  if (name==='surface')    comp.addRepresentation('surface',    {opacity:.7, color: scheme});
  if (name==='licorice')   comp.addRepresentation('licorice',   {sele:'hetero and not water', color: scheme});
}

function applyReps(comp, activeReps, plddtArray) {
  clearReps(comp);
  activeReps.forEach(r => addRep(comp, r, plddtArray));
}

// ── build UI ──────────────────────────────────────────────────────────────────
const tabContainer   = document.getElementById('tool-tabs');
const panelContainer = document.getElementById('tool-panels');

TOOLS.forEach((tool, tIdx) => {
  // tab
  const tab = document.createElement('button');
  tab.className = 'tool-tab' + (tIdx===0 ? ' active' : '');
  tab.textContent = tool.toUpperCase();
  tab.onclick = () => switchTool(tool);
  tabContainer.appendChild(tab);

  // panel
  const panel = document.createElement('div');
  panel.id = 'panel-' + tool;
  panel.className = 'tool-panel' + (tIdx===0 ? ' active' : '');

  SAMPLES.forEach(sample => {
    const d = (DATA[sample]||{})[tool];
    if (!d) return;

    const vid = (tool + '_' + sample).replace(/[^a-z0-9]/gi,'_');
    // activePlddt: pLDDT array for the currently shown model (used to colour NGL)
    viewerState[vid] = {stage:null, component:null, activeReps:new Set(['cartoon']), activeModel:'ranked_0', activePlddt:[]};

    const meanPlddt    = d.plddt.length ? (d.plddt.reduce((a,b)=>a+b,0)/d.plddt.length).toFixed(1) : '—';
    const ptmStr       = d.ptm             != null ? d.ptm.toFixed(3)             : '—';
    const iptmStr      = d.iptm            != null ? d.iptm.toFixed(3)            : '—';
    const aggStr       = d.aggregate_score != null ? d.aggregate_score.toFixed(4) : null;

    const modelNames   = Object.keys(d.models || {});
    const multiModel   = modelNames.length > 1;

    const modelBtns = multiModel ? modelNames.map((name, i) =>
      `<button class="model-btn${i===0?' active':''}" data-vid="${vid}" data-model="${name}"
               onclick="switchModel('${vid}','${name}')"
               style="${i===0?'border-color:'+RANK_COLORS[0]+';color:'+RANK_COLORS[0]+'':''}">
        ${name}
       </button>`
    ).join('') : '';

    const repBtns = REPS.map(r =>
      `<button id="rep-${vid}-${r.replace('+','-')}" class="rep-btn${r==='cartoon'?' active':''}"
               onclick="toggleRep('${vid}','${r}')">${REP_LABELS[r]}</button>`
    ).join('');

    const card = document.createElement('div');
    card.className = 'sample-card';
    card.innerHTML = `
      <div class="card-header">
        <div style="display:flex;align-items:center;gap:.75rem;">
          <span class="card-title">${sample}</span>
          <span class="card-tool">${tool.toUpperCase()}</span>
        </div>
        <div style="display:flex;align-items:center;gap:.5rem;flex-wrap:wrap;">
          <span class="metric-pill">pLDDT <strong>${meanPlddt}</strong></span>
          <span class="metric-pill">pTM <strong>${ptmStr}</strong></span>
          <span class="metric-pill">ipTM <strong>${iptmStr}</strong></span>
          ${aggStr ? `<span class="metric-pill">agg <strong>${aggStr}</strong></span>` : ''}
          ${d.pdb_b64 ? `<button class="dl-btn" onclick="downloadStruct('${d.pdb_b64}','${d.pdb_filename}')">
            <i class="bi bi-download"></i> ${d.pdb_filename.split('.').pop().toUpperCase()}
          </button>` : ''}
        </div>
      </div>
      <div class="card-body">
        <div class="col-left">
          <div class="legend">
            <span class="legend-item"><span class="legend-dot" style="background:#2563eb"></span>Very high (&ge;90)</span>
            <span class="legend-item"><span class="legend-dot" style="background:#06b6d4"></span>Confident (70&ndash;90)</span>
            <span class="legend-item"><span class="legend-dot" style="background:#eab308"></span>Low (50&ndash;70)</span>
            <span class="legend-item"><span class="legend-dot" style="background:#f97316"></span>Very low (&lt;50)</span>
          </div>
          ${multiModel ? `<div><p class="sec-label" style="margin-bottom:.4rem">Model overlay</p>
            <div class="model-selector" id="msel-${vid}">${modelBtns}</div></div>` : ''}
          <div id="plot-${vid}" style="height:420px;"></div>
          <div id="ranked-${vid}"></div>
        </div>
        <div class="col-right">
          <div><p class="sec-label" style="margin-bottom:.4rem">Representation</p>
            <div class="rep-btns">${repBtns}</div>
          </div>
          <div id="ngl-${vid}" class="ngl-stage"></div>
        </div>
      </div>`;

    panel.appendChild(card);
    setTimeout(() => initCard(vid, d), 60);
  });

  panelContainer.appendChild(panel);
});

// ── switchTool ────────────────────────────────────────────────────────────────
function switchTool(tool) {
  document.querySelectorAll('.tool-tab').forEach((t,i)   => t.classList.toggle('active', TOOLS[i]===tool));
  document.querySelectorAll('.tool-panel').forEach(p     => p.classList.toggle('active', p.id==='panel-'+tool));
}

// ── initCard ──────────────────────────────────────────────────────────────────
function initCard(vid, d) {
  const state  = viewerState[vid];
  const models = d.models && Object.keys(d.models).length ? d.models
                 : {ranked_0: {plddt: d.plddt, residues: d.residues}};

  // Seed activePlddt with ranked_0 (best model) data
  const bestModel = models['ranked_0'] || models[Object.keys(models)[0]] || {};
  state.activePlddt = bestModel.plddt || d.plddt || [];

  // plot
  Plotly.newPlot('plot-'+vid, buildTraces(models, state.activeModel),
    Object.assign({}, PLOTLY_LAYOUT), {responsive:true, displayModeBar:false});

  // ranked table
  document.getElementById('ranked-'+vid).innerHTML = buildRankedTable(d.ranked_scores);

  // NGL
  if (d.pdb_b64) {
    const stageEl = document.getElementById('ngl-'+vid);
    const stage   = new NGL.Stage(stageEl, {backgroundColor:'#050a14', quality:'medium'});
    state.stage   = stage;
    new ResizeObserver(() => stage.handleResize()).observe(stageEl);

    const ext  = guessExt(d.pdb_filename);
    const mime = ext==='cif' ? 'chemical/x-cif' : 'chemical/x-pdb';
    stage.loadFile(b64toBlob(d.pdb_b64, mime), {ext}).then(comp => {
      state.component = comp;
      applyReps(comp, state.activeReps, state.activePlddt);
      comp.autoView();
    });
  }
}

// ── switchModel ───────────────────────────────────────────────────────────────
function switchModel(vid, modelName) {
  const state = viewerState[vid];
  state.activeModel = modelName;
  const d = (() => {
    for (const sample of SAMPLES) for (const tool of TOOLS) {
      const entry = (DATA[sample]||{})[tool];
      if (entry && (tool+'_'+sample).replace(/[^a-z0-9]/gi,'_') === vid) return entry;
    }
    return null;
  })();
  if (!d) return;

  const models = d.models && Object.keys(d.models).length ? d.models
                 : {ranked_0:{plddt:d.plddt, residues:d.residues}};

  // Update active pLDDT data for the chosen model → re-colour NGL viewer
  const chosenModel = models[modelName] || {};
  state.activePlddt = chosenModel.plddt || d.plddt || [];

  Plotly.react('plot-'+vid, buildTraces(models, modelName),
    Object.assign({}, PLOTLY_LAYOUT), {responsive:true, displayModeBar:false});

  // Re-colour the structure with the new model's pLDDT values
  if (state.component) {
    applyReps(state.component, state.activeReps, state.activePlddt);
  }

  // update button styles
  const msel = document.getElementById('msel-'+vid);
  if (msel) {
    const modelNames = Object.keys(models);
    msel.querySelectorAll('.model-btn').forEach(btn => {
      const name = btn.dataset.model;
      const idx  = modelNames.indexOf(name);
      const col  = RANK_COLORS[idx % RANK_COLORS.length];
      const active = name === modelName;
      btn.classList.toggle('active', active);
      btn.style.borderColor = active ? col : '';
      btn.style.color       = active ? col : '';
    });
  }
}

// ── toggleRep ─────────────────────────────────────────────────────────────────
function toggleRep(vid, repName) {
  const state = viewerState[vid];
  if (!state || !state.component) return;
  state.activeReps.has(repName) ? state.activeReps.delete(repName) : state.activeReps.add(repName);
  applyReps(state.component, state.activeReps, state.activePlddt);
  REPS.forEach(r => {
    const btn = document.getElementById('rep-'+vid+'-'+r.replace('+','-'));
    if (btn) btn.classList.toggle('active', state.activeReps.has(r));
  });
}
</script>
</body>
</html>
"""

with open("PhiPsiProt_metrics_report.html", "w") as fh:
    fh.write(HTML)

print("Report written: PhiPsiProt_metrics_report.html")
