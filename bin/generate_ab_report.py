#!/usr/bin/env python3
"""
generate_ab_report.py — HTML report para el pipeline de diseño de anticuerpos.

Lee:
  ranked_esmfold2.tsv    (rank, id, iptm_ht, iptm_global, pae_global,
                           pae_interface, plddt_global, plddt_cdr,
                           rmsd_global, rmsd_cdr, _score)
  ranked_af3.tsv         (rank, id, iptm, ptm, pae_global, plddt_global, _score)
  ranked_haddock3.tsv    (rank, id, haddock_score, vdw_energy,
                           elec_energy, buried_sasa, _score)
  pdbs_esm/*.cif         (ESMFold2 predicted complexes)
  pdbs_af3/*.cif         (AF3 predicted complexes)

Escribe: PhiPsiProt_antibody_report.html

Tres pestañas: ESMFold2 | AlphaFold3 | HADDOCK3.
Cada una: tabla ordenada por rank + gráfico de barras de score
+ dispersión de las dos métricas principales + visor NGL por candidato.
"""

import os, glob, json, base64, csv, argparse, sys


# ── Métrica de cada etapa ─────────────────────────────────────────────────────
STAGES = {
    "esmfold2": {
        "label":   "ESMFold2 pre-filtro",
        "metrics": [
            "iptm_ht", "iptm_global", "pae_global", "pae_interface",
            "plddt_global", "plddt_cdr", "rmsd_global", "rmsd_cdr"
        ],
        "labels": {
            "iptm_ht":      "ipTM H↔T",
            "iptm_global":  "ipTM global",
            "pae_global":   "PAE global",
            "pae_interface":"PAE interfaz",
            "plddt_global": "pLDDT global",
            "plddt_cdr":    "pLDDT CDR",
            "rmsd_global":  "RMSD global",
            "rmsd_cdr":     "RMSD CDR",
        },
        "scatter_x": "pae_interface",
        "scatter_y": "iptm_ht",
    },
    "af3": {
        "label":   "AlphaFold3 filtro fino",
        "metrics": ["iptm", "ptm", "pae_global", "plddt_global"],
        "labels": {
            "iptm":         "ipTM",
            "ptm":          "pTM",
            "pae_global":   "PAE global",
            "plddt_global": "pLDDT global",
        },
        "scatter_x": "pae_global",
        "scatter_y": "iptm",
    },
    "haddock3": {
        "label":   "HADDOCK3 docking",
        "metrics": ["haddock_score", "vdw_energy", "elec_energy", "buried_sasa"],
        "labels": {
            "haddock_score": "HADDOCK score",
            "vdw_energy":    "VdW",
            "elec_energy":   "Elec",
            "buried_sasa":   "BSA (Å²)",
        },
        "scatter_x": "haddock_score",
        "scatter_y": "buried_sasa",
    },
}


# ── Helpers ───────────────────────────────────────────────────────────────────
def b64_file(path):
    if not path or not os.path.exists(path):
        return ""
    with open(path, "rb") as fh:
        return base64.b64encode(fh.read()).decode()


def read_ranked(path):
    if not path or not os.path.exists(path):
        return []
    rows = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)
    def sort_key(r):
        try:
            return (0, int(r.get("rank", 9999)))
        except (ValueError, TypeError):
            return (1, 9999)
    return sorted(rows, key=sort_key)


def find_structure(struct_dir, cid):
    if not struct_dir or not os.path.isdir(struct_dir):
        return None
    for ext in (".cif", ".pdb"):
        for p in glob.glob(os.path.join(struct_dir, f"*{cid}*{ext}")):
            return p
        for p in glob.glob(os.path.join(struct_dir, f"{cid}{ext}")):
            return p
    return None


def build_stage_data(ranked_tsv, struct_dir, stage_key):
    cfg  = STAGES[stage_key]
    rows = read_ranked(ranked_tsv)
    cands = []
    for r in rows:
        cid = r.get("id", "")
        metrics = {}
        for m in cfg["metrics"]:
            v = r.get(m, "")
            try:
                metrics[m] = float(v) if v not in ("", "NA", None) else None
            except (ValueError, TypeError):
                metrics[m] = None
        try:
            score = float(r["_score"]) if r.get("_score") not in ("", "NA", None) else None
        except (ValueError, TypeError):
            score = None

        struct_path = find_structure(struct_dir, cid)
        ext         = "cif" if (struct_path or "").endswith(".cif") else "pdb"
        cands.append({
            "id":       cid,
            "rank":     r.get("rank", ""),
            "metrics":  metrics,
            "score":    score,
            "struct_b64":  b64_file(struct_path) if struct_path else "",
            "struct_ext":  ext,
            "struct_name": os.path.basename(struct_path) if struct_path else "",
        })
    return cands


def parse_args():
    p = argparse.ArgumentParser(description="PhiPsiProt antibody design report")
    p.add_argument("--ranked-esm",     default="ranked_esmfold2.tsv")
    p.add_argument("--ranked-af3",     default="ranked_af3.tsv")
    p.add_argument("--ranked-haddock", default="ranked_haddock3.tsv")
    p.add_argument("--pdbs-esm",       default="pdbs_esm")
    p.add_argument("--pdbs-af3",       default="pdbs_af3")
    p.add_argument("--out",            default="PhiPsiProt_antibody_report.html")
    return p.parse_args()


def main():
    args = parse_args()

    data = {
        "esmfold2": build_stage_data(args.ranked_esm,     args.pdbs_esm, "esmfold2"),
        "af3":      build_stage_data(args.ranked_af3,     args.pdbs_af3, "af3"),
        "haddock3": build_stage_data(args.ranked_haddock, None,          "haddock3"),
    }

    data_json   = json.dumps(data)
    stages_json = json.dumps(STAGES)

    html = f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1.0"/>
<title>PhiPsiProt — Antibody Design Report</title>
<link href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;600&family=Syne:wght@400;600;800&display=swap" rel="stylesheet">
<script src="https://cdn.rawgit.com/arose/ngl/v2.0.0-dev.37/dist/ngl.js"></script>
<script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
<style>
  :root{{--bg:#f6f8fc;--surface:#fff;--border:#d9e2ef;--accent:#2563eb;
         --accent2:#7c3aed;--accent3:#059669;--text:#0f172a;--muted:#475569;}}
  *{{box-sizing:border-box;margin:0;padding:0;}}
  body{{background:var(--bg);color:var(--text);font-family:'Syne',sans-serif;font-size:16px;line-height:1.5;}}
  .mono{{font-family:'JetBrains Mono',monospace;}}
  .nav{{display:flex;align-items:center;gap:1rem;padding:1.2rem 2rem;
        border-bottom:1px solid var(--border);background:#fff;position:sticky;top:0;z-index:50;}}
  .nav-title{{font-size:1.6rem;font-weight:800;
               background:linear-gradient(90deg,var(--accent),var(--accent2));
               -webkit-background-clip:text;-webkit-text-fill-color:transparent;}}
  .nav-sub{{font-size:.9rem;color:var(--muted);font-family:'JetBrains Mono',monospace;}}
  .tabs{{display:flex;gap:.5rem;padding:1rem 2rem 0;}}
  .tab{{font-family:'JetBrains Mono',monospace;font-weight:700;padding:.8rem 2rem;
        border-radius:.8rem .8rem 0 0;border:1px solid var(--border);border-bottom:none;
        background:#eef2ff;color:var(--muted);cursor:pointer;transition:all .2s;}}
  .tab.active{{background:#fff;color:var(--accent);}}
  .tab.stage-af3{{}}
  .tab.stage-haddock3{{}}
  .panel{{display:none;padding:1.5rem 2rem 3rem;}}
  .panel.active{{display:block;}}
  .summary-bar{{display:flex;gap:1rem;margin-bottom:1.5rem;flex-wrap:wrap;}}
  .stat{{background:#fff;border:1px solid var(--border);border-radius:.8rem;
          padding:.8rem 1.2rem;flex:1;min-width:120px;}}
  .stat-val{{font-size:2rem;font-weight:800;color:var(--accent);font-family:'JetBrains Mono',monospace;}}
  .stat-lbl{{font-size:.8rem;color:var(--muted);}}
  .charts{{display:grid;grid-template-columns:1fr 1fr;gap:1.5rem;margin-bottom:2rem;}}
  @media(max-width:900px){{.charts{{grid-template-columns:1fr;}}}}
  .chart-box{{background:#fff;border:1px solid var(--border);border-radius:1rem;padding:1rem;}}
  .chart-box h3{{font-size:.85rem;color:var(--muted);margin-bottom:.5rem;
                  font-family:'JetBrains Mono',monospace;text-transform:uppercase;letter-spacing:.05em;}}
  table.rank{{width:100%;border-collapse:collapse;background:#fff;
               border:1px solid var(--border);border-radius:1rem;overflow:hidden;
               margin-bottom:2rem;font-family:'JetBrains Mono',monospace;font-size:.78rem;}}
  table.rank th{{background:#f1f5f9;padding:.6rem .8rem;text-align:right;
                  color:var(--muted);border-bottom:1px solid var(--border);white-space:nowrap;}}
  table.rank th:first-child,table.rank th:nth-child(2){{text-align:left;}}
  table.rank td{{padding:.5rem .8rem;text-align:right;border-bottom:1px solid #eef2f7;}}
  table.rank td:first-child,table.rank td:nth-child(2){{text-align:left;}}
  table.rank tr.top td{{background:rgba(37,99,235,.06);font-weight:600;}}
  .badge{{display:inline-block;padding:.1rem .5rem;border-radius:9999px;font-size:.65rem;
           font-family:'JetBrains Mono',monospace;}}
  .badge.top{{background:rgba(37,99,235,.15);color:var(--accent);}}
  .card{{background:#fff;border:1px solid var(--border);border-radius:1rem;
          margin-bottom:1.5rem;overflow:hidden;}}
  .card-h{{display:flex;align-items:center;justify-content:space-between;
            padding:.8rem 1.2rem;border-bottom:1px solid var(--border);flex-wrap:wrap;gap:.5rem;}}
  .card-title{{font-weight:800;font-size:1rem;font-family:'JetBrains Mono',monospace;}}
  .pills{{display:flex;gap:.4rem;flex-wrap:wrap;}}
  .pill{{font-family:'JetBrains Mono',monospace;font-size:.68rem;padding:.2rem .5rem;
          border:1px solid var(--border);border-radius:.3rem;color:var(--muted);}}
  .ngl{{width:100%;height:400px;background:#050a14;}}
  .legend{{display:flex;gap:.8rem;padding:.6rem 1.2rem;flex-wrap:wrap;font-size:.75rem;
            border-top:1px solid var(--border);background:#fafbff;}}
  .ldot{{width:12px;height:12px;border-radius:3px;display:inline-block;
          margin-right:.3rem;vertical-align:middle;}}
  .pipeline-flow{{display:flex;align-items:center;gap:.5rem;margin-bottom:1.5rem;flex-wrap:wrap;}}
  .pf-stage{{background:#fff;border:1px solid var(--border);border-radius:.5rem;
              padding:.4rem .8rem;font-family:'JetBrains Mono',monospace;font-size:.8rem;}}
  .pf-stage.active{{background:#2563eb;color:#fff;border-color:#2563eb;}}
  .pf-arrow{{color:var(--muted);font-size:1.2rem;}}
</style>
</head>
<body>
<nav class="nav">
  <span class="nav-title">PhiPsiProt</span>
  <span class="nav-sub">Antibody Design Report</span>
</nav>
<div class="tabs" id="tabs"></div>
<div id="panels"></div>

<script>
const DATA   = {data_json};
const STAGES = {stages_json};
const STAGE_KEYS = ["esmfold2","af3","haddock3"];
const STAGE_COLORS = {{esmfold2:"#2563eb", af3:"#7c3aed", haddock3:"#059669"}};

function b64toBlob(b64,mime){{
  const bin=atob(b64),a=new Uint8Array(bin.length);
  for(let i=0;i<bin.length;i++)a[i]=bin.charCodeAt(i);
  return new Blob([a],{{type:mime}});
}}

function fmt(v,d=3){{
  if(v==null||v===undefined)return '—';
  return typeof v==='number'?v.toFixed(d):v;
}}

function plddtScheme(){{
  return NGL.ColormakerRegistry.addScheme(function(){{
    this.atomColor=function(atom){{
      const v=atom.bfactor;
      if(v>=90)return 0x2563eb; if(v>=70)return 0x06b6d4;
      if(v>=50)return 0xeab308; return 0xf97316;
    }};
  }});
}}

function buildTable(cands, stageKey){{
  const cfg = STAGES[stageKey];
  const metrics = cfg.metrics;
  const labels  = cfg.labels;
  const ths = metrics.map(m=>`<th>${{labels[m]||m}}</th>`).join('');
  const rows = cands.map(c=>{{
    const isTop = c.rank && parseInt(c.rank) <= 10;
    const cls   = isTop ? 'top' : '';
    const badge = isTop ? '<span class="badge top">★</span>' : '';
    const cells = metrics.map(m=>{{
      const v = c.metrics[m];
      const d = m.includes('sasa')||m.includes('energy')||m.includes('score')?1:
                m.includes('plddt')||m.includes('iptm')||m.includes('ptm')?4:3;
      return `<td>${{fmt(v,d)}}</td>`;
    }}).join('');
    return `<tr class="${{cls}}"><td>${{c.rank||'—'}}</td><td>${{c.id}} ${{badge}}</td>${{cells}}<td>${{fmt(c.score,4)}}</td></tr>`;
  }}).join('');
  return `<table class="rank">
    <thead><tr><th>#</th><th>ID</th>${{ths}}<th>score</th></tr></thead>
    <tbody>${{rows}}</tbody>
  </table>`;
}}

function renderCharts(stageKey, cands){{
  const cfg   = STAGES[stageKey];
  const color = STAGE_COLORS[stageKey];
  const scored = cands.filter(c=>c.score!=null);

  // Barras de score
  Plotly.newPlot('bars-'+stageKey, [{{
    type:'bar',
    x: scored.map(c=>c.id),
    y: scored.map(c=>c.score),
    marker:{{color: scored.map(c=>parseInt(c.rank)<=10?color:'#94a3b8')}},
    text: scored.map(c=>`#${{c.rank}}`),
    textposition:'outside',
    textfont:{{size:9,family:'JetBrains Mono'}},
  }}],{{
    margin:{{t:10,b:100,l:40,r:10}}, height:300,
    xaxis:{{tickangle:-50,tickfont:{{size:8,family:'JetBrains Mono'}}}},
    yaxis:{{title:'score normalizado'}},
    paper_bgcolor:'rgba(0,0,0,0)', plot_bgcolor:'rgba(0,0,0,0)',
  }},{{responsive:true,displayModeBar:false}});

  // Dispersión métricas principales
  const sx = cfg.scatter_x, sy = cfg.scatter_y;
  Plotly.newPlot('scatter-'+stageKey, [{{
    type:'scatter', mode:'markers+text',
    x: scored.map(c=>c.metrics[sx]),
    y: scored.map(c=>c.metrics[sy]),
    text: scored.map(c=>c.id),
    textposition:'top center',
    textfont:{{size:8,family:'JetBrains Mono'}},
    marker:{{
      size:12,
      color: scored.map(c=>parseInt(c.rank)<=10?color:'#94a3b8'),
      line:{{color:'#fff',width:1}}
    }},
  }}],{{
    margin:{{t:10,b:50,l:60,r:10}}, height:300,
    xaxis:{{title: cfg.labels[sx]||sx}},
    yaxis:{{title: cfg.labels[sy]||sy}},
    paper_bgcolor:'rgba(0,0,0,0)', plot_bgcolor:'rgba(0,0,0,0)',
  }},{{responsive:true,displayModeBar:false}});
}}

function initNGL(vid, c){{
  if(!c.struct_b64) return;
  const el = document.getElementById('ngl-'+vid);
  if(!el) return;
  const stage = new NGL.Stage(el, {{backgroundColor:'#050a14'}});
  new ResizeObserver(()=>stage.handleResize()).observe(el);
  const mime = c.struct_ext==='cif'?'chemical/x-cif':'chemical/x-pdb';
  stage.loadFile(b64toBlob(c.struct_b64,mime),{{ext:c.struct_ext}}).then(comp=>{{
    comp.addRepresentation('cartoon',{{color:plddtScheme()}});
    comp.autoView();
  }});
}}

// ── Build tabs + panels ────────────────────────────────────────────────────
const tabsEl   = document.getElementById('tabs');
const panelsEl = document.getElementById('panels');

STAGE_KEYS.forEach((key,i)=>{{
  const cfg   = STAGES[key];
  const cands = DATA[key] || [];
  const color = STAGE_COLORS[key];

  // Tab
  const tab = document.createElement('button');
  tab.className = 'tab' + (i===0?' active':'') + ' stage-'+key;
  tab.innerHTML = `${{cfg.label}} <span style="color:${{color}};font-size:.8rem">(${{cands.length}})</span>`;
  tab.onclick = ()=>{{
    document.querySelectorAll('.tab').forEach((t,j)=>t.classList.toggle('active',j===i));
    document.querySelectorAll('.panel').forEach(p=>p.classList.toggle('active',p.id==='panel-'+key));
  }};
  tabsEl.appendChild(tab);

  // Pipeline flow indicator
  const flowHtml = STAGE_KEYS.map((k,j)=>{{
    const n = (DATA[k]||[]).length;
    const isActive = k===key;
    return `<span class="pf-stage${{isActive?' active':''}}">${{STAGES[k].label}} (${{n}})</span>${{j<STAGE_KEYS.length-1?'<span class="pf-arrow">→</span>':''}}`;
  }}).join('');

  // Summary stats
  const top10  = cands.filter(c=>parseInt(c.rank)<=10).length;
  const hasStruct = cands.filter(c=>c.struct_b64).length;
  const avgScore = cands.filter(c=>c.score!=null).reduce((s,c)=>s+c.score,0) /
                   (cands.filter(c=>c.score!=null).length||1);

  // Panel
  const panel = document.createElement('div');
  panel.id = 'panel-'+key;
  panel.className = 'panel' + (i===0?' active':'');
  panel.innerHTML = `
    <div class="pipeline-flow">${{flowHtml}}</div>
    <div class="summary-bar">
      <div class="stat"><div class="stat-val">${{cands.length}}</div><div class="stat-lbl">diseños evaluados</div></div>
      <div class="stat"><div class="stat-val" style="color:${{color}}">${{top10}}</div><div class="stat-lbl">top 10 seleccionados</div></div>
      <div class="stat"><div class="stat-val">${{hasStruct}}</div><div class="stat-lbl">con estructura</div></div>
      <div class="stat"><div class="stat-val">${{isNaN(avgScore)?'—':avgScore.toFixed(3)}}</div><div class="stat-lbl">score medio</div></div>
    </div>
    <div class="charts">
      <div class="chart-box"><h3>Score por diseño</h3><div id="bars-${{key}}"></div></div>
      <div class="chart-box"><h3>${{cfg.labels[cfg.scatter_x]||cfg.scatter_x}} vs ${{cfg.labels[cfg.scatter_y]||cfg.scatter_y}}</h3><div id="scatter-${{key}}"></div></div>
    </div>
    ${{buildTable(cands, key)}}
    <div class="legend">
      <strong style="margin-right:.5rem">pLDDT:</strong>
      <span><span class="ldot" style="background:#2563eb"></span>≥90 (muy alto)</span>
      <span><span class="ldot" style="background:#06b6d4"></span>70–90 (confiable)</span>
      <span><span class="ldot" style="background:#eab308"></span>50–70 (bajo)</span>
      <span><span class="ldot" style="background:#f97316"></span>&lt;50 (muy bajo)</span>
    </div>
    <div id="cards-${{key}}"></div>`;
  panelsEl.appendChild(panel);

  // Tarjetas NGL por candidato (solo top estructuras)
  const cardsEl = panel.querySelector('#cards-'+key);
  cands.filter(c=>c.struct_b64).forEach(c=>{{
    const vid = (key+'_'+c.id).replace(/[^a-z0-9]/gi,'_');
    const pills = cfg.metrics.map(m=>{{
      const v = c.metrics[m];
      const d = m.includes('sasa')?0:m.includes('iptm')||m.includes('ptm')||m.includes('plddt')?4:3;
      return `<span class="pill">${{cfg.labels[m]||m}} ${{fmt(v,d)}}</span>`;
    }}).join('');
    const card = document.createElement('div');
    card.className = 'card';
    card.innerHTML = `
      <div class="card-h">
        <span class="card-title">#${{c.rank||'—'}} ${{c.id}}</span>
        <div class="pills">${{pills}}<span class="pill" style="border-color:${{color}};color:${{color}}">score ${{fmt(c.score,4)}}</span></div>
      </div>
      <div id="ngl-${{vid}}" class="ngl"></div>`;
    cardsEl.appendChild(card);
    setTimeout(()=>initNGL(vid,c), 80*cands.indexOf(c));
  }});

  setTimeout(()=>renderCharts(key,cands), 100);
}});
</script>
</body>
</html>"""

    with open(args.out, "w") as fh:
        fh.write(html)
    print(f"Report escrito: {args.out}")


if __name__ == "__main__":
    main()