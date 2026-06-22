#!/usr/bin/env python3
"""
generate_cascade_report.py — report HTML para la cascada de cribado.

Lee:
  metrics_esm/*_metrics.tsv    (formato largo: id  metric  value)
  metrics_chai/*_metrics.tsv
  ranked_esm.csv               (rank, id, ...metrics, score, passed_gates, in_top_n)
  ranked_chai.csv
  pdbs_esm/*.pdb               (modelos de ESMFold, B-factor = pLDDT)
  pdbs_chai/*.pdb              (modelos de Chai, convertidos a pdb)

Escribe: PhiPsiProt_cascade_report.html

Dos pestañas: ESM y CHAI. Cada una: tabla ordenada por rank con las 6 métricas
+ score, gráfico de barras de score, dispersión ΔG vs TM-score, y por candidato
un visor NGL coloreado por pLDDT.
"""
import os, glob, json, base64, csv

METRIC_ORDER = ["rosetta_dg", "rosetta_ddg", "tm_score", "mean_plddt",
                "sap_score", "delta_sap", "n_cysteines", "net_charge", "net_charge_abs"]
METRIC_LABEL = {
    "rosetta_dg": "ΔG (REU)", "rosetta_ddg": "ΔΔG (REU)", "tm_score": "TM-score",
    "mean_plddt": "pLDDT", "sap_score": "SAP", "delta_sap": "ΔSAP",
    "n_cysteines": "Cys", "net_charge": "Carga", "net_charge_abs": "|Carga|",
}


def b64_file(path):
    if not path or not os.path.exists(path):
        return ""
    with open(path, "rb") as fh:
        return base64.b64encode(fh.read()).decode()


def read_ranked(path):
    """ranked.csv -> lista de dicts ordenada por rank (los sin rank al final)."""
    if not os.path.exists(path):
        return []
    rows = []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            rows.append(row)
    def sort_key(r):
        try:
            return (0, int(r["rank"]))
        except (ValueError, KeyError):
            return (1, 9999)
    return sorted(rows, key=sort_key)


def find_pdb(pdb_dir, cid):
    """Busca el PDB de un candidato por id (prefijo del nombre)."""
    if not os.path.isdir(pdb_dir):
        return None
    for p in glob.glob(os.path.join(pdb_dir, "*")):
        bn = os.path.basename(p)
        if bn.startswith(cid) and (bn.endswith(".pdb") or bn.endswith(".cif")):
            return p
    # fallback: contiene el id
    for p in glob.glob(os.path.join(pdb_dir, "*")):
        if cid in os.path.basename(p):
            return p
    return None


def build_round(ranked_csv, pdb_dir):
    """Construye la estructura de datos de una ronda."""
    rows = read_ranked(ranked_csv)
    cands = []
    for r in rows:
        cid = r["id"]
        metrics = {}
        for m in METRIC_ORDER:
            v = r.get(m, "NA")
            try:
                metrics[m] = float(v) if v not in ("NA", "", None) else None
            except (ValueError, TypeError):
                metrics[m] = None
        try:
            score = float(r["score"]) if r.get("score") not in ("", "NA", None) else None
        except (ValueError, TypeError):
            score = None
        pdb_path = find_pdb(pdb_dir, cid)
        cands.append({
            "id": cid,
            "rank": r.get("rank", ""),
            "metrics": metrics,
            "score": score,
            "passed_gates": r.get("passed_gates", ""),
            "in_top_n": r.get("in_top_n", ""),
            "pdb_b64": b64_file(pdb_path) if pdb_path else "",
            "pdb_filename": os.path.basename(pdb_path) if pdb_path else "",
        })
    return cands


DATA = {
    "esm":  build_round("ranked_esm.csv",  "pdbs_esm"),
    "chai": build_round("ranked_chai.csv", "pdbs_chai"),
}

DATA_JSON = json.dumps(DATA)
METRIC_ORDER_JSON = json.dumps(METRIC_ORDER)
METRIC_LABEL_JSON = json.dumps(METRIC_LABEL)

HTML = """\
<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1.0"/>
<title>PhiPsiProt — Cascade Report</title>
<link href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;600&family=Syne:wght@400;600;800&display=swap" rel="stylesheet">
<script src="https://cdn.rawgit.com/arose/ngl/v2.0.0-dev.37/dist/ngl.js"></script>
<script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
<style>
  :root{--bg:#f6f8fc;--surface:#fff;--border:#d9e2ef;--accent:#2563eb;--accent2:#7c3aed;--text:#0f172a;--muted:#475569;}
  *{box-sizing:border-box;margin:0;padding:0;}
  body{background:var(--bg);color:var(--text);font-family:'Syne',sans-serif;font-size:16px;line-height:1.5;}
  .mono{font-family:'JetBrains Mono',monospace;}
  .nav{display:flex;align-items:center;gap:1rem;padding:1.2rem 2rem;border-bottom:1px solid var(--border);background:#fff;position:sticky;top:0;z-index:50;}
  .nav-title{font-size:1.6rem;font-weight:800;background:linear-gradient(90deg,var(--accent),var(--accent2));-webkit-background-clip:text;-webkit-text-fill-color:transparent;}
  .tabs{display:flex;gap:.5rem;padding:1rem 2rem 0;}
  .tab{font-family:'JetBrains Mono',monospace;font-weight:700;padding:.8rem 2rem;border-radius:.8rem .8rem 0 0;border:1px solid var(--border);border-bottom:none;background:#eef2ff;color:var(--muted);cursor:pointer;}
  .tab.active{background:#fff;color:var(--accent);}
  .panel{display:none;padding:1.5rem 2rem 3rem;}
  .panel.active{display:block;}
  .charts{display:grid;grid-template-columns:1fr 1fr;gap:1.5rem;margin-bottom:2rem;}
  @media(max-width:900px){.charts{grid-template-columns:1fr;}}
  .chart-box{background:#fff;border:1px solid var(--border);border-radius:1rem;padding:1rem;}
  .chart-box h3{font-size:.9rem;color:var(--muted);margin-bottom:.5rem;font-family:'JetBrains Mono',monospace;}
  table.rank{width:100%;border-collapse:collapse;background:#fff;border:1px solid var(--border);border-radius:1rem;overflow:hidden;margin-bottom:2rem;font-family:'JetBrains Mono',monospace;font-size:.8rem;}
  table.rank th{background:#f1f5f9;padding:.6rem .8rem;text-align:right;color:var(--muted);border-bottom:1px solid var(--border);}
  table.rank th:first-child,table.rank th:nth-child(2){text-align:left;}
  table.rank td{padding:.5rem .8rem;text-align:right;border-bottom:1px solid #eef2f7;}
  table.rank td:first-child,table.rank td:nth-child(2){text-align:left;}
  table.rank tr.top td{background:rgba(37,99,235,.06);font-weight:600;}
  table.rank tr.gated td{opacity:.45;}
  .badge{display:inline-block;padding:.1rem .5rem;border-radius:9999px;font-size:.65rem;}
  .badge.top{background:rgba(37,99,235,.15);color:var(--accent);}
  .badge.no{background:rgba(239,68,68,.12);color:#dc2626;}
  .card{background:#fff;border:1px solid var(--border);border-radius:1rem;margin-bottom:1.5rem;overflow:hidden;}
  .card-h{display:flex;align-items:center;justify-content:space-between;padding:.8rem 1.2rem;border-bottom:1px solid var(--border);flex-wrap:wrap;gap:.5rem;}
  .card-title{font-weight:800;font-size:1.1rem;}
  .pills{display:flex;gap:.4rem;flex-wrap:wrap;}
  .pill{font-family:'JetBrains Mono',monospace;font-size:.7rem;padding:.2rem .5rem;border:1px solid var(--border);border-radius:.3rem;color:var(--muted);}
  .ngl{width:100%;height:420px;background:#050a14;border-radius:0 0 1rem 1rem;}
  .legend{display:flex;gap:.8rem;padding:.6rem 1.2rem;flex-wrap:wrap;font-size:.75rem;}
  .ldot{width:12px;height:12px;border-radius:3px;display:inline-block;margin-right:.3rem;vertical-align:middle;}
</style>
</head>
<body>
<nav class="nav"><span class="nav-title">PhiPsiProt</span><span class="mono" style="color:var(--muted)">Cascade Report</span></nav>
<div class="tabs" id="tabs"></div>
<div id="panels"></div>
<script>
const DATA = """ + DATA_JSON + """;
const METRIC_ORDER = """ + METRIC_ORDER_JSON + """;
const METRIC_LABEL = """ + METRIC_LABEL_JSON + """;
const ROUNDS = [["esm","Ronda ESMFold"],["chai","Ronda Chai-1"]];
const viewerState = {};

function b64toBlob(b64,mime){const bin=atob(b64),a=new Uint8Array(bin.length);for(let i=0;i<bin.length;i++)a[i]=bin.charCodeAt(i);return new Blob([a],{type:mime});}
function guessExt(fn){return fn&&fn.endsWith('.cif')?'cif':'pdb';}
function fmt(v,d=3){return v==null?'—':(typeof v==='number'?v.toFixed(d):v);}

function plddtScheme(){
  return NGL.ColormakerRegistry.addScheme(function(){
    this.atomColor=function(atom){
      const v=atom.bfactor;
      if(v>=90)return 0x2563eb; if(v>=70)return 0x06b6d4;
      if(v>=50)return 0xeab308; return 0xf97316;
    };
  });
}

function buildTable(cands){
  const cols = METRIC_ORDER.map(m=>`<th>${METRIC_LABEL[m]}</th>`).join('');
  const rows = cands.map(c=>{
    const cls = c.in_top_n==='yes'?'top':(c.passed_gates==='no'?'gated':'');
    const cells = METRIC_ORDER.map(m=>`<td>${fmt(c.metrics[m], m==='n_cysteines'?0:(m==='mean_plddt'?1:3))}</td>`).join('');
    const badge = c.in_top_n==='yes'?'<span class="badge top">TOP</span>':(c.passed_gates==='no'?'<span class="badge no">gated</span>':'');
    return `<tr class="${cls}"><td>${c.rank||'—'}</td><td>${c.id} ${badge}</td>${cells}<td>${fmt(c.score,4)}</td></tr>`;
  }).join('');
  return `<table class="rank"><thead><tr><th>#</th><th>ID</th>${cols}<th>score</th></tr></thead><tbody>${rows}</tbody></table>`;
}

function renderCharts(round, cands){
  // solo candidatos con score (los que pasaron gates)
  const scored = cands.filter(c=>c.score!=null);
  // barras de score
  Plotly.newPlot('bars-'+round,[{
    type:'bar', x:scored.map(c=>c.id), y:scored.map(c=>c.score),
    marker:{color:scored.map(c=>c.in_top_n==='yes'?'#2563eb':'#94a3b8')}
  }],{margin:{t:10,b:90,l:40,r:10},height:300,
      xaxis:{tickangle:-45,tickfont:{size:9,family:'JetBrains Mono'}},
      yaxis:{title:'score'},paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(0,0,0,0)'},
      {responsive:true,displayModeBar:false});
  // dispersion dG vs TM
  Plotly.newPlot('scatter-'+round,[{
    type:'scatter', mode:'markers+text',
    x:scored.map(c=>c.metrics.rosetta_dg), y:scored.map(c=>c.metrics.tm_score),
    text:scored.map(c=>c.id), textposition:'top center', textfont:{size:8},
    marker:{size:11,color:scored.map(c=>c.in_top_n==='yes'?'#2563eb':'#94a3b8')}
  }],{margin:{t:10,b:45,l:55,r:10},height:300,
      xaxis:{title:'ΔG (REU) — menor mejor'},yaxis:{title:'TM-score — mayor mejor'},
      paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(0,0,0,0)'},
      {responsive:true,displayModeBar:false});
}

function initNGL(vid,c){
  if(!c.pdb_b64)return;
  const el=document.getElementById('ngl-'+vid);
  const stage=new NGL.Stage(el,{backgroundColor:'#050a14'});
  new ResizeObserver(()=>stage.handleResize()).observe(el);
  const ext=guessExt(c.pdb_filename), mime=ext==='cif'?'chemical/x-cif':'chemical/x-pdb';
  stage.loadFile(b64toBlob(c.pdb_b64,mime),{ext}).then(comp=>{
    comp.addRepresentation('cartoon',{color:plddtScheme()});
    comp.autoView();
  });
}

const tabsEl=document.getElementById('tabs'), panelsEl=document.getElementById('panels');
ROUNDS.forEach(([key,label],i)=>{
  const cands=DATA[key]||[];
  const tab=document.createElement('button');
  tab.className='tab'+(i===0?' active':''); tab.textContent=label+` (${cands.length})`;
  tab.onclick=()=>{document.querySelectorAll('.tab').forEach((t,j)=>t.classList.toggle('active',j===i));
                   document.querySelectorAll('.panel').forEach(p=>p.classList.toggle('active',p.id==='panel-'+key));};
  tabsEl.appendChild(tab);

  const panel=document.createElement('div');
  panel.id='panel-'+key; panel.className='panel'+(i===0?' active':'');
  panel.innerHTML=`
    <div class="charts">
      <div class="chart-box"><h3>Score por candidato</h3><div id="bars-${key}"></div></div>
      <div class="chart-box"><h3>ΔG vs TM-score</h3><div id="scatter-${key}"></div></div>
    </div>
    ${buildTable(cands)}
    <div class="legend">
      <span><span class="ldot" style="background:#2563eb"></span>pLDDT≥90</span>
      <span><span class="ldot" style="background:#06b6d4"></span>70–90</span>
      <span><span class="ldot" style="background:#eab308"></span>50–70</span>
      <span><span class="ldot" style="background:#f97316"></span>&lt;50</span>
    </div>
    <div id="cards-${key}"></div>`;
  panelsEl.appendChild(panel);

  // tarjetas con visor por candidato (solo los que tienen pdb)
  const cardsEl=panel.querySelector('#cards-'+key);
  cands.forEach(c=>{
    if(!c.pdb_b64)return;
    const vid=(key+'_'+c.id).replace(/[^a-z0-9]/gi,'_');
    viewerState[vid]=c;
    const pills=METRIC_ORDER.map(m=>`<span class="pill">${METRIC_LABEL[m]} ${fmt(c.metrics[m],m==='n_cysteines'?0:2)}</span>`).join('');
    const card=document.createElement('div');
    card.className='card';
    card.innerHTML=`<div class="card-h"><span class="card-title">#${c.rank||'—'} ${c.id}</span><div class="pills">${pills}<span class="pill">score ${fmt(c.score,4)}</span></div></div><div id="ngl-${vid}" class="ngl"></div>`;
    cardsEl.appendChild(card);
    setTimeout(()=>initNGL(vid,c),60);
  });

  setTimeout(()=>renderCharts(key,cands),80);
});
</script>
</body>
</html>
"""

with open("PhiPsiProt_cascade_report.html", "w") as fh:
    fh.write(HTML)
print("Report escrito: PhiPsiProt_cascade_report.html")
