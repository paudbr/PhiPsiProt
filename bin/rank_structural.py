#!/usr/bin/env python3
"""
rank_structural.py — combina los TSV de métricas de todos los candidatos,
aplica un scoring ponderado configurable por YAML y produce un ranking.

Entradas:
  - varios  *_metrics.tsv  (formato largo: id  metric  value)
  - scoring_config.yaml

YAML:
  normalization: zscore | minmax
  gates:                       # filtros duros (solo si --apply-gates)
    n_cysteines:   { max: 2 }
    net_charge:    { min: -10, max: 10 }
  directions:                  # 'lower' invierte (menor es mejor)
    rosetta_dg:    lower
    sap_score:     lower
    net_charge_abs: lower
  weights:                     # solo las métricas que ponderan
    rosetta_dg:   0.35
    tm_score:     0.25
    mean_plddt:   0.20
    sap_score:    0.20
  top_n: 5

Flags:
  --apply-gates / --no-apply-gates   (post-ESM aplica; post-Chai no)
  --select-top-n / --no-select-top-n (post-ESM recorta; post-Chai solo ordena)

Salida: ranked.csv con columnas:
  rank, id, <cada métrica>, score, passed_gates, in_top_n
"""
import argparse
import csv
import glob
import sys
import math


def load_yaml(path):
    """Carga el YAML de scoring. Usa PyYAML si está; si no, un mini-parser
    en stdlib que entiende EXACTAMENTE el formato de scoring_config:
      - claves de nivel 0:  normalization, top_n, y bloques gates/directions/weights
      - bloques anidados con 2 espacios de indentación
      - gates con dict inline {min: x, max: y}
      - directions/weights con valor escalar
    No es un parser YAML general, pero cubre nuestro fichero sin dependencias.
    """
    try:
        import yaml
        with open(path) as fh:
            return yaml.safe_load(fh)
    except ImportError:
        return _parse_scoring_yaml_stdlib(path)


def _coerce(v):
    v = v.strip()
    if v == "" or v.lower() in ("null", "none", "~"):
        return None
    low = v.lower()
    if low in ("true", "yes"):
        return True
    if low in ("false", "no"):
        return False
    try:
        if "." in v or "e" in low:
            return float(v)
        return int(v)
    except ValueError:
        return v.strip('"').strip("'")


def _parse_inline_dict(s):
    # "{ min: 0.0, max: 10 }" -> {"min":0.0, "max":10}
    s = s.strip().lstrip("{").rstrip("}").strip()
    out = {}
    if not s:
        return out
    for part in s.split(","):
        if ":" not in part:
            continue
        k, v = part.split(":", 1)
        out[k.strip()] = _coerce(v)
    return out


def _parse_scoring_yaml_stdlib(path):
    cfg = {}
    current_block = None       # nombre del bloque de nivel 1 (gates/directions/weights)
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            # quita comentarios (fuera de comillas; nuestro fichero no usa # en valores)
            if "#" in line:
                line = line[:line.index("#")]
            if not line.strip():
                continue

            indent = len(line) - len(line.lstrip(" "))
            stripped = line.strip()

            if indent == 0:
                # clave de nivel 0
                if stripped.endswith(":"):
                    # abre un bloque (gates:, directions:, weights:)
                    key = stripped[:-1].strip()
                    cfg[key] = {}
                    current_block = key
                elif ":" in stripped:
                    key, val = stripped.split(":", 1)
                    cfg[key.strip()] = _coerce(val)
                    current_block = None
            else:
                # dentro de un bloque
                if current_block is None or ":" not in stripped:
                    continue
                key, val = stripped.split(":", 1)
                key = key.strip()
                val = val.strip()
                if val.startswith("{"):
                    cfg[current_block][key] = _parse_inline_dict(val)
                else:
                    cfg[current_block][key] = _coerce(val)
    return cfg


def load_metrics(tsv_glob):
    """Lee todos los *_metrics.tsv (largo) -> {id: {metric: value|None}}"""
    data = {}
    for f in sorted(glob.glob(tsv_glob)):
        with open(f) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                cid = row["id"]
                metric = row["metric"]
                raw = row["value"]
                try:
                    val = float(raw) if raw not in ("NA", "", "nan", "None") else None
                except ValueError:
                    val = None
                data.setdefault(cid, {})[metric] = val
    return data


def normalize(values, method):
    """values: dict id->float (sin None). Devuelve dict id->z/minmax."""
    vals = list(values.values())
    if not vals:
        return {}
    if method == "minmax":
        lo, hi = min(vals), max(vals)
        rng = hi - lo
        if rng < 1e-12:
            return {k: 0.5 for k in values}
        return {k: (v - lo) / rng for k, v in values.items()}
    else:  # zscore
        mean = sum(vals) / len(vals)
        var = sum((v - mean) ** 2 for v in vals) / len(vals)
        sd = math.sqrt(var)
        if sd < 1e-12:
            return {k: 0.0 for k in values}
        return {k: (v - mean) / sd for k, v in values.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--metrics-glob", default="*_metrics.tsv")
    ap.add_argument("--out", default="ranked.csv")
    ap.add_argument("--apply-gates", dest="apply_gates", action="store_true")
    ap.add_argument("--no-apply-gates", dest="apply_gates", action="store_false")
    ap.add_argument("--select-top-n", dest="select_top_n", action="store_true")
    ap.add_argument("--no-select-top-n", dest="select_top_n", action="store_false")
    ap.set_defaults(apply_gates=True, select_top_n=True)
    args = ap.parse_args()

    cfg = load_yaml(args.config)
    norm_method = cfg.get("normalization", "zscore")
    gates       = cfg.get("gates", {}) or {}
    directions  = cfg.get("directions", {}) or {}
    weights     = cfg.get("weights", {}) or {}
    top_n       = int(cfg.get("top_n", 5))

    data = load_metrics(args.metrics_glob)
    if not data:
        raise SystemExit("No se leyó ninguna métrica.")

    ids = sorted(data.keys())
    all_metrics = sorted({m for d in data.values() for m in d})

    # ── gates (filtros duros) ───────────────────────────────────────────────
    passed = {cid: True for cid in ids}
    if args.apply_gates:
        for cid in ids:
            for metric, lims in gates.items():
                v = data[cid].get(metric)
                if v is None:
                    passed[cid] = False  # sin valor no pasa el gate (conservador)
                    continue
                if "min" in lims and v < lims["min"]:
                    passed[cid] = False
                if "max" in lims and v > lims["max"]:
                    passed[cid] = False

    # ── normalización por métrica ponderada (solo sobre los que pasan gates) ─
    ranking_ids = [c for c in ids if passed[c]] if args.apply_gates else ids[:]
    if not ranking_ids:
        sys.stderr.write("[WARN] ningún candidato pasó los gates.\n")

    norm_cache = {}
    for metric in weights:
        vals = {c: data[c].get(metric) for c in ranking_ids
                if data[c].get(metric) is not None}
        n = normalize(vals, norm_method)
        # dirección: 'lower' -> invertir signo (menor mejor => mayor score)
        if str(directions.get(metric, "higher")).lower() == "lower":
            n = {k: -v for k, v in n.items()}
        norm_cache[metric] = n

    wsum = sum(weights.values()) or 1.0
    score = {}
    for cid in ranking_ids:
        s = 0.0
        for metric, w in weights.items():
            s += w * norm_cache.get(metric, {}).get(cid, 0.0)
        score[cid] = s / wsum

    order = sorted(ranking_ids, key=lambda c: score.get(c, float("-inf")),
                   reverse=True)

    # ── escribir ranked.csv ─────────────────────────────────────────────────
    with open(args.out, "w", newline="") as out:
        w = csv.writer(out)
        header = ["rank", "id"] + all_metrics + ["score", "passed_gates", "in_top_n"]
        w.writerow(header)

        # primero los rankeados (los que entran al score), luego los excluidos por gate
        rank = 0
        written = set()
        for cid in order:
            rank += 1
            in_top = "yes" if (not args.select_top_n or rank <= top_n) else "no"
            row = [rank, cid]
            row += [_fmt(data[cid].get(m)) for m in all_metrics]
            row += [f"{score[cid]:.4f}", "yes", in_top]
            w.writerow(row)
            written.add(cid)

        # candidatos que no pasaron gates (sin score), al final
        for cid in ids:
            if cid in written:
                continue
            row = ["", cid]
            row += [_fmt(data[cid].get(m)) for m in all_metrics]
            row += ["", "no", "no"]
            w.writerow(row)

    n_top = sum(1 for c in order if (not args.select_top_n or order.index(c) < top_n))
    sys.stderr.write(
        f"[OK] {len(ids)} candidatos, {len(ranking_ids)} pasaron gates, "
        f"top_n={'all' if not args.select_top_n else top_n}.\n"
    )


def _fmt(v):
    if v is None:
        return "NA"
    return f"{v:.6g}"


if __name__ == "__main__":
    main()