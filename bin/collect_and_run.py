#!/usr/bin/env python3
"""
Recoge metricas de GNINA y Boltz2, las combina por candidate_id, aplica un
scoring global ponderado (con normalizacion z-score o min-max) y rankea.

Uso:
    collect_and_rank.py \
        --gnina-scores gnina_scores.tsv \
        --boltz-dir boltz_output/ \
        --config scoring_config.yaml \
        --output-table combined_metrics.tsv \
        --output-ranked ranked_candidates.tsv

Notas de orientacion de metricas (todas se llevan internamente a "mas alto = mejor"):
  GNINA
    gnina_cnn_score      [0,1]   higher better  (calidad de pose)
    gnina_cnn_affinity   ~pK     higher better
    gnina_affinity       kcal/mol LOWER better  -> se invierte
  Boltz2
    boltz_confidence     [0,1]   higher better
    boltz_ptm            [0,1]   higher better
    boltz_iptm           [0,1]   higher better
    boltz_ligand_iptm    [0,1]   higher better  (interfaz proteina-ligando)
    boltz_affinity_pred_value     LOWER better  -> se invierte
    boltz_affinity_prob_binary [0,1] higher better
"""
import argparse
import csv
import glob
import json
import math
import os
import sys

try:
    import yaml
except ImportError:
    sys.exit("ERROR: falta pyyaml. Instala con: pip install pyyaml")


# ----------------------------------------------------------------------
# Parseo de GNINA
# ----------------------------------------------------------------------
def parse_gnina_scores(path):
    """Lee gnina_scores.tsv (una fila por replica) y agrega por candidate_id
    haciendo la media de cada metrica sobre las replicas."""
    by_cand = {}
    numeric = ("cnn_score", "cnn_affinity", "affinity")
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cid = row["candidate_id"]
            by_cand.setdefault(cid, {m: [] for m in numeric})
            for m in numeric:
                v = row.get(m, "")
                if v not in ("", None):
                    try:
                        by_cand[cid][m].append(float(v))
                    except ValueError:
                        pass
    out = {}
    for cid, mets in by_cand.items():
        out[cid] = {
            "gnina_cnn_score":    _mean(mets["cnn_score"]),
            "gnina_cnn_affinity": _mean(mets["cnn_affinity"]),
            "gnina_affinity":     _mean(mets["affinity"]),
            "gnina_n_replicates": len(mets["cnn_score"]),
        }
    return out


def _mean(xs):
    return sum(xs) / len(xs) if xs else None


# ----------------------------------------------------------------------
# Parseo de Boltz
# ----------------------------------------------------------------------
def find_boltz_json(boltz_dir, candidate_id, kind):
    """Busca el json de Boltz de forma tolerante a la estructura de carpetas.
    kind = 'confidence' | 'affinity'."""
    patterns = [
        os.path.join(boltz_dir, f"boltz_results_{candidate_id}",
                     "predictions", candidate_id,
                     f"{kind}_*model_0.json" if kind == "confidence"
                     else f"{kind}_*.json"),
        # por si el dir de boltz ya apunta dentro de predictions
        os.path.join(boltz_dir, "**", candidate_id,
                     f"{kind}_*model_0.json" if kind == "confidence"
                     else f"{kind}_*.json"),
        os.path.join(boltz_dir, "**",
                     f"{kind}_{candidate_id}*model_0.json" if kind == "confidence"
                     else f"{kind}_{candidate_id}*.json"),
    ]
    for pat in patterns:
        hits = sorted(glob.glob(pat, recursive=True))
        if hits:
            return hits[0]
    return None


def parse_boltz(boltz_dir, candidate_id):
    res = {
        "boltz_confidence": None, "boltz_ptm": None, "boltz_iptm": None,
        "boltz_ligand_iptm": None, "boltz_complex_plddt": None,
        "boltz_affinity_pred_value": None, "boltz_affinity_prob_binary": None,
    }
    conf = find_boltz_json(boltz_dir, candidate_id, "confidence")
    if conf:
        with open(conf) as fh:
            d = json.load(fh)
        res["boltz_confidence"]   = d.get("confidence_score")
        res["boltz_ptm"]          = d.get("ptm")
        res["boltz_iptm"]         = d.get("iptm")
        res["boltz_ligand_iptm"]  = d.get("ligand_iptm")
        res["boltz_complex_plddt"] = d.get("complex_plddt")

    aff = find_boltz_json(boltz_dir, candidate_id, "affinity")
    if aff:
        with open(aff) as fh:
            d = json.load(fh)
        res["boltz_affinity_pred_value"]  = d.get("affinity_pred_value")
        res["boltz_affinity_prob_binary"] = d.get("affinity_probability_binary")
    return res


# ----------------------------------------------------------------------
# Normalizacion y scoring
# ----------------------------------------------------------------------
def orient(values, direction):
    """Lleva a 'mas alto = mejor'. Si direction == 'lower', invierte signo."""
    if direction == "lower":
        return [(-v if v is not None else None) for v in values]
    return list(values)


def zscore(values):
    present = [v for v in values if v is not None]
    if len(present) < 2:
        # sin varianza utilizable -> todo 0 (neutro)
        return [0.0 if v is not None else None for v in values]
    m = sum(present) / len(present)
    sd = math.sqrt(sum((v - m) ** 2 for v in present) / (len(present) - 1))
    if sd == 0:
        return [0.0 if v is not None else None for v in values]
    return [((v - m) / sd if v is not None else None) for v in values]


def minmax(values):
    present = [v for v in values if v is not None]
    if not present:
        return [None for _ in values]
    lo, hi = min(present), max(present)
    if hi == lo:
        return [0.5 if v is not None else None for v in values]
    return [((v - lo) / (hi - lo) if v is not None else None) for v in values]


def normalize(values, method):
    return zscore(values) if method == "zscore" else minmax(values)


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gnina-scores", required=True)
    ap.add_argument("--boltz-dir", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--output-table", default="combined_metrics.tsv")
    ap.add_argument("--output-ranked", default="ranked_candidates.tsv")
    ap.add_argument("--include-original", action="store_true",
                    help="Incluir el ligando 'original' en el ranking "
                         "(por defecto se mantiene como referencia, fuera del top-N).")
    args = ap.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    method     = cfg.get("normalization", "zscore")
    gates      = cfg.get("gates", {}) or {}
    directions = cfg.get("directions", {}) or {}
    weights    = cfg.get("weights", {}) or {}
    top_n      = int(cfg.get("top_n", 5))

    # 1) recoger metricas crudas
    gnina = parse_gnina_scores(args.gnina_scores)
    candidates = sorted(gnina.keys())
    if not candidates:
        sys.exit("ERROR: no se encontraron candidatos en el TSV de GNINA")

    rows = {}
    for cid in candidates:
        rows[cid] = dict(gnina[cid])
        rows[cid].update(parse_boltz(args.boltz_dir, cid))
        rows[cid]["candidate_id"] = cid

    # 2) gates (filtros duros)
    for cid in candidates:
        passed, reasons = True, []
        for metric, rule in gates.items():
            val = rows[cid].get(metric)
            if val is None:
                passed = False
                reasons.append(f"{metric}=NA")
                continue
            if "min" in rule and val < rule["min"]:
                passed = False
                reasons.append(f"{metric}<{rule['min']}")
            if "max" in rule and val > rule["max"]:
                passed = False
                reasons.append(f"{metric}>{rule['max']}")
        rows[cid]["passed_gates"] = passed
        rows[cid]["gate_fail_reason"] = ";".join(reasons) if reasons else ""

    # 3) normalizar SOLO sobre los que pasan gates (el resto no compite)
    #    pero el 'original' se normaliza junto para tener referencia comparable.
    scoring_set = [c for c in candidates if rows[c]["passed_gates"]]
    if not scoring_set:
        print("WARNING: ningun candidato paso los gates; se rankea sin gates.",
              file=sys.stderr)
        scoring_set = list(candidates)

    # columnas normalizadas + score global
    for metric in weights:
        raw = [rows[c].get(metric) for c in scoring_set]
        oriented = orient(raw, directions.get(metric, "higher"))
        normed = normalize(oriented, method)
        for c, nv in zip(scoring_set, normed):
            rows[c][f"z_{metric}"] = nv

    for cid in candidates:
        if cid not in scoring_set:
            rows[cid]["global_score"] = None
            continue
        score, wsum = 0.0, 0.0
        for metric, w in weights.items():
            nv = rows[cid].get(f"z_{metric}")
            if nv is None:
                continue
            score += w * nv
            wsum += abs(w)
        rows[cid]["global_score"] = round(score, 4) if wsum > 0 else None

    # 4) ranking (excluye 'original' del top-N salvo que se pida)
    def sort_key(cid):
        s = rows[cid].get("global_score")
        return (s is not None, s if s is not None else -1e9)

    rankable = [c for c in scoring_set
                if rows[c].get("global_score") is not None
                and (args.include_original or c != "original")]
    rankable.sort(key=sort_key, reverse=True)

    for i, cid in enumerate(rankable, 1):
        rows[cid]["rank"] = i
        rows[cid]["in_top_n"] = (i <= top_n)

    # 5) escribir tabla completa
    base_cols = [
        "candidate_id", "rank", "in_top_n", "global_score",
        "passed_gates", "gate_fail_reason",
        "gnina_cnn_score", "gnina_cnn_affinity", "gnina_affinity",
        "gnina_n_replicates",
        "boltz_confidence", "boltz_ptm", "boltz_iptm", "boltz_ligand_iptm",
        "boltz_complex_plddt",
        "boltz_affinity_pred_value", "boltz_affinity_prob_binary",
    ]
    z_cols = sorted({k for r in rows.values() for k in r if k.startswith("z_")})
    all_cols = base_cols + z_cols

    def fmt(v):
        if v is None:
            return "NA"
        if isinstance(v, bool):
            return "yes" if v else "no"
        if isinstance(v, float):
            return f"{v:.4f}"
        return str(v)

    ordered = sorted(
        candidates,
        key=lambda c: (rows[c].get("rank", 1e9), c)
    )
    with open(args.output_table, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(all_cols)
        for cid in ordered:
            w.writerow([fmt(rows[cid].get(col)) for col in all_cols])

    # 6) tabla rankeada top-N
    with open(args.output_ranked, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(base_cols)
        for cid in rankable[:top_n]:
            w.writerow([fmt(rows[cid].get(col)) for col in base_cols])

    print(f"[collect_and_rank] {len(candidates)} candidatos, "
          f"{len(rankable)} rankeados, top_n={top_n}", file=sys.stderr)
    print(f"[collect_and_rank] tabla completa -> {args.output_table}", file=sys.stderr)
    print(f"[collect_and_rank] top-N          -> {args.output_ranked}", file=sys.stderr)


if __name__ == "__main__":
    main()