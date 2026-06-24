/*
 * FILTER_RANK_AB
 * ──────────────────────────────────────────────────────────────────────────
 * Generic filter + ranking module for antibody design pipeline.
 * Reads a YAML scoring config (same format as cascade scoring_config).
 * Applies hard gates, normalizes scores, weights and ranks.
 * Returns top N (absolute) or top pct (percentage) designs.
 *
 * Used after ESMFold2, AF3 and HADDOCK3.
 */
process FILTER_RANK_AB {
    tag "${meta.id}_${stage}"
    label 'process_single'

    container 'docker.io/python:3.11-slim'

    input:
    tuple val(meta),
          val(stage),              // 'esmfold2' | 'af3' | 'haddock3'
          path(metrics_tsv)        // TSV with metrics (one row per design)
    path  scoring_yaml             // YAML scoring config

    output:
    tuple val(meta), val(stage), path("${meta.id}_${stage}_ranked.tsv"),   emit: ranked
    tuple val(meta), val(stage), path("${meta.id}_${stage}_passed.txt"),   emit: passed_ids
    path "versions.yml",                                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    pip install pyyaml --quiet --break-system-packages 2>/dev/null || true

    python3 <<'PYSCRIPT'
import csv, sys, math, yaml
import numpy as np
from pathlib import Path

stage      = "${stage}"
sample_id  = "${meta.id}"
tsv_path   = "${metrics_tsv}"
yaml_path  = "${scoring_yaml}"

# ── Load config ───────────────────────────────────────────────────────────
with open(yaml_path) as fh:
    cfg = yaml.safe_load(fh)

gates      = cfg.get("gates", {})
directions = cfg.get("directions", {})
weights    = cfg.get("weights", {})
norm_type  = cfg.get("normalization", "zscore")
top_n      = cfg.get("top_n", None)
top_pct    = cfg.get("top_pct", None)

# ── Load metrics TSV ──────────────────────────────────────────────────────
rows = []
with open(tsv_path) as fh:
    reader = csv.DictReader(fh, delimiter="\\t")
    for row in reader:
        # Convert numeric columns
        parsed = {"id": row["id"]}
        for k, v in row.items():
            if k == "id":
                continue
            try:
                parsed[k] = float(v)
            except (ValueError, TypeError):
                parsed[k] = v
        rows.append(parsed)

print(f"[FILTER_RANK_AB] {stage}: {len(rows)} designs loaded", file=sys.stderr)

# ── Apply hard gates ──────────────────────────────────────────────────────
passed = []
for row in rows:
    ok = True
    for metric, bounds in gates.items():
        if metric not in row:
            continue
        val = row[metric]
        if "min" in bounds and val < bounds["min"]:
            ok = False
            break
        if "max" in bounds and val > bounds["max"]:
            ok = False
            break
    if ok:
        passed.append(row)

print(f"[FILTER_RANK_AB] {stage}: {len(passed)}/{len(rows)} pass gates", file=sys.stderr)

if not passed:
    print(f"WARNING: no designs passed gates for {sample_id} {stage}", file=sys.stderr)
    # Write empty outputs
    with open(f"{sample_id}_{stage}_ranked.tsv", "w") as fh:
        fh.write("id\\tscore\\trank\\n")
    with open(f"{sample_id}_{stage}_passed.txt", "w") as fh:
        pass
    sys.exit(0)

# ── Normalize and score ───────────────────────────────────────────────────
active_metrics = [m for m in weights if m in passed[0]]

def normalize(values, method):
    arr = np.array(values, dtype=float)
    if method == "zscore":
        mu, sigma = arr.mean(), arr.std()
        if sigma < 1e-9:
            return np.zeros_like(arr)
        return (arr - mu) / sigma
    elif method == "minmax":
        mn, mx = arr.min(), arr.max()
        if mx - mn < 1e-9:
            return np.zeros_like(arr)
        return (arr - mn) / (mx - mn)
    return arr

# Normalize each metric
norm_vals = {}
for m in active_metrics:
    raw = [row[m] for row in passed]
    normed = normalize(raw, norm_type)
    # Invert if lower is better
    if directions.get(m, "higher") == "lower":
        normed = -normed
    norm_vals[m] = normed

# Weighted sum
total_weight = sum(weights[m] for m in active_metrics)
scores = np.zeros(len(passed))
for m in active_metrics:
    w = weights[m] / total_weight
    scores += w * norm_vals[m]

# Add score to rows
for i, row in enumerate(passed):
    row["_score"] = float(scores[i])

# Sort descending (higher score = better)
passed.sort(key=lambda r: r["_score"], reverse=True)

# ── Apply top_n / top_pct cutoff ─────────────────────────────────────────
if top_n is not None:
    final = passed[:int(top_n)]
elif top_pct is not None:
    n = max(1, math.ceil(len(passed) * top_pct / 100.0))
    final = passed[:n]
else:
    final = passed

print(f"[FILTER_RANK_AB] {stage}: {len(final)} designs selected", file=sys.stderr)

# ── Write ranked TSV ──────────────────────────────────────────────────────
fieldnames = ["id", "_score"] + [k for k in passed[0].keys()
                                  if k not in ("id", "_score")]
with open(f"{sample_id}_{stage}_ranked.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["rank"] + fieldnames, delimiter="\\t",
                       extrasaction="ignore")
    w.writeheader()
    for rank, row in enumerate(final, 1):
        w.writerow({"rank": rank, **row})

# ── Write passed IDs (one per line) ──────────────────────────────────────
with open(f"{sample_id}_{stage}_passed.txt", "w") as fh:
    for row in final:
        fh.write(row["id"] + "\\n")

print(f"[FILTER_RANK_AB] Written {sample_id}_{stage}_ranked.tsv "
      f"({len(final)} designs)", file=sys.stderr)
PYSCRIPT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pyyaml: \$(python3 -c "import yaml; print(yaml.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    printf "rank\\tid\\t_score\\n" > ${meta.id}_${stage}_ranked.tsv
    printf "1\\t${meta.id}_stub\\t0.85\\n" >> ${meta.id}_${stage}_ranked.tsv
    printf "${meta.id}_stub\\n" > ${meta.id}_${stage}_passed.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
