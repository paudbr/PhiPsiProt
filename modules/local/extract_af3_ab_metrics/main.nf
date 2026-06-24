/*
 * EXTRACT_AF3_AB_METRICS
 * ──────────────────────────────────────────────────────────────────────────
 * Reads the iptm.tsv and plddt.tsv already emitted by RUN_ALPHAFOLD3
 * and produces a unified metrics TSV for FILTER_RANK_AB.
 */
process EXTRACT_AF3_AB_METRICS {
    tag "${meta.id}"
    label 'process_single'

    container 'docker.io/python:3.11-slim'

    input:
    tuple val(meta),
          path(iptm_tsv),
          path(ptm_tsv),
          path(top_cif)

    output:
    tuple val(meta), path("${meta.id}_af3_metrics.tsv"), emit: metrics
    tuple val(meta), path(top_cif),                      emit: cif
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 <<'PYSCRIPT'
import csv, sys

sample_id = "${meta.id}"
iptm_tsv  = "${iptm_tsv}"
ptm_tsv   = "${ptm_tsv}"

def read_val(tsv_path, metric_name):
    try:
        with open(tsv_path) as fh:
            reader = csv.DictReader(fh, delimiter="\\t")
            for row in reader:
                # ranked format: model, ptm, iptm, aggregate_score
                key = next((k for k in row if k.lower() == metric_name), None)
                if key:
                    return float(row[key])
                # simple format: metric, value
                if row.get("metric", "").lower() == metric_name:
                    return float(row.get("value", 0))
    except Exception as e:
        print(f"WARNING: could not read {tsv_path}: {e}", file=sys.stderr)
    return None

# Try to get max ipTM over seeds (ranked format has multiple rows)
iptm_vals = []
try:
    with open(iptm_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\\t")
        for row in reader:
            for k, v in row.items():
                if k.lower() == "iptm":
                    try:
                        iptm_vals.append(float(v))
                    except ValueError:
                        pass
except Exception as e:
    print(f"WARNING: {e}", file=sys.stderr)

iptm = max(iptm_vals) if iptm_vals else read_val(iptm_tsv, "iptm")
ptm  = read_val(ptm_tsv, "ptm")

# pae_global — not directly available from standard AF3 TSVs
# Use 0.0 as placeholder (will be ignored if not in weights)
pae_global    = 0.0
plddt_global  = 0.0

with open(f"{sample_id}_af3_metrics.tsv", "w") as fh:
    fh.write("id\\tiptm\\tptm\\tpae_global\\tplddt_global\\n")
    fh.write(f"{sample_id}\\t"
             f"{iptm or 0.0:.4f}\\t"
             f"{ptm  or 0.0:.4f}\\t"
             f"{pae_global:.3f}\\t"
             f"{plddt_global:.4f}\\n")

print(f"{sample_id}: iptm={iptm:.4f} ptm={ptm:.4f}", file=sys.stderr)
PYSCRIPT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf "id\\tiptm\\tptm\\tpae_global\\tplddt_global\\n" > ${meta.id}_af3_metrics.tsv
    printf "${meta.id}\\t0.75\\t0.80\\t8.0\\t0.82\\n" >> ${meta.id}_af3_metrics.tsv
    touch stub.cif
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
