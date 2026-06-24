/*
 * EXTRACT_HADDOCK3_AB_METRICS
 * ──────────────────────────────────────────────────────────────────────────
 * Reads HADDOCK3 capri_ss.tsv output and produces unified metrics TSV
 * for FILTER_RANK_AB.
 */
process EXTRACT_HADDOCK3_AB_METRICS {
    tag "${meta.id}"
    label 'process_single'

    container 'docker.io/python:3.11-slim'

    input:
    tuple val(meta), path(haddock_scores_tsv)

    output:
    tuple val(meta), path("${meta.id}_haddock_metrics.tsv"), emit: metrics
    path "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    pip install pandas --quiet --break-system-packages 2>/dev/null || true

    python3 <<'PYSCRIPT'
import csv, sys

sample_id = "${meta.id}"
tsv_path  = "${haddock_scores_tsv}"

rows = []
with open(tsv_path) as fh:
    reader = csv.DictReader(fh, delimiter="\\t")
    for row in reader:
        rows.append(row)

if not rows:
    print(f"WARNING: empty HADDOCK scores for {sample_id}", file=sys.stderr)
    with open(f"{sample_id}_haddock_metrics.tsv", "w") as fh:
        fh.write("id\\thaddock_score\\tvdw_energy\\telec_energy\\tburied_sasa\\n")
    sys.exit(0)

# Take best scoring model
def get_float(row, *keys):
    for k in keys:
        for col in row:
            if col.lower() == k.lower():
                try:
                    return float(row[col])
                except ValueError:
                    pass
    return 0.0

best = min(rows, key=lambda r: get_float(r, "score", "haddock-score", "haddock_score"))

haddock_score = get_float(best, "score", "haddock-score", "haddock_score")
vdw_energy    = get_float(best, "vdw", "vdw_energy", "evdw")
elec_energy   = get_float(best, "elec", "elec_energy", "eelec")
buried_sasa   = get_float(best, "bsa", "buried_sasa", "buried-sasa")

with open(f"{sample_id}_haddock_metrics.tsv", "w") as fh:
    fh.write("id\\thaddock_score\\tvdw_energy\\telec_energy\\tburied_sasa\\n")
    fh.write(f"{sample_id}\\t{haddock_score:.3f}\\t{vdw_energy:.3f}\\t"
             f"{elec_energy:.3f}\\t{buried_sasa:.3f}\\n")

print(f"{sample_id}: haddock_score={haddock_score:.3f} bsa={buried_sasa:.3f}",
      file=sys.stderr)
PYSCRIPT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf "id\\thaddock_score\\tvdw_energy\\telec_energy\\tburied_sasa\\n" \
        > ${meta.id}_haddock_metrics.tsv
    printf "${meta.id}\\t-25.3\\t-15.2\\t-8.1\\t1250.0\\n" \
        >> ${meta.id}_haddock_metrics.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
