process DDG_FILTER {

    tag "filter_ddg"

    container 'quay.io/phipsiprot_ddg_filter:dev'

    publishDir "${params.outdir}/filtering", mode: 'copy'

    input:
    path candidates_csv

    output:
    path "filtered_candidates.csv"

    script:
    """
    filter_ddg.py \
        --input ${candidates_csv} \
        --output filtered_candidates.csv \
        --threshold ${params.ddg_threshold}
    """
}
