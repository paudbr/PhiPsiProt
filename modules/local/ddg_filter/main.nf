process DDG_FILTER {

    tag "filter_ddg"

    container 'quay.io/phipsiprot_ddg_filter:dev'

    publishDir "${params.outdir}/filtering", mode: 'copy'

    input:
    path candidates_csv

    output:
    path "filtered_candidates.csv"
    path "filtered_candidates.fasta"
    path "screening_samplesheet.csv"

    script:
    """
    python $projectDir/bin/filter_ddg.py \
        --input ${candidates_csv} \
        --output filtered_candidates.csv \
        --fasta_output filtered_candidates.fasta \
        --samplesheet_output screening_samplesheet.csv \
        --threshold ${params.ddg_threshold}
    """
}
