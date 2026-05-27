process MERGE_SCREENING_RESULTS {

    tag "merge_screening_results"

    container 'quay.io/phipsiprot_ddg_filter:dev'

    publishDir "${params.outdir}/final", mode: 'copy'

    input:
    path candidates_csv
    path selection_summary_csv

    output:
    path "final_screening_results.csv"

    script:
    """
    python $projectDir/bin/merge_screening_results.py \
        --candidates ${candidates_csv} \
        --selection_summary ${selection_summary_csv} \
        --output final_screening_results.csv
    """
}
