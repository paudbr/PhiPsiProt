process SCREENING_PLOTS {

    tag "screening_plots"

    container 'quay.io/phipsiprot/screening_plots:dev'

    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path ranked_candidates_csv

    output:
    path "*.png"

    script:
    """

    python3 $projectDir/bin/plot_screening_results.py \
        --input ${ranked_candidates_csv} \
        --outdir .
    """
}
