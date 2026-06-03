process SCREENING_REPORT {

    tag "screening_report"

    container 'quay.io/phipsiprot/reporting:dev'

    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path final_results
    path plots

    output:
    path "screening_report.html"

    script:
    """
    mkdir -p plots
    cp ${plots} plots/ || true

    python3 $projectDir/bin/generate_screening_report.py \
        --results_csv ${final_results} \
        --plots_dir plots \
        --outdir .
    """
}
