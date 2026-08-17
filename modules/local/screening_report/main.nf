process SCREENING_REPORT {

    tag "screening_report"

    container 'quay.io/phipsiprot/reporting:dev'

    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path final_results
    path plots
    path pymol_images

    output:
    path "screening_report.html"

    script:
    """
    mkdir -p plots
    mkdir -p pymol

    cp ${plots} plots/ || true
    cp ${pymol_images} pymol/ || true

    python3 $projectDir/bin/generate_screening_report.py \
        --results_csv ${final_results} \
        --plots_dir plots \
        --pymol_dir pymol \
        --outdir .
    """
}
