/*
===============================================================================
SCREENING_PLOTS

Purpose:
Generate visual summaries of mutational screening results.

Inputs:
- ranked_candidates.csv

Outputs:
- ddg_ranking_barplot.png
- top10_stabilizing_mutations.png
- ddg_heatmap.png

Notes:
- Plots are based on PyRosetta ddG values.
- ddG values are clipped in some plots only for visualization readability.
===============================================================================
*/

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