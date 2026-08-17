/*
===============================================================================
MERGE_SCREENING_RESULTS

Purpose:
Merge ranked screening candidates with residue selection metadata.

Inputs:
- ranked_candidates.csv
- selection_summary.csv

Outputs:
- final_screening_results.csv
===============================================================================
*/

process MERGE_SCREENING_RESULTS {

    tag "merge_screening_results"

    container 'quay.io/phipsiprot_ddg_filter:dev'

    publishDir "${params.outdir}/final", mode: 'copy'

    input:
    path ranked_candidates_csv
    path selection_summary_csv

    output:
    path "final_screening_results.csv"

    script:
    """
    python3 $projectDir/bin/merge_screening_results.py \
        --candidates ${ranked_candidates_csv} \
        --selection_summary ${selection_summary_csv} \
        --output final_screening_results.csv
    """
}