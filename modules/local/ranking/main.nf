/*
===============================================================================
RANK_SCREENING

Purpose:
Rank annotated screening candidates using PyRosetta ddG and biophysical
developability descriptors.

Inputs:
- candidates_biophysical_annotated.csv

Outputs:
- ranked_candidates.csv

Notes:
- No candidates are removed.
- Lower final_score is better.
===============================================================================
*/

process RANK_SCREENING {

    tag "rank_screening"

    container 'quay.io/phipsiprot_ddg_filter:dev'

    publishDir "${params.outdir}/ranking", mode: 'copy'

    input:
    path annotated_candidates_csv

    output:
    path "ranked_candidates.csv"

    script:
    """
    python3 $projectDir/bin/rank_screening_candidates.py \
        --input ${annotated_candidates_csv} \
        --output ranked_candidates.csv
    """
}