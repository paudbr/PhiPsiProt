process RANK_SCREENING {

    tag "rank_screening"

    container 'quay.io/phipsiprot_ddg_filter:dev'

    publishDir "${params.outdir}/ranking", mode: 'copy'

    input:
    path filtered_candidates_csv

    output:
    path "ranked_candidates.csv"

    script:
    """
    python $projectDir/bin/rank_screening_candidates.py \
        --input ${filtered_candidates_csv} \
        --output ranked_candidates.csv
    """
}
