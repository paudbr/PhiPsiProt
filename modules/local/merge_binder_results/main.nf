process MERGE_BINDER_RESULTS {

    tag "merge_binder_results"

    publishDir "${params.outdir}/final", mode: 'copy'

    input:
    path candidates_csv
    path mpnn_fastas

    output:
    path "final_binder_candidates.csv"

    script:
    """
    python $projectDir/bin/merge_binder_results.py \
        --candidates_csv $candidates_csv \
        --mpnn_fastas ${mpnn_fastas} \
        --output_csv final_binder_candidates.csv
    """
}
