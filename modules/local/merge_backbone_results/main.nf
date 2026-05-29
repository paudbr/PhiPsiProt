process MERGE_BACKBONE_RESULTS {

    tag "merge_backbone_results"

    publishDir "${params.outdir}/final", mode: 'copy'

    input:
    path candidates_csv
    path mpnn_fastas

    output:
    path "final_backbone_results.csv"

    script:
    """
    python $projectDir/bin/merge_backbone_results.py \
        --candidates_csv $candidates_csv \
        --mpnn_fastas ${mpnn_fastas} \
        --output_csv final_backbone_results.csv
    """
}