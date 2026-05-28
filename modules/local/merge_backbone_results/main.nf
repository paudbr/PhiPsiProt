process MERGE_BACKBONE_RESULTS {

    publishDir "${params.outdir}/final", mode: 'copy'

    input:
    path candidates_csv
    path mpnn_fasta

    output:
    path "final_backbone_results.csv"

    script:
    """
    python $projectDir/bin/merge_backbone_results.py \
        --candidates_csv $candidates_csv \
        --mpnn_fasta $mpnn_fasta \
        --output_csv final_backbone_results.csv
    """
}
