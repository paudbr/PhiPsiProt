process BACKBONE_SUMMARY {

    tag "$designed_pdb"

    container 'quay.io/phipsiprot_backbone:dev'

    publishDir "${params.outdir}/backbone", mode: 'copy'

    input:
    path designed_pdb

    output:
    path "backbone_candidates.csv"
    path "backbone_candidates.fasta"
    path "backbone_samplesheet.csv"

    script:
    """
    python $projectDir/bin/summarize_backbone_design.py \
        --pdb $designed_pdb \
        --csv backbone_candidates.csv \
        --fasta backbone_candidates.fasta \
        --samplesheet backbone_samplesheet.csv
    """
}
