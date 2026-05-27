process PYROSETTA_SCREENING {

    tag "$input_pdb"

    container 'quay.io/phipsiprot_pyrosetta:dev'

    publishDir "${params.outdir}/screening", mode: 'copy'

    input:
    path input_pdb
    val target_chain
    val positions

    output:
    path "candidates.csv"
    path "candidates.fasta"

    script:
    """
    python $projectDir/bin/run_pyrosetta_screening.py \
        --input_pdb $input_pdb \
        --target_chain ${target_chain ?: 'ALL'} \
        --positions ${positions ?: 'ALL'} \
        --output candidates.csv \
        --fasta_output candidates.fasta
    """
}
