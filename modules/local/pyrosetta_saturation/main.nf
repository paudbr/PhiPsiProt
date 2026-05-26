process PYROSETTA_SATURATION {
    tag "$input_pdb"

    container 'quay.io/phipsiprot_pyrosetta:dev'

    input:
    path input_pdb
    val target_chain
    val positions

    output:
    path "candidates.csv"

    script:
    """
    python $projectDir/bin/run_pyrosetta_saturation.py \
        --input_pdb $input_pdb \
        --target_chain ${target_chain ?: 'ALL'} \
        --positions ${positions ?: 'ALL'} \
        --output candidates.csv
    """
}
