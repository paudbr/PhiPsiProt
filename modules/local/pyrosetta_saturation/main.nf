process PYROSETTA_SATURATION {
    tag "$input_pdb"

    container 'phipsiprot_pyrosetta:dev'

    input:
    path input_pdb

    output:
    path "candidates.csv"

    script:
    """
    python $projectDir/bin/run_pyrosetta_saturation.py \
        --input_pdb $input_pdb \
        --output candidates.csv
    """
}
