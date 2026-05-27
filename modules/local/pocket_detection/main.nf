process POCKET_DETECTION {

    tag "$input_pdb"

    container 'quay.io/phipsiprot_pyrosetta:dev'

    publishDir "${params.outdir}/pocket_detection", mode: 'copy'

    input:
    path input_pdb
    val ligand_resname
    val distance_cutoff

    output:
    path "pocket_positions.txt"

    script:
    """
    python $projectDir/bin/detect_pocket_residues.py \
        --input_pdb $input_pdb \
        --ligand_resname $ligand_resname \
        --distance_cutoff $distance_cutoff \
        --output pocket_positions.txt
    """
}
