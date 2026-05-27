process RESIDUE_SELECTION {

    tag "$input_pdb"

    container 'quay.io/phipsiprot_pyrosetta:dev'

    publishDir "${params.outdir}/residue_selection", mode: 'copy'

    input:
    path input_pdb
    val target_chain
    val positions
    val selection_mode
    val ligand_resname
    val interface_chain
    val distance_cutoff

    output:
    path "selected_positions.txt"

    script:
    """
    python $projectDir/bin/select_mutable_residues.py \
        --input_pdb $input_pdb \
        --target_chain ${target_chain ?: 'A'} \
        --positions ${positions ?: 'ALL'} \
        --selection_mode ${selection_mode ?: 'manual'} \
        --ligand_resname ${ligand_resname ?: ''} \
        --interface_chain ${interface_chain ?: ''} \
        --distance_cutoff ${distance_cutoff ?: 6.0} \
        --output selected_positions.txt
    """
}
