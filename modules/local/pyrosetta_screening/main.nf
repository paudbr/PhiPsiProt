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
    if [ -f "${positions}" ]; then
        POSITIONS=\$(cat ${positions})
    else
        POSITIONS="${positions ?: 'ALL'}"
    fi

    if [ -z "\$POSITIONS" ]; then
        echo "ERROR: No mutable positions selected. Check selection_mode, ligand_resname/interface_chain, and distance_cutoff." >&2
        exit 1
    fi

    python $projectDir/bin/run_pyrosetta_screening.py \
        --input_pdb $input_pdb \
        --target_chain ${target_chain ?: 'ALL'} \
        --positions "\$POSITIONS" \
        --output candidates.csv \
        --fasta_output candidates.fasta
    """
    
}
