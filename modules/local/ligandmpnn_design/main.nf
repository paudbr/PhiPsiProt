/*
===============================================================================
LIGANDMPNN_DESIGN

Purpose:
Design protein sequences in a functional molecular context using LigandMPNN.

Inputs:
- designed_pdb:
    PDB containing the protein structure and, when applicable, ligand,
    cofactor, metal, DNA, RNA, substrate, or defined pocket context.

Outputs:
- seqs/*.fa
- backbones/*.pdb

Parameters:
- params.ligandmpnn_seed:
    Random seed for LigandMPNN.
- params.ligandmpnn_temperature:
    Sampling temperature for sequence generation.

Notes:
- This module is the core sequence-design step of Mode 3.
- It is used for ligand-aware / function-aware design.
- The input PDB should already contain the functional context to preserve.
===============================================================================
*/

process LIGANDMPNN_DESIGN {

    tag "$designed_pdb"

    executor 'local'

    publishDir "${params.outdir}/designs/ligandmpnn", mode: 'copy'

    input:
    path designed_pdb

    output:
    path "seqs/*.fa"
    path "backbones/*.pdb"

    script:
    """
    mkdir -p ligandmpnn_output

    INPUT_PDB=\$(readlink -f $designed_pdb)

    docker run --rm --gpus all \
        -v \$(dirname \$INPUT_PDB):/input \
        -v \$PWD/ligandmpnn_output:/output \
        docker.io/rosettacommons/ligandmpnn:latest \
        --model_type ligand_mpnn \
        --pdb_path /input/\$(basename \$INPUT_PDB) \
        --out_folder /output \
        --seed ${params.ligandmpnn_seed ?: 111} \
        --temperature ${params.ligandmpnn_temperature ?: 0.1}

    sudo chown -R \$(id -u):\$(id -g) ligandmpnn_output

    cp -r ligandmpnn_output/seqs .
    cp -r ligandmpnn_output/backbones .
    """
}