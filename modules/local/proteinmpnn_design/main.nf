/*
===============================================================================
PROTEINMPNN_DESIGN

Purpose:
Design amino-acid sequences compatible with each RFdiffusion-generated backbone.

Inputs:
- designed_pdb

Outputs:
- seqs/*.fa

Parameters:
- params.num_seq_per_target:
    Number of sequences generated per backbone.
- params.sampling_temp:
    Sampling temperature used during sequence generation.
- params.mpnn_model_name:
    ProteinMPNN model variant.

Notes:
- This module does not modify the backbone structure.
- It only designs sequences compatible with the input PDB geometry.
- The generated FASTA files are used downstream for result merging and filtering.
===============================================================================
*/

process PROTEINMPNN_DESIGN {

    tag "$designed_pdb"

    executor 'local'

    publishDir "${params.outdir}/designs/proteinmpnn", mode: 'copy'

    input:
    path designed_pdb

    output:
    path "seqs/*.fa"

    script:
    """
    mkdir -p proteinmpnn_output

    INPUT_PDB=\$(readlink -f $designed_pdb)

    ORIGINAL_PDB="${params.input_pdb}"
    if [ ! -f "\$ORIGINAL_PDB" ]; then
        ORIGINAL_PDB="$projectDir/${params.input_pdb}"
    fi

    python3 $projectDir/bin/build_backbone_contig.py \
        --pdb "\$ORIGINAL_PDB" \
        --chain ${params.chain ?: "A"} \
        --design_strategy ${params.design_strategy ?: "local_redesign"} \
        --redesign_region "${params.redesign_region ?: ""}" \
        --output contig.txt

    CONTIG=\$(cat contig.txt)

    python3 $projectDir/bin/build_mpnn_fixed_positions.py \
        --pdb \$INPUT_PDB \
        --contig "\$CONTIG" \
        --chain ${params.chain ?: "A"} \
        --output fixed_positions.jsonl

    docker run --rm --gpus all \
        -u \$(id -u):\$(id -g) \
        -v \$(dirname \$INPUT_PDB):/input \
        -v \$PWD:/work \
        -v \$PWD/proteinmpnn_output:/output \
        docker.io/rosettacommons/proteinmpnn:latest \
        --pdb_path /input/\$(basename \$INPUT_PDB) \
        --out_folder /output \
        --num_seq_per_target ${params.num_seq_per_target ?: 1} \
        --sampling_temp ${params.sampling_temp ?: 0.1} \
        --model_name ${params.mpnn_model_name ?: 'v_48_020'} \
        --fixed_positions_jsonl /work/fixed_positions.jsonl

    cp -r proteinmpnn_output/seqs .
    """
}