/*
===============================================================================
PROTEINMPNN_UNCONDITIONAL

Purpose:
Design complete amino-acid sequences for backbones generated de novo by
RFdiffusion.

This simple test:
- uses no original PDB
- fixes no positions
- designs all chains and all residues present in the input backbone
===============================================================================
*/

process PROTEINMPNN_UNCONDITIONAL {

    tag "${designed_pdb.simpleName}"

    executor 'local'

    publishDir "${params.outdir}/proteinmpnn_unconditional", mode: 'copy'

    input:
    path designed_pdb

    output:
    path "seqs/*.fa", emit: fastas

    script:
    def numSeqs = params.num_seq_per_target ?: 1
    def temperature = params.sampling_temp ?: '0.1'
    def modelName = params.mpnn_model_name ?: 'v_48_020'
    def seed = params.mpnn_seed ?: 0
    def batchSize = params.mpnn_batch_size ?: 1

    """
    mkdir -p proteinmpnn_output

    INPUT_PDB=\$(readlink -f ${designed_pdb})

    docker run --rm --gpus all \
        -u \$(id -u):\$(id -g) \
        -v \$(dirname \$INPUT_PDB):/input \
        -v \$PWD/proteinmpnn_output:/output \
        rosettacommons/proteinmpnn:latest \
        --pdb_path /input/\$(basename \$INPUT_PDB) \
        --out_folder /output \
        --num_seq_per_target ${numSeqs} \
        --sampling_temp "${temperature}" \
        --model_name ${modelName} \
        --seed ${seed} \
        --batch_size ${batchSize}

    cp -r proteinmpnn_output/seqs .
    """
}
