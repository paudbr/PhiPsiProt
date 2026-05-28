process PROTEINMPNN_DESIGN {

    tag "$designed_pdb"

    executor 'local'

    publishDir "${params.outdir}/proteinmpnn", mode: 'copy'

    input:
    path designed_pdb

    output:
    path "seqs/*.fa"

    script:
    """
    mkdir -p proteinmpnn_output

    INPUT_PDB=\$(readlink -f $designed_pdb)

    docker run --rm --gpus all \
        -u \$(id -u):\$(id -g) \
        -v \$(dirname \$INPUT_PDB):/input \
        -v \$PWD/proteinmpnn_output:/output \
        docker.io/rosettacommons/proteinmpnn:latest \
        --pdb_path /input/\$(basename \$INPUT_PDB) \
        --out_folder /output \
        --num_seq_per_target ${params.num_seq_per_target ?: 1} \
        --sampling_temp ${params.sampling_temp ?: 0.1} \
        --model_name ${params.mpnn_model_name ?: 'v_48_020'}

    cp -r proteinmpnn_output/seqs .
    """
}
