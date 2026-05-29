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
