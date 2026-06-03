process RFANTIBODY_PROTEINMPNN {

    tag "rfantibody_proteinmpnn"

    executor 'local'

    publishDir "${params.outdir}/intermediate", mode: 'copy'

    input:
    path rfdiffusion_qv

    output:
    path "2_proteinmpnn.qv"

    script:
    """
    INPUT_QV=\$(readlink -f $rfdiffusion_qv)

    docker run --rm --gpus all \
        --entrypoint /bin/bash \
        -v \$(dirname \$INPUT_QV):/input \
        -v \$PWD:/output \
        phipsiprot/rfantibody:dev \
        -lc "cd /output && proteinmpnn \
            -q /input/\$(basename \$INPUT_QV) \
            --output-quiver 2_proteinmpnn.qv \
            -n ${params.antibody_num_seqs ?: 4} \
            -t ${params.antibody_sampling_temp ?: 0.2}"

    sudo chown -R \$(id -u):\$(id -g) .
    """
}
