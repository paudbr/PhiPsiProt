process RFANTIBODY_RF2 {

    tag "rfantibody_rf2"

    executor 'local'

    publishDir "${params.outdir}/intermediate", mode: 'copy'

    input:
    path proteinmpnn_qv

    output:
    path "3_rf2.qv"

    script:
    """
    INPUT_QV=\$(readlink -f $proteinmpnn_qv)

    docker run --rm --gpus all \
        --entrypoint /bin/bash \
        -v \$(dirname \$INPUT_QV):/input \
        -v \$PWD:/output \
        phipsiprot/rfantibody:dev \
        -lc "cd /output && rf2 \
            -q /input/\$(basename \$INPUT_QV) \
            --output-quiver 3_rf2.qv \
            -r ${params.antibody_rf2_recycles ?: 10}"

    sudo chown -R \$(id -u):\$(id -g) .
    """
}
