process RFANTIBODY_SCORES {

    tag "rfantibody_scores"

    executor 'local'

    publishDir "${params.outdir}/intermediate", mode: 'copy'

    input:
    path rf2_qv

    output:
    path "rfantibody_scores.tsv"

    script:
    """
    INPUT_QV=\$(readlink -f $rf2_qv)

    docker run --rm --gpus all \
        --entrypoint /bin/bash \
        -v \$(dirname \$INPUT_QV):/input \
        -v \$PWD:/output \
        phipsiprot/rfantibody:dev \
        -lc "cd /output && qvscorefile /input/\$(basename \$INPUT_QV) > rfantibody_scores.tsv"

    sudo chown -R \$(id -u):\$(id -g) .
    """
}
