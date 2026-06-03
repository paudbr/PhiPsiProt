process RFANTIBODY_DESIGN {

    tag "rfantibody_design"

    executor 'local'

    publishDir "${params.outdir}/intermediate", mode: 'copy'

    input:
    path target_pdb
    path framework_pdb

    output:
    path "1_rfdiffusion.qv"

    script:
    """
    TARGET_PDB=\$(readlink -f $target_pdb)
    FRAMEWORK_PDB=\$(readlink -f $framework_pdb)

    docker run --rm --gpus all \
        --entrypoint /bin/bash \
        -v \$(dirname \$TARGET_PDB):/target \
        -v \$(dirname \$FRAMEWORK_PDB):/framework \
        -v \$PWD:/output \
        phipsiprot/rfantibody:dev \
        -lc "cd /output && rfdiffusion \
            -t /target/\$(basename \$TARGET_PDB) \
            -f /framework/\$(basename \$FRAMEWORK_PDB) \
            -q 1_rfdiffusion.qv \
            -n ${params.num_designs ?: 1} \
            -l '${params.antibody_loops}' \
            -h '${params.antibody_hotspots}'"

    sudo chown -R \$(id -u):\$(id -g) .
    """
}
