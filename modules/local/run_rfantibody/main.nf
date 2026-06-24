/*
 * RFdiffusion backbone design using rfantibody_local:1.0.0
 *
 * CLI: rfdiffusion -t target.pdb -f framework.pdb -q out.qv -n N
 *                  -l "H1:7,H2:6,H3:5-13" -h "T305,T456"
 *
 * All inputs must be in HLT format (chains H/L/T, CDR REMARKs).
 */
process RUN_RFANTIBODY {
    tag "${meta.id}"
    label 'process_gpu'

    container 'quay.io/rfantibody_local:1.0.0'

    input:
    tuple val(meta),
          path(target_hlt_pdb),
          path(framework_hlt_pdb)
    path  weights_pt             // ← nuevo: fichero .pt de pesos
    val   hotspots
    val   design_loops
    val   num_designs

    output:
    tuple val(meta), path("${meta.id}_rfdiff.qv"), emit: quiver
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def n     = num_designs  ?: params.ab_num_designs  ?: 50
    def loops = design_loops ?: params.ab_design_loops
    def hots  = hotspots     ?: params.ab_hotspot_residues
    """

    TARGET=\$(readlink -f ${target_hlt_pdb})
    FRAMEWORK=\$(readlink -f ${framework_hlt_pdb})
    WEIGHTS=\$(readlink -f ${weights_pt})

    rfdiffusion \\
        -t \$TARGET \\
        -f \$FRAMEWORK \\
        -q ${meta.id}_rfdiff.qv \\
        -n ${n} \\
        -l "${loops}" \\
        -h "${hots}" \\
        -w \$WEIGHTS


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rfantibody: \$(rfdiffusion --version 2>/dev/null || echo "rfantibody_local:1.0.0")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_rfdiff.qv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rfantibody: stub
    END_VERSIONS
    """
}