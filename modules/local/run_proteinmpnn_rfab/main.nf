/*
 * ProteinMPNN — rfantibody's own wrapper, NOT generic ProteinMPNN.
 * Reads HLT REMARK annotations to know which residues are CDRs.
 * Operates on Quiver files natively.
 *
 * CLI: proteinmpnn -q in.qv --output-quiver out.qv -n 4 -t 0.2
 */
process RUN_PROTEINMPNN_RFAB {
    tag "${meta.id}"
    label 'process_medium'

    container 'quay.io/rfantibody_local:1.0.0'

    input:
    tuple val(meta), path(backbones_qv)
    path  weights_pt             // ← nuevo: ProteinMPNN_v48_noise_0.2.pt
    val   seqs_per_struct
    val   temperature

    output:
    tuple val(meta), path("${meta.id}_mpnn.qv"), emit: quiver
    path "versions.yml",                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def n         = seqs_per_struct ?: params.ab_seqs_per_struct  ?: 4
    def temp      = temperature     ?: params.ab_mpnn_temperature ?: 0.2
    def loops_arg = params.ab_mpnn_loops ? "-l \"${params.ab_mpnn_loops}\"" : ''
    """
    proteinmpnn \\
        -q ${backbones_qv} \\
        --output-quiver ${meta.id}_mpnn.qv \\
        -n ${n} \\
        -t ${temp} \\
        -w ${weights_pt} \\
        ${loops_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rfantibody: \$(proteinmpnn --version 2>/dev/null || echo "rfantibody_local:1.0.0")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_mpnn.qv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rfantibody: stub
    END_VERSIONS
    """
}