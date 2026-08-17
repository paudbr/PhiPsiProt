process GENERATE_HADDOCK_AIRS {
    tag "generate_airs"
    label 'process_medium'
    input:
    tuple val(candidate_id), path(receptor_pdb), path(ligand_pdb)
    output:
    tuple val(candidate_id),
          path("restraints.tbl"),
          path("active_residues.txt"),
          path("passive_residues.txt"),   emit: airs
    path "versions.yml",                  emit: versions
    script:
    def task_process = task.process
    """
    generate_airs.py \\
        ${receptor_pdb} \\
        ${ligand_pdb} \\
        "${params.binding_site_residues}"

    cat << END_VERSIONS > versions.yml
    "${task_process}":
        python: \$(python3 --version | sed 's/Python //g')
END_VERSIONS
    """
}