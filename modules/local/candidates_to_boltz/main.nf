process CANDIDATES_TO_BOLTZ {
    tag "candidates_to_boltz"
    label 'process_single'

    input:
    path candidates_csv
    val  receptor_sequence   // ← NUEVO (puede ser "" si hay ligando)

    output:
    path "boltz_inputs/*.yaml" , emit: yaml 
    path "boltz_inputs_manifest.tsv", emit: manifest
    path "versions.yml"         , emit: versions

    script:
    def ligand_arg = params.ligand_ccd      ? "--ligand-ccd ${params.ligand_ccd}" :
                     params.ligand_smiles   ? "--ligand-smiles '${params.ligand_smiles}'" :
                     receptor_sequence      ? "--receptor-sequence '${receptor_sequence}'" :
                     ""
    def pocket_arg = params.pocket_residues ? "--pocket-residues '${params.pocket_residues}'" : ""
    def nomsa_arg = params.binder_no_msa ? "--binder_no_msa" : ""
    def affinity_arg = params.with_affinity ? "--with_affinity" : ""
    def task_process = task.process
    """
    candidates_to_boltz.py \\
        --candidates ${candidates_csv} \\
        ${ligand_arg} ${pocket_arg} ${nomsa_arg} ${affinity_arg} \\
        --runs ${params.ligand_filter_runs}

    cat << END_VERSIONS > versions.yml
    "${task_process}":
        python: \$(python3 --version | sed 's/Python //g')
END_VERSIONS
    """
}