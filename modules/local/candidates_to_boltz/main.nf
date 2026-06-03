process CANDIDATES_TO_BOLTZ {
    tag "candidates_to_boltz"
    label 'process_single'

    input:
    path candidates_csv
    val  receptor_sequence   // ← NUEVO (puede ser "" si hay ligando)

    output:
    path "boltz_inputs/*.fasta" , emit: fasta
    path "boltz_inputs_manifest.tsv", emit: manifest
    path "versions.yml"         , emit: versions

    script:
    def ligand_arg = params.ligand_ccd      ? "--ligand-ccd ${params.ligand_ccd}" :
                     params.ligand_smiles   ? "--ligand-smiles '${params.ligand_smiles}'" :
                     receptor_sequence      ? "--receptor-sequence '${receptor_sequence}'" :
                     ""
    """
    candidates_to_boltz.py \\
        --candidates ${candidates_csv} \\
        ${ligand_arg} \\
        --runs ${params.ligand_filter_runs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}