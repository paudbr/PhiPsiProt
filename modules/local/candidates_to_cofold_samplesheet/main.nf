process CANDIDATES_TO_COFOLD_SAMPLESHEET {
    tag "candidates_to_cofold_samplesheet"
    label 'process_single'

    input:
    path boltz_summary      // boltz2_ligand_summary.tsv
    path candidates_csv     // CSV original con binder_sequence
    val  receptor_sequence  // string con la secuencia del receptor

    output:
    path "cofold_samplesheet.csv", emit: samplesheet
    path "fastas/*.fa",            emit: fastas
    path "versions.yml",           emit: versions

    script:
    """
    candidates_to_cofold_samplesheet.py \\
        --summary     ${boltz_summary} \\
        --candidates  ${candidates_csv} \\
        --receptor    "${receptor_sequence}" \\
        --outdir      fastas

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}