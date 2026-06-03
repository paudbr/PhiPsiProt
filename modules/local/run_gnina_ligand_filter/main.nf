process RUN_GNINA_LIGAND_FILTER {
    when:
    params.docking_tool == 'gnina'
    tag "gnina_ligand_filter"
    label 'process_medium'
    label 'process_gpu'

    container "quay.io/nf-core/proteinfold_gnina:1.3.2"

    input:
    path candidates_csv
    path reference_pdb
    path ligand_file

    output:
    path "gnina_scores.tsv"  , emit: scores
    path "gnina_summary.tsv" , emit: summary
    path "*.gnina.log"       , emit: logs
    path "*.sdf"             , optional: true, emit: poses
    path "versions.yml"      , emit: versions
    path "gnina_score_distributions.pdf"     , emit: plots  

    script:
    def ligand_arg = ligand_file.name == "NO_FILE" ? "" : "--ligand-file ${ligand_file}"
    """
    gnina_ligand_filter.py \\
        --candidates ${candidates_csv} \\
        --reference-pdb ${reference_pdb} \\
        ${ligand_arg} \\
        --ligand-resname ${params.ligand_resname} \\
        --pocket-residues "${params.pocket_residues}" \\
        --box-size ${params.gnina_box_size} \\
        --runs ${params.ligand_filter_runs} \\
        --gnina-exhaustiveness ${params.gnina_exhaustiveness} \\
        --gnina-num-modes ${params.gnina_num_modes} \\
        --gnina-seed ${params.gnina_seed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gnina: \$(gnina --version 2>&1 | head -n 1 || echo "unknown")
    END_VERSIONS
    """
}
