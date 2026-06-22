process COLLECT_AND_RANK {
    tag "collect_and_rank"
    label 'process_single'
    conda "conda-forge::python=3.11 conda-forge::pyyaml=6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    path gnina_scores
    path boltz_dirs           // collect() de los boltz_results_*
    path scoring_config
    path candidates_csv

    output:
    path "ranked_candidates.tsv" , emit: ranked
    path "combined_metrics.tsv"  , emit: table
    path "versions.yml"          , emit: versions

    script:
    """
    collect_and_rank.py \\
        --gnina-scores ${gnina_scores} \\
        --boltz-dir . \\
        --config ${scoring_config} \\
        --candidates-csv ${candidates_csv} \\
        --output-table combined_metrics.tsv \\
        --output-ranked ranked_candidates.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
END_VERSIONS
    """
}