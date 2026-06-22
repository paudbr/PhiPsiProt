process CANDIDATES_CSV_TO_FASTA {
    tag "csv_to_fasta"
    input:
    path candidates_csv
    output:
    path "*.fasta",      emit: fastas
    path "versions.yml", emit: versions
    script:
    def task_process = task.process
    """
    candidates_to_csv_fasta.py ${candidates_csv}

    cat << END_VERSIONS > versions.yml
    "${task_process}":
        python: \$(python3 --version | sed 's/Python //g')
END_VERSIONS
    """
}