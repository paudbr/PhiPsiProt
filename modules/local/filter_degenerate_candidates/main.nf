process FILTER_DEGENERATE_CANDIDATES {
    tag "filter_degenerate"
    label 'process_single'

    input:
    path candidates_csv

    output:
    path "candidates_clean.csv",    emit: clean
    path "candidates_rejected.csv", optional: true, emit: rejected
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def max_run     = params.binder_max_homopolymer ?: 5
    def min_entropy = params.binder_min_entropy      ?: 2.5
    def max_single  = params.binder_max_single_aa    ?: 0.40
    def col_seq     = params.seq_col ?: "final_sequence"
    """
    filter_degenerate.py \\
        --input ${candidates_csv} \\
        --output candidates_clean.csv \\
        --rejected candidates_rejected.csv \\
        --max-run ${max_run} \\
        --min-entropy ${min_entropy} \\
        --max-single-aa ${max_single} \\
        --seq-col ${col_seq}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
	END_VERSIONS
    """
}