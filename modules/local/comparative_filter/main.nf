process COMPARATIVE_FILTER {

    tag "comparative_filter"

    container 'docker.io/biopython/biopython:latest'

    publishDir "${params.outdir}/comparative_filter", mode: 'copy'

    input:
    path input_csv
    path input_pdb

    output:
    path "comparative_filtered.csv"

    script:
    """
    python3 $projectDir/modules/local/comparative_filter/bin/comparative_filter.py \
        --input_csv $input_csv \
        --input_pdb $input_pdb \
        --output_csv comparative_filtered.csv \
        --max_delta_instability ${params.max_delta_instability ?: 0.0}
    """
}
