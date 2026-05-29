process BIOPHYSICAL_FILTER {

    tag "biophysical_filter"

    container 'docker.io/biopython/biopython:latest'

    publishDir "${params.outdir}/intermediate", mode: 'copy'

    input:
    path input_csv

    output:
    path "biophysical_filtered.csv"

    script:
    """
    python3 $projectDir/modules/local/biophysical_filter/bin/biophysical_filter.py \
        --input_csv $input_csv \
        --output_csv biophysical_filtered.csv \
        --filters "${params.biophysical_filters}"
    """
}
