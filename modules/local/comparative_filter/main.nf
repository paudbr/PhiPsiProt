process COMPARATIVE_FILTER {

    tag "comparative_filter"

    container 'docker.io/biopython/biopython:latest'

    publishDir "${params.outdir}/final", mode: 'copy'

    input:
    path input_csv
    path input_pdb

    output:
    path "final_backbone_candidates.csv"

    script:
    """
    python3 $projectDir/modules/local/comparative_filter/bin/comparative_filter.py \
      --input_csv biophysical_filtered.csv \
      --input_pdb $input_pdb \
      --output_csv final_backbone_candidates.csv \
      --max_delta_instability ${params.max_delta_instability}
    """
}
