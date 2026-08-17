/*
===============================================================================
BIOPHYSICAL_ANNOTATION

Purpose:
Annotate screening candidates with sequence-based biophysical descriptors.
Annotate backbone-design candidates with sequence-based biophysical descriptors.


Inputs:
- candidates_ddg_annotated.csv

Outputs:
- candidates_biophysical_annotated.csv

Notes:
- No candidates are removed.
- The module adds developability-related descriptors used later for ranking.
- In Mode 2, annotations are calculated on ProteinMPNN-designed sequences.
===============================================================================
*/

process BIOPHYSICAL_ANNOTATION {

    tag "biophysical_annotation"

    container 'docker.io/biopython/biopython:latest'

    publishDir "${params.outdir}/biophysical_annotation", mode: 'copy'

    input:
    path input_csv

    output:
    path "candidates_biophysical_annotated.csv"

    script:
    """
    python3 $projectDir/bin/annotate_biophysics.py \
        --input_csv $input_csv \
        --output_csv candidates_biophysical_annotated.csv
    """
}