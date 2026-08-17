/*
===============================================================================
DDG_ANNOTATION

Purpose:
Annotate PyRosetta screening candidates with ddG threshold information.

Inputs:
- candidates.csv

Outputs:
- candidates_ddg_annotated.csv

Notes:
- No candidates are removed.
- ddg_pass indicates whether a candidate satisfies the configured threshold.
- More negative ddG values are interpreted as more stabilizing.
===============================================================================
*/

process DDG_ANNOTATION {

    tag "ddg_annotation"

    container 'quay.io/phipsiprot_ddg_filter:dev'

    publishDir "${params.outdir}/ddg_annotation", mode: 'copy'

    input:
    path candidates_csv

    output:
    path "candidates_ddg_annotated.csv"

    script:
    def threshold = params.ddg_threshold ?: 2.0

    """
    python3 $projectDir/bin/annotate_ddg.py \
        --input ${candidates_csv} \
        --output candidates_ddg_annotated.csv \
        --threshold ${threshold}
    """
}