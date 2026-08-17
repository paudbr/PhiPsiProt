/*
===============================================================================
BACKBONE_REPORT

Purpose:
Generate an HTML report for Mode 2 backbone stabilization results.

Inputs:
- candidates_biophysical_annotated.csv

Outputs:
- backbone_report.html

Notes:
- Summarizes RFdiffusion backbone designs.
- Summarizes ProteinMPNN sequence designs.
- Reports biophysical annotations.
===============================================================================
*/

process BACKBONE_REPORT {

    tag "backbone_report"

    container 'docker.io/library/python:3.11'

    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path annotated_csv
    path input_pdb
    path designed_pdbs

    output:
    path "backbone_report.html"

    script:
    """
    mkdir -p structures

    cp $input_pdb structures/input.pdb
    cp ${designed_pdbs} structures/

    python3 $projectDir/bin/build_backbone_contig.py \
        --pdb $input_pdb \
        --chain ${params.chain ?: "A"} \
        --design_strategy ${params.design_strategy ?: "local_redesign"} \
        --redesign_region "${params.redesign_region ?: ""}" \
        --output contig.txt

    REPORT_CONTIG=\$(cat contig.txt)

    python3 $projectDir/bin/generate_backbone_report.py \
        --input_csv $annotated_csv \
        --output_html backbone_report.html \
        --structure_dir structures \
        --contig "\$REPORT_CONTIG" \
        --design_strategy "${params.design_strategy ?: 'NA'}" \
        --redesign_region "${params.redesign_region ?: 'NA'}" \
        --chain "${params.chain ?: 'A'}"
    """
}