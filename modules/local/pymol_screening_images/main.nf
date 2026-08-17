/*
===============================================================================
PYMOL_SCREENING_IMAGES

Purpose:
Generate PyMOL structural images for the screening report.

Inputs:
- input_pdb
- final_screening_results.csv

Outputs:
- target_structure.png
- top_candidate_mutation.png
- PyMOL .pml scripts
===============================================================================
*/

process PYMOL_SCREENING_IMAGES {

    tag "pymol_screening_images"

    container 'quay.io/phipsiprot_pymol_visualization:dev'

    publishDir "${params.outdir}/pymol", mode: 'copy'

    input:
    path input_pdb
    path final_results_csv

    output:
    path "*.png", emit: images
    path "*.pml", emit: scripts

    script:
    """
    python3 $projectDir/bin/generate_screening_pymol_images.py \
        --input_pdb ${input_pdb} \
        --results_csv ${final_results_csv} \
        --outdir .
    """
}