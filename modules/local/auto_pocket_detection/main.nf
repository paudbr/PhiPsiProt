/*
Module: AUTO_POCKET_DETECTION

Purpose:
Detect binding pockets automatically using fpocket.

This module is used when no ligand is available. fpocket predicts pockets
from the protein structure, and the selected pocket rank is converted into
a list of residue positions for Mode 1 screening.

Outputs:
- selected_positions.txt
- selection_summary.csv
*/

process AUTO_POCKET_DETECTION {

    tag "$input_pdb"

    container 'quay.io/phipsiprot/fpocket:dev'

    publishDir "${params.outdir}/auto_pocket", mode: 'copy'

    input:
    path input_pdb
    val target_chain
    val pocket_rank

    output:
    path "selected_positions.txt"
    path "selection_summary.csv"

    script:
    """
    fpocket -f $input_pdb

    python3 $projectDir/bin/fpocket_to_positions.py \
        --pocket_pdb *_out/pockets/pocket${pocket_rank}_atm.pdb \
        --target_chain ${target_chain} \
        --output selected_positions.txt

    echo "selection_mode,target_chain,pocket_rank,selected_positions" > selection_summary.csv
    echo "auto_pocket,${target_chain},${pocket_rank},\$(cat selected_positions.txt)" >> selection_summary.csv
    """
}