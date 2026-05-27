#!/usr/bin/env bash -C -e -u -o pipefail
python /home/lgonlopez/storage/phipsiprot/PhiPsiProt/bin/select_mutable_residues.py         --input_pdb test_interface.pdb         --target_chain A         --positions ALL         --selection_mode interface         --ligand_resname NONE         --interface_chain B         --distance_cutoff 6.0         --output selected_positions.txt

echo "selection_mode,target_chain,positions,ligand_resname,interface_chain,distance_cutoff,selected_positions" > selection_summary.csv
echo "interface,A,ALL,NONE,B,6.0,$(cat selected_positions.txt)" >> selection_summary.csv
