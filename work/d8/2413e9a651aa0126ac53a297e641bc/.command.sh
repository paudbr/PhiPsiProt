#!/usr/bin/env bash -C -e -u -o pipefail
python3 /home/lgonlopez/storage/phipsiprot/PhiPsiProt/bin/select_mutable_residues.py         --input_pdb 3HS4.pdb         --target_chain A         --selection_mode manual         --distance_cutoff 6.0          --positions 92,215         --output selected_positions.txt

echo "selection_mode,target_chain,positions,ligand_resname,interface_chain,distance_cutoff,selected_positions" > selection_summary.csv
echo "manual,A,92,215,NONE,NONE,6.0,$(cat selected_positions.txt)" >> selection_summary.csv
