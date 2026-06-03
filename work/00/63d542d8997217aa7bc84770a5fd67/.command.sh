#!/usr/bin/env bash -C -e -u -o pipefail
python3 /home/lgonlopez/storage/phipsiprot/PhiPsiProt/bin/select_mutable_residues.py         --input_pdb 1BRS.pdb         --target_chain A         --selection_mode all         --distance_cutoff 6.0                  --output selected_positions.txt

echo "selection_mode,target_chain,positions,ligand_resname,interface_chain,distance_cutoff,selected_positions" > selection_summary.csv
echo "all,A,ALL,NONE,NONE,6.0,$(cat selected_positions.txt)" >> selection_summary.csv
