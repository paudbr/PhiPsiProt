#!/usr/bin/env bash -C -e -u -o pipefail
if [ -f "/home/lgonlopez/storage/phipsiprot/PhiPsiProt/work/98/7c79e2099a68a912653244056b9006/selected_positions.txt" ]; then
    POSITIONS=$(cat /home/lgonlopez/storage/phipsiprot/PhiPsiProt/work/98/7c79e2099a68a912653244056b9006/selected_positions.txt)
else
    POSITIONS="/home/lgonlopez/storage/phipsiprot/PhiPsiProt/work/98/7c79e2099a68a912653244056b9006/selected_positions.txt"
fi

if [ -z "$POSITIONS" ]; then
    echo "ERROR: No mutable positions selected. Check selection_mode, ligand_resname/interface_chain, pocket_rank and distance_cutoff." >&2
    exit 1
fi

python3 /home/lgonlopez/storage/phipsiprot/PhiPsiProt/bin/run_pyrosetta_screening.py         --input_pdb 3HS4.pdb         --target_chain A         --positions "$POSITIONS"         --output candidates.csv         --fasta_output candidates.fasta
