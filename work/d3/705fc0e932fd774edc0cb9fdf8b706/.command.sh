#!/usr/bin/env bash -C -e -u -o pipefail
if [ -f "/home/lgonlopez/storage/phipsiprot/PhiPsiProt/work/9b/3c580e39c1a43dd463a8f66b43743d/selected_positions.txt" ]; then
    POSITIONS=$(cat /home/lgonlopez/storage/phipsiprot/PhiPsiProt/work/9b/3c580e39c1a43dd463a8f66b43743d/selected_positions.txt)
else
    POSITIONS="/home/lgonlopez/storage/phipsiprot/PhiPsiProt/work/9b/3c580e39c1a43dd463a8f66b43743d/selected_positions.txt"
fi

if [ -z "$POSITIONS" ]; then
    echo "ERROR: No mutable positions selected. Check selection_mode, ligand_resname/interface_chain, and distance_cutoff." >&2
    exit 1
fi

python /home/lgonlopez/storage/phipsiprot/PhiPsiProt/bin/run_pyrosetta_screening.py         --input_pdb test_interface.pdb         --target_chain A         --positions "$POSITIONS"         --output candidates.csv         --fasta_output candidates.fasta
