#!/usr/bin/env bash -C -e -u -o pipefail
if [ -f "/home/lgonlopez/storage/phipsiprot/PhiPsiProt/work/ed/60c1f1e9643c9b55a53dd11a399426/selected_positions.txt" ]; then
    POSITIONS=$(cat /home/lgonlopez/storage/phipsiprot/PhiPsiProt/work/ed/60c1f1e9643c9b55a53dd11a399426/selected_positions.txt)
else
    POSITIONS="/home/lgonlopez/storage/phipsiprot/PhiPsiProt/work/ed/60c1f1e9643c9b55a53dd11a399426/selected_positions.txt"
fi

if [ -z "$POSITIONS" ]; then
    echo "ERROR: No mutable positions selected. Check selection_mode, ligand_resname/interface_chain, and distance_cutoff." >&2
    exit 1
fi

python /home/lgonlopez/storage/phipsiprot/PhiPsiProt/bin/run_pyrosetta_screening.py         --input_pdb test_interface.pdb         --target_chain A         --positions "$POSITIONS"         --output candidates.csv         --fasta_output candidates.fasta
