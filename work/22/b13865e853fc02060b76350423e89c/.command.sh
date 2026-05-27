#!/usr/bin/env bash -C -e -u -o pipefail
if [ -f "/home/lgonlopez/storage/phipsiprot/PhiPsiProt/work/ec/27e5f907fd0451941e88a5814b397a/pocket_positions.txt" ]; then
    POSITIONS=$(cat /home/lgonlopez/storage/phipsiprot/PhiPsiProt/work/ec/27e5f907fd0451941e88a5814b397a/pocket_positions.txt)
else
    POSITIONS="/home/lgonlopez/storage/phipsiprot/PhiPsiProt/work/ec/27e5f907fd0451941e88a5814b397a/pocket_positions.txt"
fi

python /home/lgonlopez/storage/phipsiprot/PhiPsiProt/bin/run_pyrosetta_screening.py         --input_pdb test_ligand.pdb         --target_chain A         --positions $POSITIONS         --output candidates.csv         --fasta_output candidates.fasta
