#!/usr/bin/env bash -C -e -u -o pipefail
python /home/lgonlopez/storage/phipsiprot/PhiPsiProt/bin/detect_pocket_residues.py         --input_pdb test_ligand.pdb         --ligand_resname ATP         --distance_cutoff 6.0         --output pocket_positions.txt
