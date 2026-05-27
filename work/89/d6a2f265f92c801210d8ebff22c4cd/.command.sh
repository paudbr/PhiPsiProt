#!/usr/bin/env bash -C -e -u -o pipefail
python /home/lgonlopez/storage/phipsiprot/PhiPsiProt/bin/run_pyrosetta_screening.py         --input_pdb test_multi.pdb         --target_chain A         --positions 2         --output candidates.csv         --fasta_output candidates.fasta
