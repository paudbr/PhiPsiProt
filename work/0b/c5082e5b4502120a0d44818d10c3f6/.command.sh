#!/usr/bin/env bash -C -e -u -o pipefail
python /home/lgonlopez/storage/phipsiprot/PhiPsiProt/bin/run_pyrosetta_saturation.py         --input_pdb test.pdb         --target_chain A         --positions 1         --output candidates.csv
