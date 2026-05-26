#!/usr/bin/env bash -C -e -u -o pipefail
python /home/lgonlopez/storage/phipsiprot/PhiPsiProt/bin/run_pyrosetta_saturation.py         --input_pdb test.pdb         --output candidates.csv
