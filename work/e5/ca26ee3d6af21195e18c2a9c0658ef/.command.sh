#!/usr/bin/env bash -C -e -u -o pipefail
echo "id,fasta,pdb,mode" > candidates.csv
echo "dummy_antibody_design,dummy.fasta,dummy.pdb,antibody_design" >> candidates.csv
