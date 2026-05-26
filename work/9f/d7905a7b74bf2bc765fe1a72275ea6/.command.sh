#!/usr/bin/env bash -C -e -u -o pipefail
echo "id,fasta,pdb,mode" > candidates.csv
echo "dummy_peptide_design,dummy.fasta,dummy.pdb,peptide_design" >> candidates.csv
