#!/usr/bin/env bash -C -e -u -o pipefail
python /home/lgonlopez/storage/phipsiprot/PhiPsiProt/bin/merge_binder_results.py         --candidates_csv binder_candidates.csv         --mpnn_fastas binder_design_0.fa binder_design_2.fa binder_design_1.fa         --output_csv final_binder_candidates.csv
