#!/usr/bin/env bash -C -e -u -o pipefail
python /home/lgonlopez/storage/phipsiprot/PhiPsiProt/bin/merge_backbone_results.py         --candidates_csv backbone_candidates.csv         --mpnn_fastas backbone_design_1.fa backbone_design_0.fa backbone_design_2.fa         --output_csv final_backbone_results.csv
