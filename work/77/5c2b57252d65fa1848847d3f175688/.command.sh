#!/usr/bin/env bash -C -e -u -o pipefail
python /home/lgonlopez/storage/phipsiprot/PhiPsiProt/bin/merge_screening_results.py         --candidates ranked_candidates.csv         --selection_summary selection_summary.csv         --output final_screening_results.csv
