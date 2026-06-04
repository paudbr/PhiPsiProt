#!/usr/bin/env bash -C -e -u -o pipefail
python /home/lgonlopez/storage/phipsiprot/PhiPsiProt/bin/rank_screening_candidates.py         --input biophysical_filtered.csv         --output ranked_candidates.csv
