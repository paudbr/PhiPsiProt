#!/usr/bin/env bash -C -e -u -o pipefail
python3 /home/lgonlopez/storage/phipsiprot/PhiPsiProt/bin/plot_screening_results.py         --input ranked_candidates.csv         --outdir .
