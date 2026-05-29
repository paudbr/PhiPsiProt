#!/usr/bin/env bash -C -e -u -o pipefail
python3 /home/lgonlopez/storage/phipsiprot/PhiPsiProt/modules/local/comparative_filter/bin/comparative_filter.py         --input_csv biophysical_filtered.csv         --input_pdb 1PPN.pdb         --output_csv comparative_filtered.csv         --max_delta_instability 10
