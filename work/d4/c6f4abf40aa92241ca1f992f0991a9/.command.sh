#!/usr/bin/env bash -C -e -u -o pipefail
python3 /home/lgonlopez/storage/phipsiprot/PhiPsiProt/modules/local/biophysical_filter/bin/biophysical_filter.py         --input_csv final_backbone_results.csv         --output_csv biophysical_filtered.csv         --filters "instability,gravy,isoelectric_point"
