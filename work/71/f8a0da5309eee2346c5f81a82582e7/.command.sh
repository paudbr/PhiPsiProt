#!/usr/bin/env bash -C -e -u -o pipefail
python /home/lgonlopez/storage/phipsiprot/PhiPsiProt/bin/filter_ddg.py         --input candidates.csv         --output filtered_candidates.csv         --fasta_output filtered_candidates.fasta         --samplesheet_output screening_samplesheet.csv         --threshold 2.0
