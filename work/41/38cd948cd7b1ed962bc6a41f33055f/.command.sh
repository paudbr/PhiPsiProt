#!/usr/bin/env bash -C -e -u -o pipefail
python /home/lgonlopez/storage/phipsiprot/PhiPsiProt/bin/summarize_backbone_design.py         --pdbs backbone_design_0.pdb backbone_design_1.pdb backbone_design_2.pdb         --csv backbone_candidates.csv         --fasta backbone_candidates.fasta         --samplesheet backbone_samplesheet.csv
