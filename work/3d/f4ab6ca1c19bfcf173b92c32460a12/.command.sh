#!/usr/bin/env bash -C -e -u -o pipefail
mkdir -p plots
cp ddg_heatmap.png ddg_ranking_barplot.png top10_stabilizing_mutations.png plots/ || true

python3 /home/lgonlopez/storage/phipsiprot/PhiPsiProt/bin/generate_screening_report.py         --results_csv final_screening_results.csv         --plots_dir plots         --outdir .
