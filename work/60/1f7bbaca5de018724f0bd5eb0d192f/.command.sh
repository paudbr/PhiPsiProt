#!/usr/bin/env bash -C -e -u -o pipefail
python - <<'PY'
import csv

final_columns = [
    "candidate_id",
    "final_sequence",
    "mpnn_score",
    "molecular_weight",
    "gravy",
    "aromaticity",
    "isoelectric_point",
    "instability_index",
    "status",
]

with open("biophysical_filtered.csv", newline="") as infile:
    reader = csv.DictReader(infile)
    rows = []


    for row in reader:
    	row["status"] = "ready_for_structural_validation"

    	clean = {}
    	for col in final_columns:
        	clean[col] = row.get(col, "NA")
    	rows.append(clean)

with open("final_peptide_candidates.csv", "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=final_columns)
    writer.writeheader()
    writer.writerows(rows)
PY
