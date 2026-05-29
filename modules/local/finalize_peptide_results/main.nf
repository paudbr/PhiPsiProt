process FINALIZE_PEPTIDE_RESULTS {

    tag "finalize_peptide_results"

    publishDir "${params.outdir}/final", mode: 'copy'

    input:
    path input_csv

    output:
    path "final_peptide_candidates.csv"

    script:
    """
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

with open("$input_csv", newline="") as infile:
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
    """
}
