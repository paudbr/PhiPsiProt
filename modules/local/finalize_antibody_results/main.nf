process FINALIZE_ANTIBODY_RESULTS {

    tag "finalize_antibody_results"

    executor 'local'

    publishDir "${params.outdir}/final", mode: 'copy'

    input:
    path scores_tsv

    output:
    path "final_antibody_candidates.csv"

    script:
    """
    python3 - <<'PY'
import csv

input_file = "$scores_tsv"
output_file = "final_antibody_candidates.csv"

with open(input_file, newline="") as infile:
    first = infile.readline()
    infile.seek(0)

    delimiter = "\\t" if "\\t" in first else ","
    reader = csv.DictReader(infile, delimiter=delimiter)

    rows = []
    for i, row in enumerate(reader, start=1):
        clean = dict(row)
        clean["candidate_id"] = clean.get("tag", clean.get("description", "antibody_design_{}".format(i)))
        clean["status"] = "ready_for_review"
        rows.append(clean)

if rows:
    fieldnames = ["candidate_id"] + [c for c in rows[0].keys() if c != "candidate_id"]
else:
    fieldnames = ["candidate_id", "status"]

with open(output_file, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
PY
    """
}
