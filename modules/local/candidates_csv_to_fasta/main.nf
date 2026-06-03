process CANDIDATES_CSV_TO_FASTA {

    tag "csv_to_fasta"

    input:
    path candidates_csv

    output:
    path "*.fasta"

    script:
    """
    python <<'PY'
import pandas as pd

df = pd.read_csv("${candidates_csv}")

for _, row in df.iterrows():
    cid = row["candidate_id"]
    seq = row["binder_sequence"]

    with open(f"{cid}.fasta", "w") as f:
        f.write(f">{cid}\\n")
        f.write(seq + "\\n")
PY
"""
}