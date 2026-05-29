#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--candidates_csv", required=True)
parser.add_argument("--mpnn_fastas", nargs="+", required=True)
parser.add_argument("--output_csv", required=True)
args = parser.parse_args()

def parse_mpnn_fasta(fasta_path):
    fasta_path = Path(fasta_path)
    candidate_id = fasta_path.stem

    sequence = "NA"
    score = "NA"

    lines = [x.strip() for x in fasta_path.read_text().splitlines() if x.strip()]

    for i, line in enumerate(lines):
        if line.startswith(">") and "id=1" in line:
            if i + 1 < len(lines):
                full_sequence = lines[i + 1]
                sequence = full_sequence.split(":")[0]

            if "overall_confidence=" in line:
                try:
                    score = line.split("overall_confidence=")[1].split(",")[0]
                except Exception:
                    score = "NA"

            break

    return candidate_id, sequence, score




mpnn_by_candidate = {}

for fasta in args.mpnn_fastas:
    candidate_id, sequence, score = parse_mpnn_fasta(fasta)
    mpnn_by_candidate[candidate_id] = {
        "binder_sequence": sequence,
        "mpnn_score": score,
    }


rows = []

with open(args.candidates_csv, newline="") as infile:
    reader = csv.DictReader(infile)

    for row in reader:
        candidate_id = row["candidate_id"]
        mpnn = mpnn_by_candidate.get(candidate_id, {})

        binder_sequence = mpnn.get("binder_sequence", "NA")

        row["binder_sequence"] = binder_sequence
        row["final_sequence"] = binder_sequence
        row["binder_length"] = len(binder_sequence) if binder_sequence != "NA" else "NA"
        row["mpnn_score"] = mpnn.get("mpnn_score", "NA")
        row["status"] = "pending_biophysical_filter"

        rows.append(row)


fieldnames = [
    "candidate_id",
    "pdb",
    "binder_length",
    "binder_sequence",
    "final_sequence",
    "mpnn_score",
    "status",
]

with open(args.output_csv, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=fieldnames)
    writer.writeheader()

    for row in rows:
        writer.writerow({col: row.get(col, "NA") for col in fieldnames})

print("Merged {} binder candidates".format(len(rows)))
