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

    designed_sequence = None
    mpnn_score = None

    lines = [x.strip() for x in fasta_path.read_text().splitlines() if x.strip()]

    for i, line in enumerate(lines):
        if line.startswith(">T="):
            if i + 1 < len(lines):
                designed_sequence = lines[i + 1]

            if "score=" in line:
                mpnn_score = line.split("score=")[1].split(",")[0]

            break

    return candidate_id, designed_sequence, mpnn_score


mpnn_by_candidate = {}

for fasta in args.mpnn_fastas:
    candidate_id, designed_sequence, mpnn_score = parse_mpnn_fasta(fasta)

    mpnn_by_candidate[candidate_id] = {
        "designed_sequence": designed_sequence,
        "mpnn_score": mpnn_score,
    }


rows = []

with open(args.candidates_csv, newline="") as infile:
    reader = csv.DictReader(infile)

    for row in reader:
        candidate_id = row["candidate_id"]

        mpnn = mpnn_by_candidate.get(candidate_id, {})

        designed_sequence = mpnn.get("designed_sequence", "NA")
        mpnn_score = mpnn.get("mpnn_score", "NA")

        row["backbone_sequence"] = row.get("sequence", "NA")
        row["designed_sequence"] = designed_sequence
        row["final_sequence"] = designed_sequence
        row["mpnn_score"] = mpnn_score
        row["pipeline_stage"] = "backbone_ready_for_structural_validation"
        row["status"] = "pending_structural_validation"

        rows.append(row)

fieldnames = list(rows[0].keys())

with open(args.output_csv, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print("Merged {} backbone candidates".format(len(rows)))

