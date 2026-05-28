#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

parser = argparse.ArgumentParser()

parser.add_argument("--candidates_csv", required=True)
parser.add_argument("--mpnn_fasta", required=True)
parser.add_argument("--output_csv", required=True)

args = parser.parse_args()

# Leer secuencia diseñada ProteinMPNN
designed_sequence = None
mpnn_score = None

with open(args.mpnn_fasta) as handle:
    lines = [x.strip() for x in handle.readlines() if x.strip()]

for i, line in enumerate(lines):
    if line.startswith(">T="):
        designed_sequence = lines[i + 1]

        if "score=" in line:
            try:
                mpnn_score = line.split("score=")[1].split(",")[0]
            except Exception:
                mpnn_score = "NA"

        break

rows = []

with open(args.candidates_csv, newline="") as infile:
    reader = csv.DictReader(infile)

    for row in reader:
        row["designed_sequence"] = designed_sequence
        row["mpnn_score"] = mpnn_score
        row["pipeline_stage"] = "backbone_ready_for_structural_validation"
        row["status"] = "pending_structural_validation"

        rows.append(row)

fieldnames = list(rows[0].keys())

with open(args.output_csv, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f"Merged {len(rows)} backbone candidates")
