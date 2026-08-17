#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

parser = argparse.ArgumentParser()

parser.add_argument("--input_pdb", required=True)
parser.add_argument("--output_csv", required=True)
parser.add_argument("--output_fasta", required=True)
parser.add_argument("--summary", required=True)

args = parser.parse_args()

pdb = Path(args.input_pdb)

with open(args.output_csv, "w", newline="") as out:
    writer = csv.writer(out)
    writer.writerow([
        "candidate_id",
        "mode",
        "method",
        "input_pdb",
        "designed_sequence",
        "status"
    ])
    writer.writerow([
        "backbone_000001",
        "backbone",
        "rfdiffusion_placeholder",
        pdb.name,
        "MKTAYIAKQRQISFVKSHFSRQ",
        "pending_rfdiffusion"
    ])

with open(args.output_fasta, "w") as out:
    out.write(">backbone_000001|rfdiffusion_placeholder\n")
    out.write("MKTAYIAKQRQISFVKSHFSRQ\n")

with open(args.summary, "w", newline="") as out:
    writer = csv.writer(out)
    writer.writerow(["input_pdb", "n_candidates", "status"])
    writer.writerow([pdb.name, 1, "placeholder_completed"])

print("Backbone design placeholder completed")
