#!/usr/bin/env python3

import argparse
import csv

parser = argparse.ArgumentParser()

parser.add_argument("--candidates", required=True)
parser.add_argument("--selection_summary", required=True)
parser.add_argument("--output", required=True)

args = parser.parse_args()

selection_metadata = {}

with open(args.selection_summary, newline="") as infile:
    reader = csv.DictReader(infile)
    for row in reader:
        selection_metadata = row
        break

rows = []

with open(args.candidates, newline="") as infile:
    reader = csv.DictReader(infile)

    fieldnames = reader.fieldnames + [
        "selection_mode",
        "selected_positions",
        "ligand_resname",
        "interface_chain",
        "distance_cutoff",
        "pipeline_stage",
    ]

    for row in reader:
        row["selection_mode"] = selection_metadata.get("selection_mode", "NA")
        row["selected_positions"] = selection_metadata.get("selected_positions", "NA")
        row["ligand_resname"] = selection_metadata.get("ligand_resname", "NA")
        row["interface_chain"] = selection_metadata.get("interface_chain", "NA")
        row["distance_cutoff"] = selection_metadata.get("distance_cutoff", "NA")
        row["pipeline_stage"] = "screening_ready_for_structural_validation"

        rows.append(row)

with open(args.output, "w", newline="") as outfile:
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f"Final screening results written: {len(rows)}")
