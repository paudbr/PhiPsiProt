#!/usr/bin/env python3

import argparse
import csv

parser = argparse.ArgumentParser()

parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)

args = parser.parse_args()

rows = []

with open(args.input, newline="") as infile:

    reader = csv.DictReader(infile)

    fieldnames = reader.fieldnames + [
        "structural_score",
        "final_score",
        "status"
    ]

    for row in reader:

        ddg = row.get("ddg", "NA")

        if ddg == "NA":
            ddg_score = "NA"
        else:
            try:
                ddg_score = float(ddg)
            except:
                ddg_score = "NA"

        row["structural_score"] = "NA"
        row["final_score"] = ddg_score
        row["status"] = "pending_structural_validation"

        rows.append(row)

with open(args.output, "w", newline="") as outfile:

    writer = csv.DictWriter(outfile, fieldnames=fieldnames)

    writer.writeheader()
    writer.writerows(rows)

print(f"Ranked candidates: {len(rows)}")
