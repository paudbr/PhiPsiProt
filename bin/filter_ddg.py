#!/usr/bin/env python3

import argparse
import csv

parser = argparse.ArgumentParser()

parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--fasta_output", required=True)
parser.add_argument("--samplesheet_output", required=True)
parser.add_argument("--threshold", type=float, default=2.0)

args = parser.parse_args()

rows = []

with open(args.input, newline='') as infile:
    reader = csv.DictReader(infile)

    for row in reader:

        ddg = row["ddg"]

        if ddg == "NA" or ddg == "":
            rows.append(row)
            continue

        try:
            if float(ddg) <= args.threshold:
                rows.append(row)
        except:
            pass

with open(args.output, "w", newline='') as outfile:
    writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)

    writer.writeheader()
    writer.writerows(rows)

with open(args.fasta_output, "w") as fasta_out:

    for row in rows:

        fasta_id = row["fasta_id"]

        wt = row["wildtype"]
        mut = row["mutant"]

        sequence = f"{mut}"

        fasta_out.write(f">{fasta_id}\n")
        fasta_out.write(f"{sequence}\n")

with open(args.samplesheet_output, "w", newline='') as sample_out:

    writer = csv.writer(sample_out)

    writer.writerow(["id", "fasta"])

    for row in rows:

        fasta_id = row["fasta_id"]

        writer.writerow([
            fasta_id,
            "filtered_candidates.fasta"
        ])

print(f"Filtered candidates: {len(rows)}")
