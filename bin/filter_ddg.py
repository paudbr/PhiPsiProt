#!/usr/bin/env python3

import argparse
import csv

parser = argparse.ArgumentParser()

parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
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

print(f"Filtered candidates: {len(rows)}")
