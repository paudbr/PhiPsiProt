#!/usr/bin/env python3

import argparse
import csv
from Bio.SeqUtils.ProtParam import ProteinAnalysis

parser = argparse.ArgumentParser()
parser.add_argument("--input_csv", required=True)
parser.add_argument("--output_csv", required=True)
parser.add_argument("--filters", default="all")

args = parser.parse_args()
selected_filters = [
    f.strip() for f in args.filters.split(",")
]

rows = []

with open(args.input_csv, newline="") as infile:
    reader = csv.DictReader(infile)

    for row in reader:

        sequence = (
            row.get("final_sequence")
            or row.get("designed_sequence")
            or row.get("sequence")
            or ""
        )

        sequence = sequence.strip()

        try:
            analysis = ProteinAnalysis(sequence)

            if "all" in selected_filters or "gravy" in selected_filters:
                row["gravy"] = round(analysis.gravy(), 4)

            if "all" in selected_filters or "aromaticity" in selected_filters:
                row["aromaticity"] = round(analysis.aromaticity(), 4)

            if "all" in selected_filters or "instability" in selected_filters:
                row["instability_index"] = round(analysis.instability_index(), 4)

            if "all" in selected_filters or "molecular_weight" in selected_filters:
                row["molecular_weight"] = round(analysis.molecular_weight(), 2)

            if "all" in selected_filters or "isoelectric_point" in selected_filters:
                row["isoelectric_point"] = round(analysis.isoelectric_point(), 2)
                        
        except Exception:

            row["gravy"] = "NA"
            row["aromaticity"] = "NA"
            row["instability_index"] = "NA"
            row["molecular_weight"] = "NA"
            row["isoelectric_point"] = "NA"
            row["biophysical_status"] = "failed"

        rows.append(row)

fieldnames = list(rows[0].keys())

with open(args.output_csv, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print("Processed {} sequences".format(len(rows)))
