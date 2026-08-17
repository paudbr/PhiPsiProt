#!/usr/bin/env python3

"""
Annotate candidate protein sequences with biophysical descriptors.

Purpose
-------
Compute sequence-derived developability metrics and append them
to the candidate table without removing any candidate.

Computed descriptors
--------------------
- GRAVY
- aromaticity
- instability index
- molecular weight
- isoelectric point
- solubility proxy

Additional outputs
------------------
- delta values relative to the WT/reference sequence
- biophysical_status
"""

import argparse
import csv
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def clean_sequence(sequence):
    sequence = (sequence or "").strip()

    if sequence.startswith("[") and sequence.endswith("]"):
        sequence = sequence.strip("[]").strip().strip("'").strip('"')

    sequence = sequence.replace("'", "").replace('"', "").replace(" ", "")
    return sequence


def solubility_proxy(gravy):
    return 1 / (1 + abs(gravy))


def compute_properties(sequence):
    analysis = ProteinAnalysis(sequence)

    gravy = analysis.gravy()
    aromaticity = analysis.aromaticity()
    instability = analysis.instability_index()
    molecular_weight = analysis.molecular_weight()
    isoelectric_point = analysis.isoelectric_point()
    solubility = solubility_proxy(gravy)

    return {
        "gravy": gravy,
        "aromaticity": aromaticity,
        "instability_index": instability,
        "molecular_weight": molecular_weight,
        "isoelectric_point": isoelectric_point,
        "solubility_score": solubility,
    }


parser = argparse.ArgumentParser()
parser.add_argument("--input_csv", required=True)
parser.add_argument("--output_csv", required=True)
parser.add_argument("--filters", default="all")
args = parser.parse_args()

selected_filters = [f.strip() for f in args.filters.split(",")]

rows = []

with open(args.input_csv, newline="") as infile:
    reader = csv.DictReader(infile)
    input_fieldnames = list(reader.fieldnames or [])

    raw_rows = list(reader)

if not raw_rows:
    fieldnames = input_fieldnames[:]
    for extra in [
        "gravy",
        "aromaticity",
        "instability_index",
        "molecular_weight",
        "isoelectric_point",
        "solubility_score",
        "delta_gravy",
        "delta_aromaticity",
        "delta_instability_index",
        "delta_molecular_weight",
        "delta_isoelectric_point",
        "delta_solubility_score",
        "biophysical_status",
    ]:
        if extra not in fieldnames:
            fieldnames.append(extra)

    with open(args.output_csv, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()

    print("Processed 0 sequences")
    raise SystemExit(0)


# Use the first sequence in the CSV as WT-like reference for delta metrics.
# In Mode 1, all candidates come from the same input protein, so this provides
# a consistent baseline for comparing biophysical changes.
first_sequence = clean_sequence(
    raw_rows[0].get("wt_sequence")
    or raw_rows[0].get("original_sequence")
    or raw_rows[0].get("final_sequence")
    or raw_rows[0].get("designed_sequence")
    or raw_rows[0].get("sequence")
    or ""
)

try:
    wt_properties = compute_properties(first_sequence)
except Exception:
    wt_properties = None


for row in raw_rows:

    sequence = clean_sequence(
        row.get("final_sequence")
        or row.get("designed_sequence")
        or row.get("sequence")
        or ""
    )

    try:
        props = compute_properties(sequence)

        if "all" in selected_filters or "gravy" in selected_filters:
            row["gravy"] = round(props["gravy"], 4)

        if "all" in selected_filters or "aromaticity" in selected_filters:
            row["aromaticity"] = round(props["aromaticity"], 4)

        if "all" in selected_filters or "instability" in selected_filters:
            row["instability_index"] = round(props["instability_index"], 4)

        if "all" in selected_filters or "molecular_weight" in selected_filters:
            row["molecular_weight"] = round(props["molecular_weight"], 2)

        if "all" in selected_filters or "isoelectric_point" in selected_filters:
            row["isoelectric_point"] = round(props["isoelectric_point"], 2)

        if "all" in selected_filters or "solubility" in selected_filters:
            row["solubility_score"] = round(props["solubility_score"], 4)

        if wt_properties is not None:
            row["delta_gravy"] = round(
                props["gravy"] - wt_properties["gravy"], 4
            )
            row["delta_aromaticity"] = round(
                props["aromaticity"] - wt_properties["aromaticity"], 4
            )
            row["delta_instability_index"] = round(
                props["instability_index"] - wt_properties["instability_index"], 4
            )
            row["delta_molecular_weight"] = round(
                props["molecular_weight"] - wt_properties["molecular_weight"], 2
            )
            row["delta_isoelectric_point"] = round(
                props["isoelectric_point"] - wt_properties["isoelectric_point"], 2
            )
            row["delta_solubility_score"] = round(
                props["solubility_score"] - wt_properties["solubility_score"], 4
            )
        else:
            row["delta_gravy"] = "NA"
            row["delta_aromaticity"] = "NA"
            row["delta_instability_index"] = "NA"
            row["delta_molecular_weight"] = "NA"
            row["delta_isoelectric_point"] = "NA"
            row["delta_solubility_score"] = "NA"

        row["biophysical_status"] = "ok"

    except Exception:
        row["gravy"] = "NA"
        row["aromaticity"] = "NA"
        row["instability_index"] = "NA"
        row["molecular_weight"] = "NA"
        row["isoelectric_point"] = "NA"
        row["solubility_score"] = "NA"
        row["delta_gravy"] = "NA"
        row["delta_aromaticity"] = "NA"
        row["delta_instability_index"] = "NA"
        row["delta_molecular_weight"] = "NA"
        row["delta_isoelectric_point"] = "NA"
        row["delta_solubility_score"] = "NA"
        row["biophysical_status"] = "failed"

    rows.append(row)


fieldnames = list(rows[0].keys())

with open(args.output_csv, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print("Processed {} sequences".format(len(rows)))