#!/usr/bin/env python3

"""
Annotate PyRosetta screening candidates with ddG threshold information.

Purpose
-------
Read candidates.csv and add ddG-based annotation columns without removing
any candidate.

Inputs
------
candidates.csv:
    Candidate table produced by run_pyrosetta_screening.py.

Outputs
-------
candidates_ddg_annotated.csv:
    Same candidates as input, with additional ddG annotation columns.

Added columns
-------------
ddg_threshold:
    Threshold used to evaluate candidate stability.

ddg_pass:
    True if ddg <= threshold, False otherwise.

ddg_filter_status:
    Human-readable status describing the ddG evaluation.
"""

import argparse
import csv


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Annotate screening candidates with ddG threshold information."
    )

    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--threshold", type=float, default=2.0)

    return parser.parse_args()


def annotate_row(row, threshold):
    """Add ddG annotation fields to one candidate row."""
    ddg_value = row.get("ddg", "")

    row["ddg_threshold"] = threshold

    if ddg_value in ("", "NA"):
        row["ddg_pass"] = "NA"
        row["ddg_filter_status"] = "missing_ddg"
        return row

    try:
        ddg = float(ddg_value)
    except ValueError:
        row["ddg_pass"] = "NA"
        row["ddg_filter_status"] = "invalid_ddg"
        return row

    if ddg <= threshold:
        row["ddg_pass"] = "True"
        row["ddg_filter_status"] = "pass"
    else:
        row["ddg_pass"] = "False"
        row["ddg_filter_status"] = "fail"

    return row


def main():
    """Run ddG annotation."""
    args = parse_args()

    rows = []

    with open(args.input, newline="") as infile:
        reader = csv.DictReader(infile)
        original_fieldnames = reader.fieldnames or []

        for row in reader:
            rows.append(annotate_row(row, args.threshold))

    added_fieldnames = [
        "ddg_threshold",
        "ddg_pass",
        "ddg_filter_status",
    ]

    fieldnames = original_fieldnames + [
        field for field in added_fieldnames
        if field not in original_fieldnames
    ]

    with open(args.output, "w", newline="") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    total = len(rows)
    passed = sum(1 for row in rows if row.get("ddg_pass") == "True")

    print(f"Annotated candidates: {total}")
    print(f"Candidates passing ddG threshold: {passed}")


if __name__ == "__main__":
    main()