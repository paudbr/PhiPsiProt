#!/usr/bin/env python3


"""
Rank PhiPsiProt screening candidates.

Purpose
-------
Rank candidates using PyRosetta ddG and sequence-based biophysical
annotations.

Scoring
-------
Lower final_score is better.

Current formula
---------------
final_score = ddg

Notes
-----
- No candidates are removed.
- Candidates are sorted by final_score.
- structural_score is reserved for future structural validation methods.
- status is currently set to 'pending_structural_validation'.
"""

import argparse
import csv
import math


def to_float(value):
    """Convert values to float, returning None for missing or invalid values."""
    try:
        if value in ["", "NA", None]:
            return None
        return float(value)
    except Exception:
        return None


def compute_final_score(row):
    """
    Compute candidate ranking score.

    Lower final_score indicates a more favorable candidate.
    """
    ddg = to_float(row.get("ddg"))

    if ddg is None:
        return "NA"

    return round(ddg, 4)

def sort_key(row):
    """Sort candidates by final_score, placing NA values last."""
    final_score = row.get("final_score", "NA")

    if final_score == "NA":
        return math.inf

    return float(final_score)


def main():
    """Run candidate ranking."""
    parser = argparse.ArgumentParser(
        description="Rank PhiPsiProt screening candidates."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rows = []

    with open(args.input, newline="") as infile:
        reader = csv.DictReader(infile)
        base_fieldnames = reader.fieldnames or []

        fieldnames = list(base_fieldnames)
        for extra in ["structural_score", "final_score", "status"]:
            if extra not in fieldnames:
                fieldnames.append(extra)

        for row in reader:
            row["structural_score"] = row.get("structural_score", "NA")
            row["final_score"] = compute_final_score(row)
            row["status"] = "pending_structural_validation"
            rows.append(row)

    rows.sort(key=sort_key)

    with open(args.output, "w", newline="") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Ranked candidates: {len(rows)}")


if __name__ == "__main__":
    main()