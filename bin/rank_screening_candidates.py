#!/usr/bin/env python3

import argparse
import csv
import math


def to_float(value):
    try:
        if value in ["", "NA", None]:
            return None
        return float(value)
    except Exception:
        return None

def compute_final_score(row):
    """
    Lower final_score is better.

    The score combines:
    - PyRosetta ddG as the main stability term
    - delta_solubility_score relative to WT
    - delta_instability_index relative to WT
    """

    ddg = to_float(row.get("ddg"))
    delta_solubility = to_float(row.get("delta_solubility_score"))
    delta_instability = to_float(row.get("delta_instability_index"))

    if ddg is None:
        return "NA"

    score = ddg

    # Reward increased predicted solubility relative to WT.
    if delta_solubility is not None:
        score -= 5.0 * delta_solubility

    # Penalize increased instability relative to WT.
    if delta_instability is not None and delta_instability > 0:
        score += 0.05 * delta_instability

    return round(score, 4)




parser = argparse.ArgumentParser()
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
        final_score = compute_final_score(row)

        row["structural_score"] = row.get("structural_score", "NA")
        row["final_score"] = final_score
        row["status"] = "pending_structural_validation"

        rows.append(row)

rows.sort(
    key=lambda r: (
        math.inf if r["final_score"] == "NA" else float(r["final_score"])
    )
)

with open(args.output, "w", newline="") as outfile:
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f"Ranked candidates: {len(rows)}")