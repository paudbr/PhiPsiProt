#!/usr/bin/env python3

"""
Filter PyRosetta screening candidates by ΔΔG.

Purpose
-------
Read candidates.csv from the PyRosetta screening step and retain candidates
with ddG values below or equal to the selected threshold.

Inputs
------
candidates.csv:
    Candidate table produced by run_pyrosetta_screening.py.

Outputs
-------
filtered_candidates.csv:
    Filtered candidate table.

filtered_candidates.fasta:
    FASTA file containing the full mutant sequence for each retained candidate.

screening_samplesheet.csv:
    Simple samplesheet linking candidate IDs to the filtered FASTA file.
"""

import argparse
import csv


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Filter mutational screening candidates by ddG."
    )

    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--fasta_output", required=True)
    parser.add_argument("--samplesheet_output", required=True)
    parser.add_argument("--threshold", type=float, default=2.0)

    return parser.parse_args()


def keep_candidate(row, threshold):
    """Return True if a candidate passes the ddG threshold."""
    ddg = row.get("ddg", "")

    if ddg in ("", "NA"):
        return True

    try:
        return float(ddg) <= threshold
    except ValueError:
        return False


def write_fasta_record(handle, header, sequence):
    """Write one FASTA record with 80-character line wrapping."""
    handle.write(f">{header}\n")
    for i in range(0, len(sequence), 80):
        handle.write(sequence[i:i + 80] + "\n")


def main():
    """Run ddG filtering."""
    args = parse_args()

    rows = []

    with open(args.input, newline="") as infile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames

        for row in reader:
            if keep_candidate(row, args.threshold):
                rows.append(row)

    with open(args.output, "w", newline="") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    with open(args.fasta_output, "w") as fasta_out:
        for row in rows:
            fasta_id = row["fasta_id"]
            sequence = row["sequence"]
            write_fasta_record(fasta_out, fasta_id, sequence)

    with open(args.samplesheet_output, "w", newline="") as sample_out:
        writer = csv.writer(sample_out)
        writer.writerow(["id", "fasta"])

        for row in rows:
            writer.writerow([
                row["fasta_id"],
                "filtered_candidates.fasta",
            ])

    print(f"Filtered candidates: {len(rows)}")


if __name__ == "__main__":
    main()