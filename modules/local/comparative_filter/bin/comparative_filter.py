#!/usr/bin/env python3

import argparse
import csv
from Bio.SeqUtils.ProtParam import ProteinAnalysis


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def extract_pdb_sequence(pdb_file, chain_id=None):
    sequence = []
    seen = set()

    with open(pdb_file) as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue

            resname = line[17:20].strip()
            chain = line[21].strip()
            resnum = line[22:26].strip()

            if chain_id and chain != chain_id:
                continue

            key = (chain, resnum)

            if key in seen:
                continue

            seen.add(key)
            sequence.append(AA3_TO_1.get(resname, "X"))

    return "".join(sequence)


def clean_sequence(seq):
    valid = set("ACDEFGHIKLMNPQRSTVWY")
    return "".join([aa for aa in seq.upper().strip() if aa in valid])


def instability(seq):
    seq = clean_sequence(seq)
    if not seq:
        return None
    return ProteinAnalysis(seq).instability_index()


parser = argparse.ArgumentParser()
parser.add_argument("--input_csv", required=True)
parser.add_argument("--input_pdb", required=True)
parser.add_argument("--output_csv", required=True)
parser.add_argument("--max_delta_instability", type=float, default=0.0)
args = parser.parse_args()

original_pdb_sequence = extract_pdb_sequence(args.input_pdb)

rows = []

with open(args.input_csv, newline="") as infile:
    reader = csv.DictReader(infile)

    for row in reader:
        original_sequence = original_pdb_sequence

        designed_sequence = (
            row.get("final_sequence")
            or row.get("designed_sequence")
            or ""
        )

        original_instability = instability(original_sequence)
        designed_instability = instability(designed_sequence)

        row["original_sequence"] = original_sequence

        if original_instability is None or designed_instability is None:
            row["original_instability_index"] = "NA"
            row["designed_instability_index"] = "NA"
            row["delta_instability_index"] = "NA"
            row["comparative_status"] = "failed"
        else:
            delta = designed_instability - original_instability

            row["original_instability_index"] = round(original_instability, 4)
            row["designed_instability_index"] = round(designed_instability, 4)
            row["delta_instability_index"] = round(delta, 4)

            if delta <= args.max_delta_instability:
                row["comparative_status"] = "pass"
            else:
                row["comparative_status"] = "fail"

        rows.append(row)

fieldnames = list(rows[0].keys())

with open(args.output_csv, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print("Processed {} comparative candidates".format(len(rows)))