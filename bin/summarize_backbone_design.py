#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def extract_sequence(pdb_file):
    sequence = []
    seen = set()

    with open(pdb_file) as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue

            resname = line[17:20].strip()
            chain = line[21].strip()
            resnum = line[22:26].strip()

            key = (chain, resnum)

            if key in seen:
                continue

            seen.add(key)
            sequence.append(AA3_TO_1.get(resname, "X"))

    return "".join(sequence)


parser = argparse.ArgumentParser()
parser.add_argument("--pdb", required=True)
parser.add_argument("--csv", required=True)
parser.add_argument("--fasta", required=True)
parser.add_argument("--samplesheet", required=True)
args = parser.parse_args()

pdb = Path(args.pdb)
candidate_id = pdb.stem
sequence = extract_sequence(pdb)

with open(args.csv, "w", newline="") as out:
    writer = csv.writer(out)
    writer.writerow([
        "candidate_id",
        "mode",
        "method",
        "pdb",
        "sequence",
        "status",
    ])
    writer.writerow([
        candidate_id,
        "backbone",
        "rfdiffusion",
        pdb.name,
        sequence,
        "ready_for_sequence_design",
    ])

with open(args.fasta, "w") as out:
    out.write(f">{candidate_id}|rfdiffusion\n")
    out.write(sequence + "\n")

with open(args.samplesheet, "w", newline="") as out:
    writer = csv.writer(out)
    writer.writerow(["id", "fasta"])
    writer.writerow([candidate_id, Path(args.fasta).name])

print(f"Backbone summary written for {candidate_id}")
