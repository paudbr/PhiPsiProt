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
parser.add_argument("--pdbs", nargs="+", required=True)
parser.add_argument("--csv", required=True)
parser.add_argument("--fasta", required=True)
parser.add_argument("--samplesheet", required=True)
args = parser.parse_args()

pdbs = [Path(p) for p in args.pdbs]

with open(args.csv, "w", newline="") as out_csv, \
     open(args.fasta, "w") as out_fasta, \
     open(args.samplesheet, "w", newline="") as out_samplesheet:

    csv_writer = csv.writer(out_csv)
    sample_writer = csv.writer(out_samplesheet)

    csv_writer.writerow([
        "candidate_id",
        "mode",
        "method",
        "pdb",
        "sequence",
        "status",
    ])

    sample_writer.writerow(["id", "fasta"])

    for pdb in sorted(pdbs):
        candidate_id = pdb.stem
        sequence = extract_sequence(pdb)

        csv_writer.writerow([
            candidate_id,
            "backbone",
            "rfdiffusion",
            pdb.name,
            sequence,
            "ready_for_sequence_design",
        ])

        out_fasta.write(f">{candidate_id}|rfdiffusion\n")
        out_fasta.write(sequence + "\n")

        sample_writer.writerow([
            candidate_id,
            Path(args.fasta).name,
        ])

print(f"Backbone summary written for {len(pdbs)} candidates")