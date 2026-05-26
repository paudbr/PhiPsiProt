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

STANDARD_AA = list("ACDEFGHIKLMNPQRSTVWY")


def parse_pdb_residues(pdb_path):
    residues = []
    seen = set()

    with open(pdb_path) as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue

            resname = line[17:20].strip()
            chain = line[21].strip()
            resnum = line[22:26].strip()
            icode = line[26].strip()

            if resname not in AA3_TO_1:
                continue

            key = (chain, resnum, icode)

            if key in seen:
                continue

            seen.add(key)

            residues.append({
                "chain": chain,
                "position": resnum,
                "icode": icode,
                "wildtype": AA3_TO_1[resname],
                "resname": resname,
            })

    return residues


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--input_pdb", required=True)
    parser.add_argument("--target_chain", default="ALL")
    parser.add_argument("--positions", default="ALL")
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    pdb = Path(args.input_pdb)

    residues = parse_pdb_residues(pdb)

    # Filter by chain
    if args.target_chain != "ALL":
        residues = [
            r for r in residues
            if r["chain"] == args.target_chain
        ]

    # Filter by positions
    if args.positions != "ALL":

        selected_positions = set(
            p.strip() for p in args.positions.split(",")
        )

        residues = [
            r for r in residues
            if r["position"] in selected_positions
        ]

    with open(args.output, "w", newline="") as out:

        writer = csv.writer(out)

        writer.writerow([
            "candidate_id",
            "chain",
            "position",
            "wildtype",
            "mutant",
            "mutation",
            "ddg",
            "input_pdb",
        ])

        counter = 1

        for res in residues:

            wt = res["wildtype"]

            for mut in STANDARD_AA:

                if mut == wt:
                    continue

                mutation = (
                    f"{wt}"
                    f"{res['chain']}"
                    f"{res['position']}"
                    f"{mut}"
                )

                candidate_id = f"sat_{counter:06d}"

                writer.writerow([
                    candidate_id,
                    res["chain"],
                    res["position"],
                    wt,
                    mut,
                    mutation,
                    "NA",
                    pdb.name,
                ])

                counter += 1


if __name__ == "__main__":
    main()
