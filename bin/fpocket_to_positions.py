#!/usr/bin/env python3

"""
Convert fpocket pocket atom files into residue position lists.

fpocket writes pocket*_atm.pdb files containing protein atoms that define each
predicted pocket. This script extracts residue numbers for a selected chain and
writes them as a comma-separated list for PhiPsiProt Mode 1 screening.
"""

import argparse
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("--pocket_pdb", required=True)
parser.add_argument("--target_chain", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()


positions = set()

with open(args.pocket_pdb) as handle:
    for line in handle:
        if not line.startswith("ATOM"):
            continue

        chain = line[21].strip()
        resnum = line[22:26].strip()

        if chain == args.target_chain:
            positions.add(resnum)


positions = sorted(positions, key=lambda x: int(x))

with open(args.output, "w") as out:
    out.write(",".join(positions))
    out.write("\n")
