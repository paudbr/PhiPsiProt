#!/usr/bin/env python3

import argparse
import math

parser = argparse.ArgumentParser()
parser.add_argument("--input_pdb", required=True)
parser.add_argument("--ligand_resname", required=True)
parser.add_argument("--distance_cutoff", type=float, default=6.0)
parser.add_argument("--output", required=True)
args = parser.parse_args()

protein_atoms = []
ligand_atoms = []

with open(args.input_pdb) as handle:
    for line in handle:
        if not line.startswith(("ATOM", "HETATM")):
            continue

        resname = line[17:20].strip()
        chain = line[21].strip()
        resnum = line[22:26].strip()

        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

        atom = {
            "resname": resname,
            "chain": chain,
            "resnum": resnum,
            "coords": (x, y, z),
        }

        if line.startswith("ATOM"):
            protein_atoms.append(atom)
        elif resname == args.ligand_resname:
            ligand_atoms.append(atom)

selected = set()

for p in protein_atoms:
    px, py, pz = p["coords"]

    for l in ligand_atoms:
        lx, ly, lz = l["coords"]

        dist = math.sqrt((px - lx) ** 2 + (py - ly) ** 2 + (pz - lz) ** 2)

        if dist <= args.distance_cutoff:
            selected.add(p["resnum"])
            break

with open(args.output, "w") as out:
    out.write(",".join(sorted(selected, key=lambda x: int(x))))
    out.write("\n")
