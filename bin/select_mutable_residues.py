#!/usr/bin/env python3

import argparse
import math

parser = argparse.ArgumentParser()

parser.add_argument("--input_pdb", required=True)
parser.add_argument("--selection_mode", default="manual",
                    choices=["manual", "all", "pocket", "interface"])
parser.add_argument("--target_chain", required=True)
parser.add_argument("--positions", default="ALL")
parser.add_argument("--ligand_resname", default="")
parser.add_argument("--interface_chain", default="")
parser.add_argument("--distance_cutoff", type=float, default=6.0)
parser.add_argument("--output", required=True)

args = parser.parse_args()


def dist(a, b):
    return math.sqrt(
        (a[0] - b[0]) ** 2 +
        (a[1] - b[1]) ** 2 +
        (a[2] - b[2]) ** 2
    )


protein_atoms = []
ligand_atoms = []

with open(args.input_pdb) as handle:
    for line in handle:
        if not line.startswith(("ATOM", "HETATM")):
            continue

        resname = line[17:20].strip()
        chain = line[21].strip()
        resnum = line[22:26].strip()

        coords = (
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54]),
        )

        atom = {
            "resname": resname,
            "chain": chain,
            "resnum": resnum,
            "coords": coords,
        }

        if line.startswith("ATOM"):
            protein_atoms.append(atom)
        elif line.startswith("HETATM"):
            ligand_atoms.append(atom)


target_positions = sorted(
    {a["resnum"] for a in protein_atoms if a["chain"] == args.target_chain},
    key=lambda x: int(x)
)

selected = set()

if args.selection_mode == "all":
    selected = set(target_positions)

elif args.selection_mode == "manual":
    if args.positions == "ALL":
        selected = set(target_positions)
    else:
        selected = set(p.strip() for p in args.positions.split(","))

elif args.selection_mode == "pocket":
    if not args.ligand_resname:
        raise ValueError("--ligand_resname is required for selection_mode pocket")

    ligand_selected = [
        a for a in ligand_atoms
        if a["resname"] == args.ligand_resname
    ]

    for p in protein_atoms:
        if p["chain"] != args.target_chain:
            continue

        for l in ligand_selected:
            if dist(p["coords"], l["coords"]) <= args.distance_cutoff:
                selected.add(p["resnum"])
                break

elif args.selection_mode == "interface":
    if not args.interface_chain:
        raise ValueError("--interface_chain is required for selection_mode interface")

    interface_atoms = [
        a for a in protein_atoms
        if a["chain"] == args.interface_chain
    ]

    for p in protein_atoms:
        if p["chain"] != args.target_chain:
            continue

        for q in interface_atoms:
            if dist(p["coords"], q["coords"]) <= args.distance_cutoff:
                selected.add(p["resnum"])
                break

selected = sorted(selected, key=lambda x: int(x))

with open(args.output, "w") as out:
    out.write(",".join(selected))
    out.write("\n")
