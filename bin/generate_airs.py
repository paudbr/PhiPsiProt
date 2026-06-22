#!/usr/bin/env python3
import math
import sys

receptor_pdb   = sys.argv[1]
ligand_pdb     = sys.argv[2]
binding_site   = sys.argv[3]  # "168,232,58,241,67,369,97"

active_rec = [int(r.strip()) for r in binding_site.split(",") if r.strip()]

def parse_ca(pdb):
    residues = {}
    with open(pdb) as f:
        for l in f:
            if l.startswith("ATOM") and l[12:16].strip() == "CA":
                resi = int(l[22:26])
                coord = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
                residues[resi] = coord
    return residues

def dist(a, b):
    return math.sqrt(sum((a[i]-b[i])**2 for i in range(3)))

rec_ca = parse_ca(receptor_pdb)
lig_ca = parse_ca(ligand_pdb)

passive_rec = set()
for resi, coord in rec_ca.items():
    if resi in active_rec:
        continue
    for act in active_rec:
        if act in rec_ca and dist(coord, rec_ca[act]) <= 6.5:
            passive_rec.add(resi)
            break

passive_lig = sorted(lig_ca.keys())

with open("active_residues.txt", "w") as f:
    for i in sorted(active_rec):
        f.write(f"A{i}\n")

with open("passive_residues.txt", "w") as f:
    for i in sorted(passive_rec):
        f.write(f"A{i}\n")
    for i in passive_lig:
        f.write(f"B{i}\n")

with open("restraints.tbl", "w") as out:
    out.write("! HADDOCK AIRs - receptor active vs binder passive\n")

    lig_sel = " or\n        ".join(
        f"(segid B and resid {i} and name CA)" for i in passive_lig
    )

    for i in sorted(active_rec):
        out.write(
            f"assign (segid A and resid {i} and name CA)\n"
            f"       (\n        {lig_sel}\n       ) 2.0 2.0 0.0\n"
        )

    rec_sel = " or\n        ".join(
        f"(segid A and resid {i} and name CA)" for i in sorted(passive_rec)
    )
    if rec_sel:
        for i in passive_lig:
            out.write(
                f"assign (segid B and resid {i} and name CA)\n"
                f"       (\n        {rec_sel}\n       ) 2.0 2.0 0.0\n"
            )

print(f"Active receptor residues: {len(active_rec)}")
print(f"Passive receptor residues: {len(passive_rec)}")
print(f"Passive binder residues: {len(passive_lig)}")