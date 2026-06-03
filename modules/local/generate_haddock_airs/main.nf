process GENERATE_HADDOCK_AIRS {

    tag "generate_airs"
    label 'process_medium'

    input:
    tuple val(candidate_id), path(receptor_pdb), path(ligand_pdb)

    output:
    tuple val(candidate_id),
          path("restraints.tbl"),
          path("active_residues.txt"),
          path("passive_residues.txt")

    script:
    """
    python <<'PY'
import math

CUT_OFF = 4.5

def dist(a,b):
    return math.sqrt(sum((a[i]-b[i])**2 for i in range(3)))

def parse(pdb):
    atoms = []
    with open(pdb) as f:
        for l in f:
            if l.startswith(("ATOM","HETATM")):
                atoms.append({
                    "coord": (
                        float(l[30:38]),
                        float(l[38:46]),
                        float(l[46:54])
                    ),
                    "chain": l[21].strip() or "A",
                    "resi": int(l[22:26])
                })
    return atoms

rec = parse("${receptor_pdb}")
lig = parse("${ligand_pdb}")

active = set()

for r in rec:
    for l in lig:
        if dist(r["coord"], l["coord"]) <= CUT_OFF:
            active.add((r["chain"], r["resi"]))

passive = set()

for r in rec:
    for (c,i) in active:
        if r["chain"] == c and abs(r["resi"]-i) <= 1:
            passive.add((r["chain"], r["resi"]))

with open("active_residues.txt","w") as f:
    for c,i in sorted(active):
        f.write(f"{c}{i}\\n")

with open("passive_residues.txt","w") as f:
    for c,i in sorted(passive):
        f.write(f"{c}{i}\\n")

with open("restraints.tbl","w") as out:
    out.write("! HADDOCK AIRs\\n")
    for c,i in active:
        out.write(
            f"assign (segid PROT and chain {c} and resid {i} and name CA) "
            f"(segid LIG and chain {c} and resid {i} and name CA) "
            f"2.0 2.0 0.0\\n"
        )

print("ACTIVE:", len(active))
PY
"""
}