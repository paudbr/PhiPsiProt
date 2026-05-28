#!/usr/bin/env python3

import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--input_pdb", required=True)
parser.add_argument("--positions", required=True)
parser.add_argument("--target_chain", default="A")
parser.add_argument("--output_prefix", default="selected_residues")
args = parser.parse_args()

pdb = Path(args.input_pdb)
positions_file = Path(args.positions)

positions = positions_file.read_text().strip().replace(",", "+")
positions_label = positions.replace("+", ",")

pml = f"""
load /data/{pdb.name}, protein

hide everything
show cartoon, protein
color gray80, protein

select selected_residues, chain {args.target_chain} and resi {positions}
show sticks, selected_residues
color red, selected_residues
label selected_residues and name CA, resn + resi

set cartoon_transparency, 0.15
set ray_opaque_background, off
bg_color white

zoom selected_residues, 12

png /data/{args.output_prefix}.png, dpi=150, ray=0
save /data/{args.output_prefix}.pse
quit
"""

Path(f"{args.output_prefix}.pml").write_text(pml)

html = f"""<!DOCTYPE html>
<html>
<head>
<title>PyMOL residue selection report</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 40px; }}
img {{ max-width: 900px; border: 1px solid #ddd; }}
</style>
</head>
<body>
<h1>Residue selection visualization</h1>
<p><b>Input PDB:</b> {pdb.name}</p>
<p><b>Target chain:</b> {args.target_chain}</p>
<p><b>Selected residues:</b> {positions_label}</p>
<img src="{args.output_prefix}.png">
</body>
</html>
"""

Path(f"{args.output_prefix}.html").write_text(html)
