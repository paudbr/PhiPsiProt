#!/usr/bin/env python3

"""
Generate PyMOL images for PhiPsiProt Mode 1 screening.

Outputs
-------
target_structure.png:
    Overview of the input protein structure.

top_candidate_mutation.png:
    Input protein structure with the top-ranked mutation position highlighted.
"""

import argparse
import csv
import subprocess
from pathlib import Path


def read_top_candidate(results_csv):
    """Read the first ranked candidate from final_screening_results.csv."""
    with open(results_csv, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            return row

    raise ValueError("No candidates found in results CSV")


def write_pml(path, commands):
    """Write a PyMOL command file."""
    with open(path, "w") as handle:
        handle.write("\n".join(commands))
        handle.write("\n")


def run_pymol(pml_file):
    """Run PyMOL in quiet batch mode."""
    subprocess.run(
        ["pymol", "-cq", str(pml_file)],
        check=True,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Generate PyMOL images for screening report."
    )
    parser.add_argument("--input_pdb", required=True)
    parser.add_argument("--results_csv", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    input_pdb = Path(args.input_pdb).resolve()
    results_csv = Path(args.results_csv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    top = read_top_candidate(results_csv)

    chain = top.get("chain", "A")
    position = top.get("position", "")
    mutation = top.get("mutation", "top_candidate")

    target_pml = outdir / "target_structure.pml"
    mutation_pml = outdir / "top_candidate_mutation.pml"

    target_png = outdir / "target_structure.png"
    mutation_png = outdir / "top_candidate_mutation.png"

    write_pml(
        target_pml,
        [
            f"load {input_pdb}, protein",
            "hide everything",
            "show cartoon, protein",
            "color slate, protein",
            "bg_color white",
            "set ray_opaque_background, off",
            "orient protein",
            "zoom protein, 8",
            "ray 1400, 1000",
            f"png {target_png}, dpi=300",
            "quit",
        ],
    )

    write_pml(
        mutation_pml,
        [
            f"load {input_pdb}, protein",
            "hide everything",
            "show cartoon, protein",
            "color gray80, protein",
            f"select mut_site, chain {chain} and resi {position}",
            "show sticks, mut_site",
            "color red, mut_site",
            "label mut_site and name CA, \"%s\" % resi",
            "set label_size, 24",
            "set label_color, black",
            "bg_color white",
            "set ray_opaque_background, off",
            "orient protein",
            "zoom mut_site, 10",
            "ray 1400, 1000",
            f"png {mutation_png}, dpi=300",
            "quit",
        ],
    )

    run_pymol(target_pml)
    run_pymol(mutation_pml)

    print(f"Generated PyMOL images for top candidate: {mutation}")


if __name__ == "__main__":
    main()
