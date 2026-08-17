#!/usr/bin/env python3

"""
Run PyRosetta mutational screening for PhiPsiProt Mode 1.

Purpose
-------
Generate all single amino-acid substitutions for the selected residue
positions and estimate their energetic effect using a PyRosetta score
function.

Inputs
------
input_pdb:
    Input PDB structure.

target_chain:
    Chain to mutate.

positions:
    Comma-separated residue positions or ALL.

Outputs
-------
candidates.csv:
    Table with one row per single-point mutant.

candidates.fasta:
    FASTA file containing the mutant sequence for each candidate.

Notes
-----
- The script generates single mutants only.
- ΔΔG is approximated as mutant_score - wildtype_score.
- More negative ΔΔG values are interpreted as more stabilizing.
- HETATM records are preserved during PDB cleaning so cofactors such as
  metal ions or ligands are not removed before scoring.
"""

import argparse
import csv
from pathlib import Path

import pyrosetta
from pyrosetta.toolbox import mutate_residue


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

STANDARD_AA = list("ACDEFGHIKLMNPQRSTVWY")


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run PyRosetta mutational screening for PhiPsiProt."
    )

    parser.add_argument("--input_pdb", required=True)
    parser.add_argument("--target_chain", required=True)
    parser.add_argument("--positions", default="ALL")
    parser.add_argument("--output", required=True)
    parser.add_argument("--fasta_output", required=True)

    return parser.parse_args()


def parse_pdb_residues(pdb_path):
    """Parse standard protein residues from ATOM records."""
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
            residues.append(
                {
                    "chain": chain,
                    "position": resnum,
                    "icode": icode,
                    "wildtype": AA3_TO_1[resname],
                    "resname": resname,
                }
            )

    return residues


def clean_pdb_for_pyrosetta(input_pdb, output_pdb):
    """
    Create a PyRosetta-readable PDB.

    Keep ATOM and HETATM records so that cofactors, metal ions and ligands
    are not removed before scoring.
    """
    with open(input_pdb) as fin, open(output_pdb, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                fout.write(line)
        fout.write("END\n")


def write_fasta_record(handle, header, sequence):
    """Write one FASTA record using 80-character line wrapping."""
    handle.write(f">{header}\n")
    for i in range(0, len(sequence), 80):
        handle.write(sequence[i:i + 80] + "\n")


def get_chain_residues(all_residues, target_chain):
    """Return residues belonging to the requested target chain."""
    chain_residues = [
        residue
        for residue in all_residues
        if residue["chain"] == target_chain
    ]

    if not chain_residues:
        raise ValueError(f"No residues found for chain '{target_chain}'")

    return chain_residues


def get_mutable_residues(chain_residues, positions):
    """Filter chain residues according to selected mutable positions."""
    if positions == "ALL":
        return chain_residues

    selected_positions = {
        position.strip()
        for position in positions.split(",")
        if position.strip()
    }

    mutable_residues = [
        residue
        for residue in chain_residues
        if residue["position"] in selected_positions
    ]

    if not mutable_residues:
        raise ValueError(
            f"No mutable residues matched selected positions: {positions}"
        )

    return mutable_residues


def main():
    """Run mutational screening."""
    args = parse_args()

    input_pdb = Path(args.input_pdb)
    clean_pdb = Path("clean_input_for_pyrosetta.pdb")

    if args.target_chain == "ALL":
        raise ValueError("Please provide a concrete --target_chain, e.g. A")

    clean_pdb_for_pyrosetta(input_pdb, clean_pdb)

    pyrosetta.init("-mute all -ignore_unrecognized_res true")
    pose = pyrosetta.pose_from_pdb(str(clean_pdb))

    scorefxn = pyrosetta.get_fa_scorefxn()
    wt_score = scorefxn(pose)

    all_residues = parse_pdb_residues(input_pdb)
    chain_residues = get_chain_residues(all_residues, args.target_chain)
    mutable_residues = get_mutable_residues(chain_residues, args.positions)

    wt_sequence = "".join(residue["wildtype"] for residue in chain_residues)

    position_to_index = {
        residue["position"]: index
        for index, residue in enumerate(chain_residues)
    }

    with open(args.output, "w", newline="") as out_csv, open(args.fasta_output, "w") as out_fasta:
        writer = csv.writer(out_csv)

        writer.writerow([
            "candidate_id",
            "chain",
            "position",
            "sequence_index",
            "wildtype",
            "mutant",
            "mutation",
            "num_mutations",
            "ddg",
            "wt_score",
            "mutant_score",
            "input_pdb",
            "fasta_id",
            "sequence",
            "wt_sequence",
        ])

        counter = 1

        for residue in mutable_residues:
            wildtype = residue["wildtype"]
            chain = residue["chain"]
            pdb_resnum = int(residue["position"])
            sequence_index = position_to_index[residue["position"]]

            pose_index = pose.pdb_info().pdb2pose(chain, pdb_resnum)

            if pose_index == 0:
                print(
                    f"WARNING: could not map residue {chain}{pdb_resnum} "
                    "to PyRosetta pose index"
                )
                continue

            for mutant in STANDARD_AA:
                if mutant == wildtype:
                    continue

                mutation = f"{wildtype}{chain}{residue['position']}{mutant}"
                candidate_id = f"screening_{counter:06d}"

                mutant_pose = pose.clone()
                mutate_residue(mutant_pose, pose_index, mutant)

                mutant_score = scorefxn(mutant_pose)
                ddg = mutant_score - wt_score

                mutant_sequence = list(wt_sequence)
                mutant_sequence[sequence_index] = mutant
                mutant_sequence = "".join(mutant_sequence)

                fasta_id = f"{candidate_id}|{mutation}"

                writer.writerow([
                    candidate_id,
                    chain,
                    residue["position"],
                    sequence_index + 1,
                    wildtype,
                    mutant,
                    mutation,
                    1,
                    round(ddg, 4),
                    round(wt_score, 4),
                    round(mutant_score, 4),
                    input_pdb.name,
                    fasta_id,
                    mutant_sequence,
                    wt_sequence,
                ])

                write_fasta_record(out_fasta, fasta_id, mutant_sequence)
                counter += 1


if __name__ == "__main__":
    main()