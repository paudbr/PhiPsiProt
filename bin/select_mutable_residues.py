#!/usr/bin/env python3

"""
Select mutable residues for PhiPsiProt Mode 1.

Purpose
-------
Parse a PDB structure and select residue positions from a target chain
according to the requested residue selection strategy.

Supported modes
---------------
manual:
    User-defined residue positions.

all:
    All residues belonging to the target chain.

pocket:
    Target-chain residues located within a distance cutoff of a ligand.

interface:
    Target-chain residues located within a distance cutoff of another
    protein chain.

Notes
-----
- auto_pocket is handled upstream by the AUTO_POCKET_DETECTION Nextflow module.
- The output is passed to the PyRosetta mutational screening step.

Output
------
selected_positions.txt:
    Comma-separated list of selected residue numbers.
"""

import argparse
import math


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Select mutable residues for PhiPsiProt Mode 1 screening."
    )

    parser.add_argument("--input_pdb", required=True)
    parser.add_argument(
        "--selection_mode",
        default="manual",
        choices=["manual", "all", "pocket", "interface"],
    )
    parser.add_argument("--target_chain", required=True)
    parser.add_argument("--positions", default="ALL")
    parser.add_argument("--ligand_resname", default="")
    parser.add_argument("--interface_chain", default="")
    parser.add_argument("--distance_cutoff", type=float, default=6.0)
    parser.add_argument("--output", required=True)

    return parser.parse_args()


def distance(coord_a, coord_b):
    """Calculate Euclidean distance between two atomic coordinates."""
    return math.sqrt(
        (coord_a[0] - coord_b[0]) ** 2
        + (coord_a[1] - coord_b[1]) ** 2
        + (coord_a[2] - coord_b[2]) ** 2
    )


def parse_pdb(pdb_file):
    """Parse protein atoms and ligand atoms from a PDB file."""
    protein_atoms = []
    ligand_atoms = []

    with open(pdb_file) as handle:
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

    return protein_atoms, ligand_atoms


def sort_residue_ids(residue_ids):
    """Sort residue identifiers numerically."""
    return sorted(residue_ids, key=lambda x: int(x))


def get_target_positions(protein_atoms, target_chain):
    """Collect residue identifiers belonging to the target chain."""
    positions = {
        atom["resnum"]
        for atom in protein_atoms
        if atom["chain"] == target_chain
    }

    return sort_residue_ids(positions)


def select_all_positions(target_positions):
    """Select all residues from the target chain."""
    return set(target_positions)


def select_manual_positions(positions, target_positions):
    """Select user-defined residue positions."""
    if positions == "ALL":
        return set(target_positions)

    return set(
        position.strip()
        for position in positions.split(",")
        if position.strip()
    )


def select_pocket_positions(
    protein_atoms,
    ligand_atoms,
    target_chain,
    ligand_resname,
    distance_cutoff,
):
    """
    Select target-chain residues having at least one atom within the
    specified cutoff distance of the ligand.
    """
    if not ligand_resname:
        raise ValueError("--ligand_resname is required for selection_mode pocket")

    ligand_selected = [
        atom
        for atom in ligand_atoms
        if atom["resname"] == ligand_resname
    ]

    if not ligand_selected:
        raise ValueError(f"No ligand atoms found with resname '{ligand_resname}'")

    selected = set()

    for protein_atom in protein_atoms:
        if protein_atom["chain"] != target_chain:
            continue

        for ligand_atom in ligand_selected:
            if distance(protein_atom["coords"], ligand_atom["coords"]) <= distance_cutoff:
                selected.add(protein_atom["resnum"])
                break

    return selected


def select_interface_positions(
    protein_atoms,
    target_chain,
    interface_chain,
    distance_cutoff,
):
    """Select target-chain residues located near a partner protein chain."""
    if not interface_chain:
        raise ValueError("--interface_chain is required for selection_mode interface")

    interface_atoms = [
        atom
        for atom in protein_atoms
        if atom["chain"] == interface_chain
    ]

    if not interface_atoms:
        raise ValueError(f"No atoms found for interface chain '{interface_chain}'")

    selected = set()

    for protein_atom in protein_atoms:
        if protein_atom["chain"] != target_chain:
            continue

        for interface_atom in interface_atoms:
            if distance(protein_atom["coords"], interface_atom["coords"]) <= distance_cutoff:
                selected.add(protein_atom["resnum"])
                break

    return selected


def write_selected_positions(selected_positions, output_file):
    """Write selected residue positions as a single comma-separated line."""
    with open(output_file, "w") as out:
        out.write(",".join(selected_positions))
        out.write("\n")


def main():
    """Run residue selection."""
    args = parse_args()

    protein_atoms, ligand_atoms = parse_pdb(args.input_pdb)

    target_positions = get_target_positions(
        protein_atoms,
        args.target_chain,
    )

    if not target_positions:
        raise ValueError(
            f"No residues found for target chain '{args.target_chain}' "
            f"in {args.input_pdb}"
        )

    if args.selection_mode == "all":
        selected = select_all_positions(target_positions)

    elif args.selection_mode == "manual":
        selected = select_manual_positions(
            args.positions,
            target_positions,
        )

    elif args.selection_mode == "pocket":
        selected = select_pocket_positions(
            protein_atoms,
            ligand_atoms,
            args.target_chain,
            args.ligand_resname,
            args.distance_cutoff,
        )

    elif args.selection_mode == "interface":
        selected = select_interface_positions(
            protein_atoms,
            args.target_chain,
            args.interface_chain,
            args.distance_cutoff,
        )

    else:
        raise ValueError(f"Unsupported selection mode: {args.selection_mode}")

    selected = sort_residue_ids(selected)

    if not selected:
        raise ValueError(
            f"No residues selected using mode '{args.selection_mode}'. "
            "Check target_chain, ligand/interface settings and distance_cutoff."
        )

    write_selected_positions(selected, args.output)


if __name__ == "__main__":
    main()