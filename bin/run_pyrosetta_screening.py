#!/usr/bin/env python3

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


def write_fasta_record(handle, header, sequence):
    handle.write(f">{header}\n")
    for i in range(0, len(sequence), 80):
        handle.write(sequence[i:i + 80] + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_pdb", required=True)
    parser.add_argument("--target_chain", default="ALL")
    parser.add_argument("--positions", default="ALL")
    parser.add_argument("--output", required=True)
    parser.add_argument("--fasta_output", required=True)
    args = parser.parse_args()

    pdb = Path(args.input_pdb)

    clean_pdb = Path("clean_input_for_pyrosetta.pdb")
    with open(pdb) as fin, open(clean_pdb, "w") as fout:
        for line in fin:
            if line.startswith("ATOM"):
                fout.write(line)
        fout.write("END\n")

    pyrosetta.init("-mute all -ignore_unrecognized_res true")
    pose = pyrosetta.pose_from_pdb(str(clean_pdb))


    scorefxn = pyrosetta.get_fa_scorefxn()

    wt_score = scorefxn(pose)

    all_residues = parse_pdb_residues(pdb)

    if args.target_chain == "ALL":
        raise ValueError("Please provide --target_chain, e.g. --target_chain A")

    chain_residues = [r for r in all_residues if r["chain"] == args.target_chain]

    if not chain_residues:
        raise ValueError(f"No residues found for chain {args.target_chain}")

    mutable_residues = chain_residues

    if args.positions != "ALL":
        selected_positions = set(p.strip() for p in args.positions.split(","))
        mutable_residues = [
            r for r in chain_residues
            if r["position"] in selected_positions
        ]

    wt_sequence = "".join(r["wildtype"] for r in chain_residues)

    position_to_index = {
        r["position"]: idx
        for idx, r in enumerate(chain_residues)
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
            "wt_sequence"
        ])

        counter = 1

        for res in mutable_residues:
            wt = res["wildtype"]
            chain = res["chain"]
            pdb_resnum = int(res["position"])
            seq_index = position_to_index[res["position"]]

            pose_index = pose.pdb_info().pdb2pose(chain, pdb_resnum)

            if pose_index == 0:
                print(f"WARNING: could not map {chain}{pdb_resnum} to pose index")
                continue

            for mut in STANDARD_AA:
                if mut == wt:
                    continue

                mutation = f"{wt}{chain}{res['position']}{mut}"
                candidate_id = f"screening_{counter:06d}"

                mutant_pose = pose.clone()
                mutate_residue(mutant_pose, pose_index, mut)

                mutant_score = scorefxn(mutant_pose)
                ddg = mutant_score - wt_score

                mutant_sequence = list(wt_sequence)
                mutant_sequence[seq_index] = mut
                mutant_sequence = "".join(mutant_sequence)

                fasta_id = f"{candidate_id}|{mutation}"

                writer.writerow([
                    candidate_id,
                    chain,
                    res["position"],
                    seq_index + 1,
                    wt,
                    mut,
                    mutation,
                    1,
                    round(ddg, 4),
                    round(wt_score, 4),
                    round(mutant_score, 4),
                    pdb.name,
                    fasta_id,
                    mutant_sequence,
                    wt_sequence,
                ])

                write_fasta_record(out_fasta, fasta_id, mutant_sequence)
                counter += 1


if __name__ == "__main__":
    main()
