#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

def read_sequence(path_or_seq: str) -> str:
    """Si es una ruta a fasta, lee la secuencia. Si ya es la secuencia, la devuelve."""
    p = Path(path_or_seq)
    if p.exists():
        seq_lines = []
        with open(p) as fh:
            for line in fh:
                line = line.strip()
                if line and not line.startswith(">"):
                    seq_lines.append(line)
        return "".join(seq_lines)
    return path_or_seq  # ya es la secuencia directamente

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates",        required=True)
    parser.add_argument("--ligand-ccd")
    parser.add_argument("--ligand-smiles")
    parser.add_argument("--receptor-sequence")
    parser.add_argument("--runs", type=int, default=3)
    args = parser.parse_args()

    # Determinar modo
    if args.ligand_ccd:
        ligand_type  = "ccd"
        ligand_value = args.ligand_ccd
        mode = "small_molecule"
    elif args.ligand_smiles:
        ligand_type  = "smiles"
        ligand_value = args.ligand_smiles
        mode = "small_molecule"
    elif args.receptor_sequence:
        mode = "protein_protein"
        receptor_seq = read_sequence(args.receptor_sequence)  # ← AQUÍ el fix
    else:
        raise ValueError("Needs --ligand-ccd, --ligand-smiles, or --receptor-sequence.")

    Path("boltz_inputs").mkdir(exist_ok=True)
    manifest_rows = []

    with open(args.candidates, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            candidate_id = row.get("candidate_id") or row.get("id")
            sequence     = row.get("final_sequence") or row.get("binder_sequence") or row.get("sequence")
            if not candidate_id or not sequence:
                continue

            for rep in range(args.runs):
                run_id = f"{candidate_id}_rep{rep}"
                out    = Path("boltz_inputs") / f"{run_id}.fasta"

                if mode == "small_molecule":
                    out.write_text(
                        f">A|protein\n{sequence}\n"
                        f">B|{ligand_type}\n{ligand_value}\n"
                    )
                else:  # protein_protein
                    out.write_text(
                        f">A|protein\n{receptor_seq}\n"  # ← ya es la secuencia
                        f">B|protein\n{sequence}\n"
                    )

                manifest_rows.append({
                    "candidate_id": candidate_id,
                    "replicate":    rep,
                    "run_id":       run_id,
                    "fasta":        str(out),
                })

    with open("boltz_inputs_manifest.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["candidate_id", "replicate", "run_id", "fasta"],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(manifest_rows)

if __name__ == "__main__":
    main()