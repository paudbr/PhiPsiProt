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


def parse_pocket_residues(pocket_str: str):
    """
    "A5,A4,A10,A238" -> [("A", 5), ("A", 4), ("A", 10), ("A", 238)]
    """
    contacts = []
    if not pocket_str:
        return contacts
    for tok in pocket_str.split(","):
        tok = tok.strip()
        if not tok:
            continue
        i = 0
        while i < len(tok) and tok[i].isalpha():
            i += 1
        chain = tok[:i]
        resnum = tok[i:]
        if not chain or not resnum.isdigit():
            raise ValueError(
                f"Residuo de pocket mal formado: '{tok}'. "
                f"Esperado formato tipo 'A5' (cadena + numero)."
            )
        contacts.append((chain, int(resnum)))
    return contacts


def yaml_bool(b: bool) -> str:
    return "true" if b else "false"


def build_yaml_text(receptor_seq, binder_seq, pocket_contacts,
                    max_distance, force_pocket,
                    binder_no_msa=False, with_affinity=False):
    """
    YAML de Boltz-2 (proteina-proteina), escrito a mano (sin PyYAML).
    Cadena A = receptor, Cadena B = binder.
    - binder_no_msa: si True, anade 'msa: empty' a la cadena B (single-sequence).
    - with_affinity: si True, anade properties: affinity: binder B.
    """
    lines = []
    lines.append("version: 1")
    lines.append("sequences:")
    lines.append("  - protein:")
    lines.append("      id: A")
    lines.append(f"      sequence: {receptor_seq}")
    lines.append("  - protein:")
    lines.append("      id: B")
    lines.append(f"      sequence: {binder_seq}")
    if binder_no_msa:
        lines.append("      msa: empty")

    if pocket_contacts:
        contacts_str = ", ".join(
            f"[{chain}, {resnum}]" for chain, resnum in pocket_contacts
        )
        lines.append("constraints:")
        lines.append("  - pocket:")
        lines.append("      binder: B")
        lines.append(f"      contacts: [{contacts_str}]")
        lines.append(f"      max_distance: {max_distance}")
        lines.append(f"      force: {yaml_bool(force_pocket)}")

    if with_affinity:
        lines.append("properties:")
        lines.append("  - affinity:")
        lines.append("      binder: B")

    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates",        required=True)
    parser.add_argument("--ligand-ccd")
    parser.add_argument("--ligand-smiles")
    parser.add_argument("--receptor-sequence")
    parser.add_argument("--pocket-residues", default="",
                        help="Residuos del pocket, ej: 'A5,A4,A10'")
    parser.add_argument("--max-distance", type=float, default=6.0)
    parser.add_argument("--force-pocket", action="store_true")
    parser.add_argument("--binder_no_msa", action="store_true",
                        help="Cadena B (binder) en single-sequence: msa: empty")
    parser.add_argument("--with_affinity", action="store_true",
                        help="Anade properties: affinity: binder B")
    parser.add_argument("--runs", type=int, default=3)
    args = parser.parse_args()

    if not args.receptor_sequence:
        raise ValueError(
            "Este generador YAML espera --receptor-sequence (modo proteina-proteina)."
        )
    receptor_seq = read_sequence(args.receptor_sequence)
    pocket_contacts = parse_pocket_residues(args.pocket_residues)

    Path("boltz_inputs").mkdir(exist_ok=True)
    manifest_rows = []

    with open(args.candidates, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            candidate_id = row.get("candidate_id") or row.get("id")
            sequence = (row.get("final_sequence")
                        or row.get("binder_sequence")
                        or row.get("sequence"))
            if not candidate_id or not sequence:
                continue

            for rep in range(args.runs):
                run_id = f"{candidate_id}_rep{rep}"
                out = Path("boltz_inputs") / f"{run_id}.yaml"

                yaml_text = build_yaml_text(
                    receptor_seq=receptor_seq,
                    binder_seq=sequence,
                    pocket_contacts=pocket_contacts,
                    max_distance=args.max_distance,
                    force_pocket=args.force_pocket,
                    binder_no_msa=args.binder_no_msa,
                    with_affinity=args.with_affinity,
                )
                out.write_text(yaml_text)

                manifest_rows.append({
                    "candidate_id": candidate_id,
                    "replicate": rep,
                    "run_id": run_id,
                    "yaml": str(out),
                })

    with open("boltz_inputs_manifest.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["candidate_id", "replicate", "run_id", "yaml"],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(manifest_rows)


if __name__ == "__main__":
    main()
