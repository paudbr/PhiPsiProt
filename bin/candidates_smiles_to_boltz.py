#!/usr/bin/env python3
"""Genera un YAML de Boltz2 (proteina + ligando SMILES + affinity + pocket)
sin dependencias externas: escribe el YAML a mano.

Soporta homo-oligomeros: con --protein-chains "A,B" escribe la misma
secuencia en varias cadenas (p. ej. la proteasa del HIV como homodimero).
"""
import argparse
import sys


def parse_args():
    p = argparse.ArgumentParser(description="SMILES -> Boltz2 YAML")
    p.add_argument("--candidate-id", required=True)
    p.add_argument("--smiles", required=True)
    p.add_argument("--receptor-sequence", required=True)
    p.add_argument("--protein-id", default="A",
                   help="ID de cadena por defecto si no se usa --protein-chains.")
    p.add_argument("--protein-chains", default=None,
                   help="IDs de cadena de proteina separados por coma, ej: 'A,B' "
                        "para homodimero. Todas comparten la misma secuencia. "
                        "Si no se da, usa --protein-id (una sola cadena).")
    p.add_argument("--ligand-id", default="L")
    p.add_argument("--output", required=True)
    p.add_argument("--no-affinity", action="store_true")
    p.add_argument("--pocket-residues", default="",
                   help="Lista de residuos del pocket, ej: 'A8,A23,B27'. "
                        "Si no se da, no se escribe el bloque de constraints.")
    p.add_argument("--pocket-max-distance", type=float, default=None,
                   help="max_distance opcional para la pocket constraint (Angstrom).")
    return p.parse_args()


def parse_pocket_residues(raw, default_chain):
    """Convierte 'A8,A23,B27' o '8,23,27' en [(chain, resid), ...].

    - Si el token empieza por letra, esa letra es la cadena.
    - Si es solo un numero, se usa default_chain como cadena.
    """
    contacts = []
    for tok in raw.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if tok[0].isalpha():
            chain = tok[0]
            resid = tok[1:]
        else:
            chain = default_chain
            resid = tok
        if not resid.isdigit():
            sys.exit(f"ERROR: residuo de pocket no valido: '{tok}'")
        contacts.append((chain, int(resid)))
    return contacts


def main():
    args = parse_args()
    seq = args.receptor_sequence.strip().upper()
    smiles = args.smiles.strip()

    if not seq:
        sys.exit("ERROR: receptor sequence vacia")
    if not smiles:
        sys.exit("ERROR: SMILES vacio")

    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    bad = set(seq) - valid_aa
    if bad:
        sys.exit(f"ERROR: secuencia con caracteres no validos: {sorted(bad)}")

    # cadenas de proteina (homo-oligomero: misma secuencia, distinto id)
    if args.protein_chains:
        chain_ids = [c.strip() for c in args.protein_chains.split(",") if c.strip()]
    else:
        chain_ids = [args.protein_id]
    if not chain_ids:
        sys.exit("ERROR: no se definieron cadenas de proteina")
    if len(chain_ids) != len(set(chain_ids)):
        sys.exit(f"ERROR: ids de cadena duplicados en --protein-chains: {chain_ids}")

    default_chain = chain_ids[0]
    contacts = parse_pocket_residues(args.pocket_residues, default_chain)

    # validacion: toda cadena referenciada en los contactos debe existir
    used_chains = {c for c, _ in contacts}
    missing = used_chains - set(chain_ids)
    if missing:
        sys.exit(
            f"ERROR: los contactos del pocket usan cadenas {sorted(missing)} "
            f"pero solo se declararon {chain_ids}. "
            f"Pasa --protein-chains '{','.join(sorted(set(chain_ids) | used_chains))}'."
        )

    # ── construir YAML ────────────────────────────────────────────
    lines = ["version: 1", "sequences:"]
    for cid in chain_ids:
        lines += [
            "  - protein:",
            f"      id: {cid}",
            f"      sequence: {seq}",
        ]
    lines += [
        "  - ligand:",
        f"      id: {args.ligand_id}",
        f'      smiles: "{smiles}"',
    ]

    # bloque de constraints (pocket) — solo si hay residuos
    if contacts:
        lines += [
            "constraints:",
            "  - pocket:",
            f"      binder: {args.ligand_id}",
            "      contacts:",
        ]
        for chain, resid in contacts:
            lines.append(f"        - [{chain}, {resid}]")
        if args.pocket_max_distance is not None:
            lines.append(f"      max_distance: {args.pocket_max_distance}")

    if not args.no_affinity:
        lines += [
            "properties:",
            "  - affinity:",
            f"      binder: {args.ligand_id}",
        ]

    with open(args.output, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    print(f"[candidates_smiles_to_boltz] escrito {args.output} "
          f"para candidato {args.candidate_id} "
          f"({len(chain_ids)} cadena(s) de proteina {chain_ids}, "
          f"{len(contacts)} contactos de pocket)", file=sys.stderr)


if __name__ == "__main__":
    main()