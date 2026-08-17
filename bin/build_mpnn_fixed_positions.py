#!/usr/bin/env python3

import argparse
import json
import re
from pathlib import Path


def fixed_positions_from_contig(contig):
    contig = contig.strip().replace("[", "").replace("]", "").replace("'", "").replace('"', "")
    parts = [p.strip() for p in contig.split("/") if p.strip()]

    fixed = []
    output_pos = 1

    for part in parts:
        motif = re.match(r"^[A-Za-z](\d+)-(\d+)$", part)
        generated = re.match(r"^(\d+)-(\d+)$", part)

        if motif:
            start = int(motif.group(1))
            end = int(motif.group(2))
            length = abs(end - start) + 1

            fixed.extend(range(output_pos, output_pos + length))
            output_pos += length

        elif generated:
            length = int(generated.group(2))
            output_pos += length

    return fixed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", required=True)
    parser.add_argument("--contig", required=True)
    parser.add_argument("--chain", default="A")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    pdb_name = Path(args.pdb).stem
    fixed = fixed_positions_from_contig(args.contig)

    data = {
        pdb_name: {
            args.chain: fixed
        }
    }

    with open(args.output, "w") as out:
        out.write(json.dumps(data) + "\n")

    print(f"Fixed {len(fixed)} positions for ProteinMPNN")
    print(json.dumps(data))


if __name__ == "__main__":
    main()
