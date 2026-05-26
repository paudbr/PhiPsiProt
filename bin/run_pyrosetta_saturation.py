#!/usr/bin/env python3

import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_pdb", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    pdb = Path(args.input_pdb)

    if not pdb.exists():
        raise FileNotFoundError(f"Input PDB not found: {pdb}")

    with open(args.output, "w") as out:
        out.write("candidate_id,mutation,ddg,input_pdb\n")
        out.write(f"test_candidate,NA,0.0,{pdb.name}\n")

if __name__ == "__main__":
    main()
