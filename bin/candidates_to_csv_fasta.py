#!/usr/bin/env python3
import sys
import pandas as pd

def main():
    csv_file = sys.argv[1]

    df = pd.read_csv(csv_file)

    for _, row in df.iterrows():
        cid = row["candidate_id"]
        seq = row["binder_sequence"]

        with open(f"{cid}.fasta", "w") as f:
            f.write(f">{cid}\n")
            f.write(seq + "\n")

if __name__ == "__main__":
    main()