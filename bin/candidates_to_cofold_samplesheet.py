#!/usr/bin/env python3
import argparse
import csv
import os

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--summary",    required=True)
    p.add_argument("--candidates", required=True)
    p.add_argument("--receptor",   required=True)
    p.add_argument("--outdir",     default="fastas")
    return p.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Leer secuencias del CSV original
    sequences = {}
    with open(args.candidates) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            sequences[row["candidate_id"]] = row["binder_sequence"]

    # Leer summary rankeado por affinity
    ranked = []
    with open(args.summary) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                aff = float(row["affinity_pred_value_mean"]) if row["affinity_pred_value_mean"] else float("inf")
            except ValueError:
                aff = float("inf")
            ranked.append((aff, row["candidate_id"]))

    # Ordenar: menor affinity_pred_value = mejor binder
    ranked.sort(key=lambda x: x[0])

    # Generar fastas y samplesheet
    samplesheet_rows = []
    for _, candidate_id in ranked:
        if candidate_id not in sequences:
            print(f"WARNING: {candidate_id} not found in candidates CSV, skipping")
            continue

        binder_seq   = sequences[candidate_id]
        receptor_seq = args.receptor

        fasta_path = os.path.join(args.outdir, f"{candidate_id}.fa")
        with open(fasta_path, "w") as fh:
            fh.write(f">receptor\n{receptor_seq}\n")
            fh.write(f">binder\n{binder_seq}\n")

        samplesheet_rows.append({
            "id":    candidate_id,
            "fasta": os.path.abspath(fasta_path)
        })

    with open("cofold_samplesheet.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "fasta"])
        writer.writeheader()
        writer.writerows(samplesheet_rows)

    print(f"Generated {len(samplesheet_rows)} cofold entries")

if __name__ == "__main__":
    main()