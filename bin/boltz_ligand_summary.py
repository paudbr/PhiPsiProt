#!/usr/bin/env python3

import csv
import glob
import json
import os
import re


def parse_run_id(path):
    basename = os.path.basename(path)
    match = re.search(r"confidence_(.+)_model_(\d+)\.json$", basename)
    if match:
        return match.group(1), int(match.group(2))
    parent = os.path.basename(os.path.dirname(path))
    return parent, 0


def main():
    rows = []
    for path in sorted(glob.glob("confidence_jsons/**/*.json", recursive=True)):
        run_id, model = parse_run_id(path)
        candidate_id = re.sub(r"_rep\d+$", "", run_id)
        rep_match = re.search(r"_rep(\d+)$", run_id)
        replicate = int(rep_match.group(1)) if rep_match else 0
        with open(path) as handle:
            data = json.load(handle)
        rows.append({
            "candidate_id": candidate_id,
            "replicate": replicate,
            "model": model,
            "probability_binary": data.get("probability_binary", data.get("complex_plddt", "")),
            "affinity_pred_value": data.get("affinity_pred_value", data.get("affinity", "")),
            "ptm": data.get("ptm", ""),
            "iptm": data.get("iptm", ""),
        })

    with open("boltz2_ligand_scores.tsv", "w", newline="") as handle:
        fields = ["candidate_id", "replicate", "model", "probability_binary", "affinity_pred_value", "ptm", "iptm"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    grouped = {}
    for row in rows:
        grouped.setdefault(row["candidate_id"], []).append(row)

    with open("boltz2_ligand_summary.tsv", "w", newline="") as handle:
        fields = ["candidate_id", "n", "probability_binary_mean", "affinity_pred_value_mean"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for candidate_id, vals in grouped.items():
            out = {"candidate_id": candidate_id, "n": len(vals)}
            for metric in ["probability_binary", "affinity_pred_value"]:
                xs = []
                for val in vals:
                    try:
                        xs.append(float(val[metric]))
                    except (TypeError, ValueError):
                        pass
                out[f"{metric}_mean"] = sum(xs) / len(xs) if xs else ""
            writer.writerow(out)


if __name__ == "__main__":
    main()
