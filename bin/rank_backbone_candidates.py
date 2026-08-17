#!/usr/bin/env python3

import argparse
import csv


LOWER_IS_BETTER = {
    "mpnn_score",
    "mpnn_global_score",
    "instability_index",
    "gravy",
}

HIGHER_IS_BETTER = {
    "solubility_score",
    "mpnn_seq_recovery",
}


def to_float(value):
    try:
        return float(value)
    except Exception:
        return None


def parse_weights(text):
    weights = {}
    for item in text.split(","):
        if not item.strip():
            continue
        name, value = item.split(":")
        weights[name.strip()] = float(value)

    total = sum(weights.values())
    if total <= 0:
        raise ValueError("Ranking weights must sum to > 0")

    return {k: v / total for k, v in weights.items()}


def normalize(values, lower_is_better=True):
    numeric = [v for v in values if v is not None]

    if not numeric:
        return [0.0 for _ in values]

    min_v = min(numeric)
    max_v = max(numeric)

    if min_v == max_v:
        return [1.0 if v is not None else 0.0 for v in values]

    normalized = []

    for v in values:
        if v is None:
            normalized.append(0.0)
            continue

        if lower_is_better:
            score = (max_v - v) / (max_v - min_v)
        else:
            score = (v - min_v) / (max_v - min_v)

        normalized.append(score)

    return normalized


def write_fasta(rows, output_fasta):
    with open(output_fasta, "w") as out:
        for row in rows:
            seq = (
                row.get("final_sequence")
                or row.get("designed_sequence")
                or row.get("biophysical_sequence")
                or row.get("sequence")
                or ""
            ).replace(" ", "").replace("\n", "")

            if not seq:
                continue

            name = row.get("design_id") or row.get("candidate_id") or "candidate"
            out.write(f">{name}\n")

            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_csv", required=True)
    parser.add_argument("--output_csv", required=True)
    parser.add_argument("--output_fasta", required=True)
    parser.add_argument(
        "--weights",
        default="mpnn_score:0.4,instability_index:0.4,solubility_score:0.2",
    )
    parser.add_argument("--top_n", type=int, default=20)

    args = parser.parse_args()

    weights = parse_weights(args.weights)

    with open(args.input_csv, newline="") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise SystemExit("Input CSV is empty")

    for metric in weights:
        if metric not in rows[0]:
            raise ValueError(f"Metric '{metric}' not found in input CSV")

    metric_norms = {}

    for metric in weights:
        values = [to_float(row.get(metric)) for row in rows]

        if metric in LOWER_IS_BETTER:
            metric_norms[metric] = normalize(values, lower_is_better=True)
        elif metric in HIGHER_IS_BETTER:
            metric_norms[metric] = normalize(values, lower_is_better=False)
        else:
            raise ValueError(
                f"Metric '{metric}' has no direction. "
                f"Add it to LOWER_IS_BETTER or HIGHER_IS_BETTER."
            )

    for i, row in enumerate(rows):
        final_score = 0.0

        for metric, weight in weights.items():
            norm_value = metric_norms[metric][i]
            row[f"{metric}_norm"] = round(norm_value, 4)
            final_score += weight * norm_value

        row["ranking_score"] = round(final_score, 4)
        row["ranking_weights"] = args.weights

    ranked = sorted(
        rows,
        key=lambda r: to_float(r.get("ranking_score")) or -1,
        reverse=True,
    )

    top_rows = ranked[:args.top_n]

    fieldnames = list(top_rows[0].keys())

    priority = [
        "rank",
        "ranking_score",
        "design_id",
        "candidate_id",
        "mpnn_score",
        "mpnn_score_norm",
        "instability_index",
        "instability_index_norm",
        "solubility_score",
        "solubility_score_norm",
        "gravy",
        "isoelectric_point",
        "mpnn_seq_recovery",
        "pdb",
        "final_sequence",
    ]

    ordered = []
    for col in priority:
        if col in fieldnames or col == "rank":
            ordered.append(col)

    for col in fieldnames:
        if col not in ordered:
            ordered.append(col)

    for idx, row in enumerate(top_rows, start=1):
        row["rank"] = idx

    with open(args.output_csv, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=ordered)
        writer.writeheader()
        writer.writerows(top_rows)

    write_fasta(top_rows, args.output_fasta)

    print(f"Wrote {len(top_rows)} ranked candidates to {args.output_csv}")
    print(f"Wrote FASTA to {args.output_fasta}")


if __name__ == "__main__":
    main()
