#!/usr/bin/env python3
"""
Filtra candidatos con secuencias de binder degeneradas (baja complejidad).
Detecta: tractos homopoliméricos, baja diversidad de aminoácidos, y
sesgo extremo de composición — típicos de diseños fallidos.
"""
import argparse
import re
from collections import Counter
import math
import pandas as pd


def max_homopolymer_run(seq):
    """Longitud de la tirada más larga del mismo aminoácido (p.ej. AAAAA -> 5)."""
    if not seq:
        return 0
    runs = [len(m.group()) for m in re.finditer(r'(.)\1*', seq)]
    return max(runs)


def shannon_entropy(seq):
    """Entropía de Shannon de la composición (bits). Baja = poco diversa."""
    if not seq:
        return 0.0
    counts = Counter(seq)
    n = len(seq)
    return -sum((c / n) * math.log2(c / n) for c in counts.values())


def max_aa_fraction(seq):
    """Fracción del aminoácido más frecuente (1.0 = todo el mismo)."""
    if not seq:
        return 1.0
    counts = Counter(seq)
    return max(counts.values()) / len(seq)


def is_degenerate(seq, max_run=5, min_entropy=2.5, max_single_aa=0.40):
    """
    Devuelve (es_degenerada, motivos).
    Umbrales ajustables:
      - max_run: tirada homopolimérica máxima permitida (5 = descarta AAAAAA)
      - min_entropy: entropía mínima de composición en bits
      - max_single_aa: fracción máxima de un solo aminoácido
    """
    reasons = []
    run = max_homopolymer_run(seq)
    if run > max_run:
        reasons.append(f"homopolymer_run={run}(>{max_run})")

    ent = shannon_entropy(seq)
    if ent < min_entropy:
        reasons.append(f"entropy={ent:.2f}(<{min_entropy})")

    frac = max_aa_fraction(seq)
    if frac > max_single_aa:
        reasons.append(f"max_aa_frac={frac:.2f}(>{max_single_aa})")

    return (len(reasons) > 0, reasons)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="CSV de candidatos")
    ap.add_argument("--output", required=True, help="CSV filtrado (solo los que pasan)")
    ap.add_argument("--rejected", default=None, help="CSV opcional con los descartados y motivos")
    ap.add_argument("--seq-col", default="fasta")
    ap.add_argument("--id-col", default="candidate_id")
    ap.add_argument("--max-run", type=int, default=5)
    ap.add_argument("--min-entropy", type=float, default=2.5)
    ap.add_argument("--max-single-aa", type=float, default=0.40)
    args = ap.parse_args()

    df = pd.read_csv(args.input)

    if args.seq_col not in df.columns:
        raise SystemExit(f"No existe la columna '{args.seq_col}'. Columnas: {list(df.columns)}")

    flags = df[args.seq_col].astype(str).apply(
        lambda s: is_degenerate(s, args.max_run, args.min_entropy, args.max_single_aa)
    )
    df["_degenerate"] = flags.apply(lambda t: t[0])
    df["_reasons"] = flags.apply(lambda t: ";".join(t[1]))

    kept = df[~df["_degenerate"]].drop(columns=["_degenerate", "_reasons"])
    rejected = df[df["_degenerate"]]

    kept.to_csv(args.output, index=False)
    print(f"Total: {len(df)} | Conservados: {len(kept)} | Descartados: {len(rejected)}")

    if len(rejected):
        print("\nDescartados:")
        for _, row in rejected.iterrows():
            print(f"  {row[args.id_col]}: {row['_reasons']}  seq={row[args.seq_col]}")

    if args.rejected and len(rejected):
        rejected.to_csv(args.rejected, index=False)


if __name__ == "__main__":
    main()