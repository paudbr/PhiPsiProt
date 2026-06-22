/*
 * Ranking combinado + selección top-N para la ronda final (AF3).
 *
 * Recoge los _metrics.tsv de TODOS los supervivientes del filtro post-Chai
 * (una fila por candidato: id, mean_plddt, ptm, ...), calcula un score
 * combinado normalizado (pesos configurables), aplica un umbral mínimo de
 * score y recorta a top-N. Emite un CSV con la columna in_top_n.
 *
 * Nota: recibe TODOS los metrics y TODOS los fastas supervivientes;
 * decide cuáles van a AF3.
 */
process RANK_AND_SELECT {
    tag "rank_and_select"
    label 'process_single'
    container "nf-core/proteinfold_esmfold:2.0.0"

    input:
    path metrics_tsvs   // muchos *_metrics.tsv (collect)
    val  top_n
    val  w_plddt        // peso pLDDT (p.ej. 0.5)
    val  w_ptm          // peso pTM   (p.ej. 0.5)
    val  score_min      // umbral mínimo de score combinado [0-1]

    output:
    path "ranked.csv", emit: ranked
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    import csv, glob, sys

    rows = []
    for f in sorted(glob.glob("*_metrics.tsv")):
        with open(f) as fh:
            for row in csv.DictReader(fh, delimiter="\\t"):
                try:
                    plddt = float(row["mean_plddt"])
                except (KeyError, ValueError):
                    continue
                try:
                    ptm = float(row["ptm"])
                except (KeyError, ValueError):
                    ptm = None
                rows.append({"id": row["id"], "plddt": plddt, "ptm": ptm})

    if not rows:
        raise SystemExit("No hay métricas para rankear.")

    # Normalización min-max por métrica (robusta a que ptm falte)
    def norm(vals):
        v = [x for x in vals if x is not None]
        if not v:
            return lambda x: 0.0
        lo, hi = min(v), max(v)
        if hi - lo < 1e-9:
            return lambda x: 1.0 if x is not None else 0.0
        return lambda x: (x - lo) / (hi - lo) if x is not None else 0.0

    # pLDDT viene en escala 0-100; lo normalizamos igual por consistencia
    n_plddt = norm([r["plddt"] for r in rows])
    n_ptm   = norm([r["ptm"]   for r in rows])

    wp, wt = float(${w_plddt}), float(${w_ptm})
    wsum = wp + wt if (wp + wt) > 0 else 1.0

    for r in rows:
        r["score"] = (wp * n_plddt(r["plddt"]) + wt * n_ptm(r["ptm"])) / wsum

    rows.sort(key=lambda r: r["score"], reverse=True)

    top_n = int(${top_n})
    smin  = float(${score_min})
    for i, r in enumerate(rows):
        in_top = (i < top_n) and (r["score"] >= smin)
        r["rank"] = i + 1
        r["in_top_n"] = "yes" if in_top else "no"

    with open("ranked.csv", "w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["rank", "id", "mean_plddt", "ptm", "score", "in_top_n"])
        for r in rows:
            w.writerow([
                r["rank"], r["id"], f"{r['plddt']:.3f}",
                f"{r['ptm']:.4f}" if r["ptm"] is not None else "NA",
                f"{r['score']:.4f}", r["in_top_n"],
            ])

    with open("versions.yml", "w") as v:
        v.write('"${task.process}":\\n')
        v.write(f"    python: {sys.version.split()[0]}\\n")
    print(f"Rankeados {len(rows)}; top_n={top_n}, score_min={smin}")
    """
}
