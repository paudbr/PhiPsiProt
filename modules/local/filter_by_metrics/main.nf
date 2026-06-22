/*
 * Filtro por métricas de plegamiento. Sirve para los dos cortes de la cascada:
 *
 *   - Tras ESMFold:  pLDDT medio + scRMSD vs original (pasa ptm_tsv = NO_FILE)
 *   - Tras Chai-1:   pLDDT medio + pTM  (pasa scrmsd_tsv = NO_FILE)
 *
 * Entradas (tuple): meta, plddt_tsv, ptm_tsv, scrmsd_tsv, fasta
 *
 * IMPORTANTE - el TSV de pLDDT tiene cabecera DISTINTA según la herramienta:
 *   - ESMFold (extract_metrics.py): columnas  'Positions'  'rank_0'   <- pLDDT en rank_0
 *   - Chai-1:                       columnas  'chain' 'residue' 'plddt'
 * Por eso la columna de pLDDT se autodetecta:
 *   1) si existe 'plddt' -> esa
 *   2) si existe 'rank_0' -> esa
 *   3) si no, la ÚLTIMA columna que sea numérica
 *
 * Umbrales: plddt_min, ptm_min, scrmsd_max  (<=0 desactiva el de pTM / scRMSD).
 * El candidato que NO pasa borra su fasta -> no aparece en el canal 'passed'.
 */
process FILTER_BY_METRICS {
    tag "$meta.id"
    label 'process_single'
    container "nf-core/proteinfold_esmfold:2.0.0"

    input:
    tuple val(meta), path(plddt_tsv), path(ptm_tsv), path(scrmsd_tsv), path(fasta)
    val   plddt_min
    val   ptm_min
    val   scrmsd_max

    output:
    tuple val(meta), path(fasta)                     , emit: passed, optional: true
    tuple val(meta), path("${meta.id}_metrics.tsv")  , emit: report
    path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    import csv, os, sys

    def is_real(p):
        return os.path.basename(p) not in ("NO_FILE", "NO_FILE_PAE") and os.path.exists(p)

    # --- pLDDT medio: autodetecta la columna segun la herramienta ---
    plddt_vals = []
    with open("${plddt_tsv}") as fh:
        reader = csv.DictReader(fh, delimiter="\\t")
        fields = reader.fieldnames or []

        # elegir columna de pLDDT
        plddt_col = None
        for cand in ("plddt", "rank_0"):
            if cand in fields:
                plddt_col = cand
                break
        if plddt_col is None:
            # fallback: ultima columna que parsee como float en la 1a fila
            rows = list(reader)
            for col in reversed(fields):
                try:
                    float(rows[0][col]); plddt_col = col; break
                except (ValueError, KeyError, IndexError):
                    continue
            data_rows = rows
        else:
            data_rows = reader

        if plddt_col is None:
            sys.stderr.write(f"[WARN] no encuentro columna pLDDT en {fields}\\n")
        else:
            for row in data_rows:
                try:
                    plddt_vals.append(float(row[plddt_col]))
                except (KeyError, ValueError, TypeError):
                    pass

    mean_plddt = sum(plddt_vals) / len(plddt_vals) if plddt_vals else 0.0

    # --- pTM (opcional; modelo top = primera fila de datos) ---
    ptm = None
    if is_real("${ptm_tsv}"):
        with open("${ptm_tsv}") as fh:
            for row in csv.DictReader(fh, delimiter="\\t"):
                try:
                    ptm = float(row["ptm"])
                except (KeyError, ValueError, TypeError):
                    ptm = None
                break

    # --- scRMSD (opcional) ---
    scrmsd = None
    if is_real("${scrmsd_tsv}"):
        with open("${scrmsd_tsv}") as fh:
            for row in csv.DictReader(fh, delimiter="\\t"):
                try:
                    scrmsd = float(row["scrmsd"])
                except (KeyError, ValueError, TypeError):
                    scrmsd = None
                break

    plddt_ok  = mean_plddt >= float(${plddt_min})
    ptm_ok    = True if (${ptm_min} <= 0 or ptm is None) else (ptm >= float(${ptm_min}))
    if ${scrmsd_max} <= 0:
        scrmsd_ok = True
    elif scrmsd is None:
        scrmsd_ok = False
    else:
        scrmsd_ok = scrmsd <= float(${scrmsd_max})

    passed = plddt_ok and ptm_ok and scrmsd_ok

    with open("${meta.id}_metrics.tsv", "w") as out:
        out.write("id\\tmean_plddt\\tptm\\tscrmsd\\tplddt_min\\tptm_min\\tscrmsd_max\\tpassed\\n")
        out.write(
            f"${meta.id}\\t{mean_plddt:.3f}\\t"
            f"{('%.4f' % ptm) if ptm is not None else 'NA'}\\t"
            f"{('%.3f' % scrmsd) if scrmsd is not None else 'NA'}\\t"
            f"${plddt_min}\\t${ptm_min}\\t${scrmsd_max}\\t{'yes' if passed else 'no'}\\n"
        )

    if not passed:
        try:
            os.remove("${fasta}")
        except OSError:
            pass

    with open("versions.yml", "w") as v:
        v.write('"${task.process}":\\n')
        v.write(f"    python: {sys.version.split()[0]}\\n")
    """
}
