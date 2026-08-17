/*
 * Explota un CSV (id,sequence,...) en N ficheros FASTA individuales,
 * uno por candidato. Cada FASTA es la entrada por-muestra que esperan
 * ESMFOLD / CHAI1 / ALPHAFOLD3 (tuple(meta, path(fasta))).
 */
process CSV_TO_FASTAS {
    tag "csv_to_fastas"
    label 'process_single'
    container "nf-core/proteinfold_esmfold:2.0.0"

    input:
    path candidates_csv

    output:
    path "fastas/*.fasta", emit: fastas
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    import csv, os, re

    os.makedirs("fastas", exist_ok=True)

    def safe(s):
        return re.sub(r"[^A-Za-z0-9_.-]", "_", str(s).strip())

    with open("${candidates_csv}") as fh:
        reader = csv.DictReader(fh)
        if "design_id" not in reader.fieldnames or "final_sequence" not in reader.fieldnames:
            raise SystemExit(
                f"CSV debe tener columnas 'design_id' y 'final_sequence'. Encontradas: {reader.fieldnames}"
            )
        n = 0
        for row in reader:
            cid = safe(row["design_id"])
            seq = str(row["final_sequence"]).strip().upper()
            if not cid or not seq:
                continue
            with open(os.path.join("fastas", f"{cid}.fasta"), "w") as out:
                out.write(f">{cid}\\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i+60] + "\\n")
            n += 1

    if n == 0:
        raise SystemExit("No se escribió ningún FASTA: CSV vacío o columnas mal.")

    with open("versions.yml", "w") as v:
        import sys
        v.write('"${task.process}":\\n')
        v.write(f"    python: {sys.version.split()[0]}\\n")
    print(f"Escritos {n} fastas.")
    """
}
