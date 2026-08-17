/*
 * scRMSD con Biopython (sin binarios externos ni paquetes nuevos).
 *
 * Superpone el modelo de ESMFold sobre un PDB de referencia único
 * (params.input_pdb) emparejando los CA por (cadena, número de residuo)
 * comunes a ambas estructuras, aplica Kabsch (Bio.SVDSuperimposer) y
 * reporta el RMSD global sobre esos CA.
 *
 * El corte por umbral lo aplica FILTER_BY_METRICS aguas abajo (lee este TSV).
 *
 * Requiere SOLO biopython, que ya está en la imagen de ESMFold.
 *
 * Salida TSV: id, scrmsd, n_ca_aligned, n_ca_model, n_ca_ref
 *   - scrmsd = NA si algo falla o si no hay CA emparejables.
 *   - n_ca_aligned permite detectar emparejamientos pobres (RMSD engañoso).
 */
process SCRMSD_USALIGN {
    tag "$meta.id"
    label 'process_single'
    container "nf-core/proteinfold_esmfold:2.0.0"

    input:
    tuple val(meta), path(model_pdb)
    path  reference_pdb

    output:
    tuple val(meta), path("${meta.id}_scrmsd.tsv"), emit: scrmsd
    path  "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    import sys
    import warnings
    warnings.filterwarnings("ignore")

    scrmsd      = "NA"
    n_aligned   = 0
    n_ca_model  = 0
    n_ca_ref    = 0

    try:
        import numpy as np
        from Bio.PDB import PDBParser
        from Bio.SVDSuperimposer import SVDSuperimposer

        parser = PDBParser(QUIET=True)

        def ca_map(path):
            # Devuelve dict {(chain_id, resseq, icode): coord_CA}
            # Solo primer modelo. Ignora HETATM (residuos no estándar via id[0]).
            s = parser.get_structure("s", path)
            model = next(s.get_models())
            out = {}
            for chain in model:
                for res in chain:
                    if res.id[0] != " ":      # salta HETATM/agua
                        continue
                    if "CA" not in res:
                        continue
                    key = (chain.id, res.id[1], res.id[2])
                    out[key] = res["CA"].get_coord()
            return out

        m = ca_map("${model_pdb}")
        r = ca_map("${reference_pdb}")
        n_ca_model = len(m)
        n_ca_ref   = len(r)

        # Emparejar por clave (cadena, resseq, icode) comun
        common = sorted(set(m.keys()) & set(r.keys()))
        n_aligned = len(common)

        if n_aligned >= 3:
            coords_m = np.array([m[k] for k in common], dtype=float)
            coords_r = np.array([r[k] for k in common], dtype=float)
            sup = SVDSuperimposer()
            # set(reference, moving) -> superpone 'moving' sobre 'reference'
            sup.set(coords_r, coords_m)
            sup.run()
            scrmsd = "%.3f" % sup.get_rms()
        else:
            print(f"[WARN] solo {n_aligned} CA comunes; RMSD no fiable", file=sys.stderr)

    except Exception as e:
        print(f"[WARN] scRMSD falló para ${meta.id}: {e}", file=sys.stderr)

    with open("${meta.id}_scrmsd.tsv", "w") as out:
        out.write("id\\tscrmsd\\tn_ca_aligned\\tn_ca_model\\tn_ca_ref\\n")
        out.write(f"${meta.id}\\t{scrmsd}\\t{n_aligned}\\t{n_ca_model}\\t{n_ca_ref}\\n")

    with open("versions.yml", "w") as v:
        v.write('"${task.process}":\\n')
        v.write(f"    python: {sys.version.split()[0]}\\n")
        try:
            import Bio
            v.write(f"    biopython: {Bio.__version__}\\n")
        except Exception:
            pass
    """
}
