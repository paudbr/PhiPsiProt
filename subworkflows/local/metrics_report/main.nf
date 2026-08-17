include { COLLECT_METRICS         } from '../../../modules/local/collect_metrics/main'
include { GENERATE_METRICS_REPORT } from '../../../modules/local/generate_metrics_report/main'


process EXTRACT_CHAI_ALL_MODELS {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(raw_files, stageAs: "raw/*")

    output:
    tuple val(meta), path("${meta.id}_chai1_allmodels_plddt.tsv"), emit: plddt
    tuple val(meta), path("${meta.id}_chai1_scores.tsv"),          emit: scores

    script:
    """
    SAMPLE_ID="${meta.id}" python3 ${moduleDir}/process_chai.py
    """
}

process MERGE_CHAI_SCORES {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(ptm_file), path(iptm_file)

    output:
    path "${meta.id}_chai_scores.tsv"

    script:
    // Usamos printf para escribir el TSV — evita el problema de \t\n
    // siendo interpolados por Groovy dentro del bloque """ """
    """
    python3 - << 'EOF'
import os

prefix = "${meta.id}_chai"
ptm_f  = "${ptm_file}"
iptm_f = "${iptm_file}"

TAB = chr(9)
NL  = chr(10)

def read_val(path):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return None
    with open(path) as f:
        content = f.read().strip()
    lines = content.splitlines()
    if len(lines) == 1:
        try:
            return float(lines[0])
        except ValueError:
            return None
    for line in lines[1:]:
        parts = line.split()
        if parts:
            try:
                return float(parts[-1])
            except ValueError:
                pass
    return None

ptm  = read_val(ptm_f)
iptm = read_val(iptm_f)

with open(f"{prefix}_scores.tsv", "w") as fh:
    fh.write("metric" + TAB + "value" + NL)
    if ptm  is not None: fh.write("ptm"  + TAB + f"{ptm:.4f}"  + NL)
    if iptm is not None: fh.write("iptm" + TAB + f"{iptm:.4f}" + NL)
EOF
    """
}

process MERGE_AF3_SCORES {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(ptm_file), path(iptm_file)

    output:
    path "${meta.id}_af3_scores.tsv"

    script:
    """
    python3 - << 'EOF'
import os

prefix = "${meta.id}_af3"
ptm_f  = "${ptm_file}"
iptm_f = "${iptm_file}"

TAB = chr(9)
NL  = chr(10)

def read_val(path):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return None
    with open(path) as f:
        content = f.read().strip()
    lines = content.splitlines()
    if len(lines) == 1:
        try:
            return float(lines[0])
        except ValueError:
            return None
    for line in lines[1:]:
        parts = line.split()
        if parts:
            try:
                return float(parts[-1])
            except ValueError:
                pass
    return None

ptm  = read_val(ptm_f)
iptm = read_val(iptm_f)

with open(f"{prefix}_scores.tsv", "w") as fh:
    fh.write("metric" + TAB + "value" + NL)
    if ptm  is not None: fh.write("ptm"  + TAB + f"{ptm:.4f}"  + NL)
    if iptm is not None: fh.write("iptm" + TAB + f"{iptm:.4f}" + NL)
EOF
    """
}

process TAG_AF3_PLDDT {
    tag "$meta.id"
    input:  tuple val(meta), path(tsv)
    output: path("${meta.id}_af3_plddt.tsv")
    script: "cp ${tsv} ${meta.id}_af3_plddt.tsv"
}

process TAG_AF3_PDB {
    tag "$meta.id"
    input:  tuple val(meta), path(pdb)
    output: path("${meta.id}_af3_${pdb.baseName}.pdb")
    script: "cp ${pdb} ${meta.id}_af3_${pdb.baseName}.pdb"
}

workflow METRICS_REPORT {

    take:
    ch_af2_pdb    // tuple val(meta), path(pdb)     — o Channel.empty()
    ch_af3_pdb    // tuple val(meta), path(pdb/cif) — o Channel.empty()
    ch_esm_pdb    // tuple val(meta), path(pdb)     — o Channel.empty()
    ch_chai_plddt // tuple val(meta), path(tsv)     — o Channel.empty()
    ch_chai_ptm   // tuple val(meta), path(tsv/txt) — o Channel.empty()
    ch_chai_iptm  // tuple val(meta), path(tsv/txt) — o Channel.empty()
    ch_chai_pdb 
    ch_chai_raw
    ch_af3_plddt
    ch_af3_ptm
    ch_af3_iptm

    main:

    def NO_FILE = file("$projectDir/assets/NO_FILE")

    // ── AF2: adjunta el json usando crossValue por id; si no hay json usa NO_FILE
    // Estrategia: tag con tool antes de cualquier join, así evitamos operar
    // sobre Channel.empty() con map/join — usamos mix+filter en su lugar.

    ch_af2_tagged = ch_af2_pdb.map { meta, pdb -> [ meta.id, meta, pdb, "af2" ] }
    ch_af3_tagged = ch_af3_pdb.map { meta, pdb -> [ meta.id, meta, pdb, "af3" ] }
    ch_esm_tagged = ch_esm_pdb.map { meta, pdb -> [ meta.id, meta, pdb, "esm" ] }
    ch_chai_tagged = ch_chai_pdb.map { meta, pdb -> [ meta.id, meta, pdb, "chai1" ] }
    // Join con remainder:true SOLO sobre canales del mismo tipo (pdb+json del mismo tool)
    // Nunca operamos sobre Channel.empty() directamente con map
    ch_af2_ready = ch_af2_tagged
        .map { meta_id, meta, pdb, tool -> [ meta, pdb, tool, NO_FILE ] }

    ch_af3_ready = ch_af3_tagged
        .flatMap { meta_id, meta, pdb, tool ->
            def fileList = (pdb instanceof List) ? pdb.flatten() : [ pdb ]
            fileList.collect { f ->
                def m = meta.clone()
                [ m, f, tool, NO_FILE ]
            }
        }

    ch_esm_ready = ch_esm_tagged
        .map { meta_id, meta, pdb, tool -> [ meta, pdb, tool, NO_FILE ] }

    ch_chai_ready = ch_chai_tagged
        .flatMap { meta_id, meta, pdb, tool ->
            def fileList = (pdb instanceof List) ? pdb.flatten() : [ pdb ]
            fileList.collect { f -> f }   // path suelto — ch_all_pdbs lo mezcla directamente
        }
    // Mezcla todo y llama COLLECT_METRICS una sola vez
    COLLECT_METRICS(
        ch_af2_ready.mix(ch_esm_ready )
    )


    ch_all_plddt  = COLLECT_METRICS.out.plddt .map { meta, f -> f }
    ch_all_scores = COLLECT_METRICS.out.scores.map { meta, f -> f }

    TAG_AF3_PLDDT( ch_af3_plddt )

    ch_all_plddt = ch_all_plddt.mix( TAG_AF3_PLDDT.out )

    MERGE_AF3_SCORES(
        ch_af3_ptm
            .map  { meta, f -> [ meta.id, meta, f ] }
            .join ( ch_af3_iptm.map { meta, f -> [ meta.id, f ] } )
            .map  { id, meta, ptm, iptm -> [ meta, ptm, iptm ] }
    )
    ch_all_scores = ch_all_scores.mix( MERGE_AF3_SCORES.out )

    EXTRACT_CHAI_ALL_MODELS( ch_chai_raw )

    ch_all_plddt  = ch_all_plddt .mix( EXTRACT_CHAI_ALL_MODELS.out.plddt .map { meta, f -> f } )
    ch_all_scores = ch_all_scores.mix( EXTRACT_CHAI_ALL_MODELS.out.scores.map { meta, f -> f } )

    // ── Chai scores: join ptm + iptm por id ───────────────────────────────
    MERGE_CHAI_SCORES(
        ch_chai_ptm
            .map  { meta, f -> [ meta.id, meta, f ] }
            .join ( ch_chai_iptm.map { meta, f -> [ meta.id, f ] } )
            .map  { id, meta, ptm, iptm -> [ meta, ptm, iptm ] }
    )
    ch_all_scores = ch_all_scores.mix( MERGE_CHAI_SCORES.out )

        // PDBs de AF3 renombrados con meta.id para evitar colisión en stageAs
    TAG_AF3_PDB(
        ch_af3_ready.map { meta, pdb, tool, json -> [ meta, pdb ] }
    )

    ch_all_pdbs = ch_af2_ready
        .mix( ch_esm_ready )
        .map { meta, pdb, tool, json -> pdb }   // af2 + esm
        .mix( ch_chai_ready )                    // chai (path suelto)
        .mix( TAG_AF3_PDB.out )                  // af3 renombrado (path suelto)

    // ── Reporte único ─────────────────────────────────────────────────────
    GENERATE_METRICS_REPORT(
        ch_all_plddt .collect(),
        ch_all_scores.collect(),
        ch_all_pdbs  .flatten().collect()
    )

    emit:
    report = GENERATE_METRICS_REPORT.out.report
}
