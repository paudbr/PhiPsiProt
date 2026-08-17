/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANTIBODY_DESIGN subworkflow v4 — pipeline completo secuencial
    ─────────────────────────────────────────────────────────────────────────────────
      Step 1  RFdiffusion    → N backbones (Quiver)
      Step 2  ProteinMPNN    → N × seqs_per_struct secuencias (Quiver)
      Step 3  AB_QV_TO_AF3_JSON → extrae PDBs + construye JSONs AF3
      Step 4  ESMFold2-Fast  → métricas por diseño
              FILTER_RANK_ESM → top 20% pasan a AF3
      Step 5  RUN_ALPHAFOLD3 → solo diseños que pasaron ESMFold2
              EXTRACT_AF3_AB_METRICS
              FILTER_RANK_AF3 → top 10 pasan a HADDOCK3
      Step 6  GENERATE_AB_HADDOCK_AIRS + RUN_HADDOCK3
              EXTRACT_HADDOCK3_AB_METRICS
              FILTER_RANK_HADDOCK → top 10 final
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RUN_RFANTIBODY              } from '../../modules/local/run_rfantibody/main'
include { RUN_PROTEINMPNN_RFAB        } from '../../modules/local/run_proteinmpnn_rfab/main'
include { AB_QV_TO_AF3_JSON           } from '../../modules/local/ab_qv_to_af3_json/main'
include { RUN_ESMFOLD2_AB             } from '../../modules/local/run_esmfold2_ab/main'
include { FILTER_RANK_AB as FILTER_RANK_ESM     } from '../../modules/local/filter_rank_ab/main'
include { FILTER_RANK_AB as FILTER_RANK_AF3     } from '../../modules/local/filter_rank_ab/main'
include { FILTER_RANK_AB as FILTER_RANK_HADDOCK } from '../../modules/local/filter_rank_ab/main'
include { PREPARE_ALPHAFOLD3_DBS      } from '../local/prepare_alphafold3_dbs'
include { RUN_ALPHAFOLD3              } from '../../modules/local/run_alphafold3/main'
include { EXTRACT_AF3_AB_METRICS      } from '../../modules/local/extract_af3_ab_metrics/main'
include { GENERATE_AB_HADDOCK_AIRS    } from '../../modules/local/generate_ab_haddock_airs/main'
include { RUN_HADDOCK3                } from '../../modules/local/run_haddock3/main'
include { EXTRACT_HADDOCK3_AB_METRICS } from '../../modules/local/extract_haddock3_ab_metrics/main'
include { GENERATE_AB_REPORT          } from '../../modules/local/generate_ab_report/main'
include { EXTRACT_PDB_FROM_QV         } from '../../modules/local/extract_pdb_from_qv/main'

workflow ANTIBODY_DESIGN {

    take:
    mode

    main:

    ch_versions = channel.empty()

    // ─── Validate required params ──────────────────────────────────────────
    if (!params.ab_target_pdb)       error "[ANTIBODY_DESIGN] --ab_target_pdb required"
    if (!params.ab_framework_pdb)    error "[ANTIBODY_DESIGN] --ab_framework_pdb required"
    if (!params.ab_hotspot_residues) error "[ANTIBODY_DESIGN] --ab_hotspot_residues required"
    if (!params.ab_design_loops)     error "[ANTIBODY_DESIGN] --ab_design_loops required"
    if (!params.ab_rfdiff_weights)   error "[ANTIBODY_DESIGN] --ab_rfdiff_weights required"
    if (!params.ab_mpnn_weights)     error "[ANTIBODY_DESIGN] --ab_mpnn_weights required"
    if (!params.ab_esmfold2_weights) error "[ANTIBODY_DESIGN] --ab_esmfold2_weights required"
    if (!params.alphafold3_db)       error "[ANTIBODY_DESIGN] --alphafold3_db required"

    def ab_type     = params.ab_type        ?: 'VHH'
    def n_seeds     = params.ab_af3_seeds   ?: 10
    def run_haddock = params.ab_run_haddock != null ? params.ab_run_haddock : true

    // ─── Scoring YAML configs ─────────────────────────────────────────────
    ch_esm_yaml     = channel.value(file(params.ab_scoring_esm,      checkIfExists: true))
    ch_af3_yaml     = channel.value(file(params.ab_scoring_af3,      checkIfExists: true))
    ch_haddock_yaml = channel.value(file(params.ab_scoring_haddock3, checkIfExists: true))

    // ─── Static input channels ────────────────────────────────────────────
    ch_target_pdb    = channel.value(file(params.ab_target_pdb,       checkIfExists: true))
    ch_framework_pdb = channel.value(file(params.ab_framework_pdb,    checkIfExists: true))
    ch_rfdiff_wts    = channel.value(file(params.ab_rfdiff_weights,   checkIfExists: true))
    ch_mpnn_wts      = channel.value(file(params.ab_mpnn_weights,     checkIfExists: true))
    ch_esm2_wts      = channel.value(file(params.ab_esmfold2_weights, checkIfExists: true))

    ch_input = channel.of([
        [ id: "ab_design_${ab_type}" ],
        file(params.ab_target_pdb,    checkIfExists: true),
        file(params.ab_framework_pdb, checkIfExists: true)
    ])

    // ─── STEP 1: RFdiffusion ──────────────────────────────────────────────
    RUN_RFANTIBODY(
        ch_input,
        ch_rfdiff_wts,
        params.ab_hotspot_residues,
        params.ab_design_loops,
        params.ab_num_designs ?: 50
    )
    ch_versions = ch_versions.mix(RUN_RFANTIBODY.out.versions)

    // ─── STEP 2: ProteinMPNN ─────────────────────────────────────────────
    RUN_PROTEINMPNN_RFAB(
        RUN_RFANTIBODY.out.quiver,
        ch_mpnn_wts,
        params.ab_seqs_per_struct  ?: 4,
        params.ab_mpnn_temperature ?: 0.2
    )
    ch_versions = ch_versions.mix(RUN_PROTEINMPNN_RFAB.out.versions)

    // ─── STEP 3: Quiver → PDBs + AF3 JSONs ───────────────────────────────
    AB_QV_TO_AF3_JSON(
        RUN_PROTEINMPNN_RFAB.out.quiver,
        ch_target_pdb,
        ab_type,
        n_seeds
    )
    ch_versions = ch_versions.mix(AB_QV_TO_AF3_JSON.out.versions)

    // Un canal por diseño: [ meta(id=design_id), json_file ]
    ch_jsons = AB_QV_TO_AF3_JSON.out.jsons
        .transpose()
        .map { meta, json_file ->
            def m = meta.clone()
            m.id  = json_file.baseName
            [ m, json_file ]
        }

    // ─── STEP 4: ESMFold2 pre-filtro — BATCH MODE ────────────────────────
    // EXTRACT_PDB_FROM_QV saca todos los PDBs del Quiver.
    // RUN_ESMFOLD2_AB recibe el directorio completo, carga el modelo
    // UNA SOLA VEZ y procesa todos en un loop → evita OOM de ESMC-6B.

    ch_qv_and_ids = ch_jsons
        .map { meta, json -> json.baseName }
        .combine( RUN_PROTEINMPNN_RFAB.out.quiver.map { meta, qv -> qv } )
        .map { design_id, qv ->
            [ [ id: design_id ], qv, design_id ]
        }

    EXTRACT_PDB_FROM_QV( ch_qv_and_ids )

    ch_designs_dir = EXTRACT_PDB_FROM_QV.out.pdb
            .map { meta, pdb -> pdb }
            .collect()
            .map { pdbs ->
                // Crear directorio staging y copiar todos los PDBs
                def staging = file("${workDir}/esm2_staging")
                staging.mkdirs()
                pdbs.each { pdb ->
                    def dest = staging.resolve(pdb.name)
                    if (!dest.exists()) pdb.copyTo(dest)
                }
                [ [ id: "ab_design_${ab_type}" ], staging ]
            }

    ch_for_esm = ch_designs_dir
        .combine(ch_target_pdb)
        .combine(ch_framework_pdb)
        .map { meta, ddir, tgt, fw -> [ meta, ddir, tgt, fw ] }

    RUN_ESMFOLD2_AB( ch_for_esm )
    ch_versions = ch_versions.mix(RUN_ESMFOLD2_AB.out.versions)

    // El batch ya produce un solo TSV con todas las métricas
    ch_esm_metrics_all = RUN_ESMFOLD2_AB.out.metrics
        .map { meta, tsv -> [ meta, "esmfold2", tsv ] }

    // Filtrar y rankear — top 20%
    FILTER_RANK_ESM( ch_esm_metrics_all, ch_esm_yaml )
    ch_versions = ch_versions.mix(FILTER_RANK_ESM.out.versions)

    // ─── STEP 5: AF3 — SOLO diseños que pasaron ESMFold2 ─────────────────
    ch_esm_passed_set = FILTER_RANK_ESM.out.passed_ids
        .map { meta, stage, txt ->
            txt.readLines().collect { it.trim() }.findAll { it } as Set
        }

    // Filtrar JSONs: solo los cuyo meta.id está en el set de pasados
    ch_af3_input = ch_jsons
        .combine(ch_esm_passed_set)
        .filter { meta, json, passed_set -> meta.id in passed_set }
        .map    { meta, json, passed_set -> [ meta, json ] }

    PREPARE_ALPHAFOLD3_DBS(
        params.alphafold3_db,
        params.alphafold3_params_path,
        params.alphafold3_small_bfd_path,
        params.alphafold3_mgnify_path,
        params.alphafold3_pdb_mmcif_path,
        params.alphafold3_uniref90_path,
        params.alphafold3_pdb_seqres_path,
        params.alphafold3_uniprot_path,
        params.alphafold3_rnacentral_path,
        params.alphafold3_nt_rna_path,
        params.alphafold3_rfam_path,
        params.alphafold3_small_bfd_link,
        params.alphafold3_mgnify_link,
        params.alphafold3_pdb_mmcif_link,
        params.alphafold3_uniref90_link,
        params.alphafold3_pdb_seqres_link,
        params.alphafold3_uniprot_link,
        params.alphafold3_rnacentral_link,
        params.alphafold3_nt_rna_link,
        params.alphafold3_rfam_link
    )
    ch_versions = ch_versions.mix(PREPARE_ALPHAFOLD3_DBS.out.versions)

    RUN_ALPHAFOLD3(
        ch_af3_input,
        PREPARE_ALPHAFOLD3_DBS.out.params,
        PREPARE_ALPHAFOLD3_DBS.out.small_bfd,
        PREPARE_ALPHAFOLD3_DBS.out.mgnify,
        PREPARE_ALPHAFOLD3_DBS.out.pdb_mmcif,
        PREPARE_ALPHAFOLD3_DBS.out.uniref90,
        PREPARE_ALPHAFOLD3_DBS.out.pdb_seqres,
        PREPARE_ALPHAFOLD3_DBS.out.uniprot
    )
    ch_versions = ch_versions.mix(RUN_ALPHAFOLD3.out.versions)

    EXTRACT_AF3_AB_METRICS(
        RUN_ALPHAFOLD3.out.iptms
            .join(RUN_ALPHAFOLD3.out.ptms)
            .join(RUN_ALPHAFOLD3.out.top_ranked_cif)
    )
    ch_versions = ch_versions.mix(EXTRACT_AF3_AB_METRICS.out.versions)

    ch_af3_metrics_all = EXTRACT_AF3_AB_METRICS.out.metrics
        .map    { meta, tsv -> tsv }
        .collectFile(
            name:       "af3_all_metrics.tsv",
            keepHeader: true,
            skip:       1
        )
        .map { tsv -> [ [ id: "ab_design_${ab_type}" ], "af3", tsv ] }

    FILTER_RANK_AF3( ch_af3_metrics_all, ch_af3_yaml )
    ch_versions = ch_versions.mix(FILTER_RANK_AF3.out.versions)

    // ─── STEP 6: HADDOCK3 — SOLO diseños que pasaron AF3 ─────────────────
    ch_haddock_scores = channel.empty()
    ch_haddock_ranked = channel.empty()

    if (run_haddock) {

        ch_af3_passed_set = FILTER_RANK_AF3.out.passed_ids
            .map { meta, stage, txt ->
                txt.readLines().collect { it.trim() }.findAll { it } as Set
            }

        ch_passed_cifs = EXTRACT_AF3_AB_METRICS.out.cif
            .combine(ch_af3_passed_set)
            .filter { meta, cif, passed_set -> meta.id in passed_set }
            .map    { meta, cif, passed_set -> [ meta, cif ] }

        ch_for_airs = ch_passed_cifs
            .combine(ch_target_pdb)
            .map { meta, cif, tgt -> [ meta, cif, tgt ] }

        GENERATE_AB_HADDOCK_AIRS(
            ch_for_airs,
            params.ab_hotspot_residues
        )

        RUN_HADDOCK3( GENERATE_AB_HADDOCK_AIRS.out.docking_input )
        ch_versions = ch_versions.mix(RUN_HADDOCK3.out.versions)

        EXTRACT_HADDOCK3_AB_METRICS( RUN_HADDOCK3.out.scores )
        ch_versions = ch_versions.mix(EXTRACT_HADDOCK3_AB_METRICS.out.versions)

        ch_haddock_metrics_all = EXTRACT_HADDOCK3_AB_METRICS.out.metrics
            .map    { meta, tsv -> tsv }
            .collectFile(
                name:       "haddock_all_metrics.tsv",
                keepHeader: true,
                skip:       1
            )
            .map { tsv -> [ [ id: "ab_design_${ab_type}" ], "haddock3", tsv ] }

        FILTER_RANK_HADDOCK( ch_haddock_metrics_all, ch_haddock_yaml )
        ch_versions       = ch_versions.mix(FILTER_RANK_HADDOCK.out.versions)
        ch_haddock_scores = RUN_HADDOCK3.out.scores
        ch_haddock_ranked = FILTER_RANK_HADDOCK.out.ranked
    }

    // ─── Summary + Report ─────────────────────────────────────────────────
    ch_iptm_summary = FILTER_RANK_AF3.out.ranked
        .map { meta, stage, tsv -> tsv }
        .collectFile(
            name:     "antibody_design_af3_ranking.tsv",
            storeDir: "${params.outdir}/antibody_design/"
        )

    // Collect ranked TSVs for the report
    ch_ranked_esm_tsv = FILTER_RANK_ESM.out.ranked
        .map { meta, stage, tsv -> tsv }
        .collectFile(name: "ranked_esmfold2.tsv", keepHeader: true, skip: 1)

    ch_ranked_af3_tsv = FILTER_RANK_AF3.out.ranked
        .map { meta, stage, tsv -> tsv }
        .collectFile(name: "ranked_af3.tsv", keepHeader: true, skip: 1)

    ch_ranked_haddock_tsv = run_haddock ?
        FILTER_RANK_HADDOCK.out.ranked
            .map { meta, stage, tsv -> tsv }
            .collectFile(name: "ranked_haddock3.tsv", keepHeader: true, skip: 1) :
        channel.value(file("$projectDir/assets/NO_FILE"))

    // Collect CIF dirs for NGL viewer
    ch_pdbs_esm = RUN_ESMFOLD2_AB.out.cifs
        .map { meta, cif_dir -> cif_dir }

    ch_pdbs_af3 = EXTRACT_AF3_AB_METRICS.out.cif
        .map { meta, cif -> cif }
        .collect()
        .map { files -> files }

    // Generate HTML report
    GENERATE_AB_REPORT(
        ch_ranked_esm_tsv,
        ch_ranked_af3_tsv,
        ch_ranked_haddock_tsv,
        ch_pdbs_esm.ifEmpty(file("$projectDir/assets/NO_FILE")),
        ch_pdbs_af3.ifEmpty(file("$projectDir/assets/NO_FILE"))
    )
    ch_versions = ch_versions.mix(GENERATE_AB_REPORT.out.versions)

    emit:
    ranked_esm     = FILTER_RANK_ESM.out.ranked
    ranked_af3     = FILTER_RANK_AF3.out.ranked
    ranked_haddock = ch_haddock_ranked
    iptm_summary   = ch_iptm_summary
    haddock_scores = ch_haddock_scores
    report         = GENERATE_AB_REPORT.out.report
    versions       = ch_versions
}

