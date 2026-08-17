/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CASCADE  —  cribado por métricas estructurales + ranking ponderado (YAML)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    candidatos (CSV id,sequence,...)
      -> FILTER_DEGENERATE_CANDIDATES   (limpia homopolimeros/entropia)
      -> CSV_TO_FASTAS                  (N fastas individuales)
      -> RUN_ESMFOLD                    (cribado barato)
      -> STRUCTURAL_METRICS             (dG, TM-score, pLDDT, Cys, carga, SAP)
      -> COLLECT_AND_RANK_STRUCTURAL    (gates + weights + top_n)   [ronda ESM]
      -> RUN_CHAI1                      (sobre el top-N de ESM)
      -> STRUCTURAL_METRICS             (mismas 6 metricas)
      -> COLLECT_AND_RANK_STRUCTURAL    (solo ordena, sin gates/top_n) [ronda Chai]
      => ranked.csv final (NO hay AF3)
*/

include { FILTER_DEGENERATE_CANDIDATES } from '../../modules/local/filter_degenerate_candidates/main'
include { CSV_TO_FASTAS                } from '../../modules/local/csv_to_fastas/main'
include { RUN_ESMFOLD                  } from '../../modules/local/run_esmfold'
include { RUN_CHAI1                    } from '../../modules/local/run_chai1'
include { STRUCTURAL_METRICS as METRICS_ESM  } from '../../modules/local/structural_metrics/main'
include { STRUCTURAL_METRICS as METRICS_CHAI } from '../../modules/local/structural_metrics/main'
include { REFERENCE_METRICS         } from '../../modules/local/structural_metrics/main'
include { COLLECT_AND_RANK_STRUCTURAL as RANK_ESM  } from '../../modules/local/collect_and_rank_structural/main'
include { COLLECT_AND_RANK_STRUCTURAL as RANK_CHAI } from '../../modules/local/collect_and_rank_structural/main'
include { CASCADE_REPORT             } from '../../modules/local/cascade_report/main'
include { PREPARE_ESMFOLD_DBS         } from '../local/prepare_esmfold_dbs'

workflow CASCADE {

    take:
    ch_candidates_csv
    ch_versions

    main:
    ch_metrics_reports = channel.empty()

    // 0. Limpiar degenerados
    FILTER_DEGENERATE_CANDIDATES(ch_candidates_csv)
    ch_versions = ch_versions.mix(FILTER_DEGENERATE_CANDIDATES.out.versions)

    // 1. CSV limpio -> N fastas
    CSV_TO_FASTAS(FILTER_DEGENERATE_CANDIDATES.out.clean)
    ch_versions = ch_versions.mix(CSV_TO_FASTAS.out.versions)

    ch_fastas = CSV_TO_FASTAS.out.fastas
        .flatten()
        .map { f -> tuple([id: f.baseName], f) }

    // 2. ESMFold (proceso directo)
    PREPARE_ESMFOLD_DBS(
        params.esmfold_db,
        params.esmfold_params_path,
        params.esmfold_3B_v1,
        params.esm2_t36_3B_UR50D,
        params.esm2_t36_3B_UR50D_contact_regression
    )
    RUN_ESMFOLD(
        ch_fastas,
        PREPARE_ESMFOLD_DBS.out.params,
        params.esmfold_num_recycles
    )
    ch_versions = ch_versions.mix(RUN_ESMFOLD.out.versions)

    // 3. Metricas estructurales sobre el PDB de ESMFold
    ch_reference_pdb  = channel.value(file(params.input_pdb, checkIfExists: true))
    ch_scoring_config = channel.value(file(params.scoring_config, checkIfExists: true))

    // 3a. ΔG y SAP del original UNA sola vez -> referencia de los deltas
    REFERENCE_METRICS(ch_reference_pdb, params.metrics_relax)
    ch_versions = ch_versions.mix(REFERENCE_METRICS.out.versions)

    // El TSV de referencia tiene: metric  value  (filas ref_dg y ref_sap).
    // Lo parseamos a un mapa y extraemos cada valor como escalar reutilizable.
    ch_ref_map = REFERENCE_METRICS.out.ref
        .splitCsv(header: true, sep: '\t')
        .map { row -> [ row.metric, row.value ] }
        .collect()                               // [[ref_dg, v], [ref_sap, v]] aplanado
        .map { items ->
            def m = [:]
            items.collate(2).each { pair -> m[pair[0]] = pair[1] }
            m
        }
    ch_ref_dg  = ch_ref_map.map { m -> m.ref_dg  ?: 'NA' }
    ch_ref_sap = ch_ref_map.map { m -> m.ref_sap ?: 'NA' }

    METRICS_ESM(
        RUN_ESMFOLD.out.top_ranked_pdb,
        ch_reference_pdb,
        params.metrics_ph,
        params.metrics_relax,
        ch_ref_dg,
        ch_ref_sap
    )
    ch_versions        = ch_versions.mix(METRICS_ESM.out.versions)
    ch_metrics_reports = ch_metrics_reports.mix(METRICS_ESM.out.metrics)

    // 4. Ranking ronda ESM: gates + top_n
    RANK_ESM(
        METRICS_ESM.out.metrics.map { meta, tsv -> tsv }.collect(),
        ch_scoring_config,
        true,
        true,
        'esm'
    )
    ch_versions = ch_versions.mix(RANK_ESM.out.versions)

    // 5. Top-N de ESM -> entrada de Chai
    ch_top_ids_esm = RANK_ESM.out.ranked
        .splitCsv(header: true)
        .filter { row -> row.in_top_n == 'yes' }
        .map { row -> row.id }

    ch_chai_in = ch_fastas
        .map { meta, fasta -> tuple(meta.id, meta, fasta) }
        .combine(ch_top_ids_esm.map { id -> tuple(id, true) }, by: 0)
        .map { id, meta, fasta, _flag -> tuple(meta, fasta) }

    // 6. Chai-1 (proceso directo) sobre supervivientes
    ch_chai_weights = channel.value(file(params.chai1_weights_path, checkIfExists: true))
    RUN_CHAI1(ch_chai_in, ch_chai_weights)
    ch_versions = ch_versions.mix(RUN_CHAI1.out.versions)

    // Chai emite el modelo top como CIF; el script de metricas lee CIF/PDB.
    METRICS_CHAI(
        RUN_CHAI1.out.top_ranked_cif,
        ch_reference_pdb,
        params.metrics_ph,
        params.metrics_relax,
        ch_ref_dg,
        ch_ref_sap
    )
    ch_versions        = ch_versions.mix(METRICS_CHAI.out.versions)
    ch_metrics_reports = ch_metrics_reports.mix(METRICS_CHAI.out.metrics)

    // 7. Ranking ronda Chai: SOLO ordenar
    RANK_CHAI(
        METRICS_CHAI.out.metrics.map { meta, tsv -> tsv }.collect(),
        ch_scoring_config,
        false,
        false,
        'chai'
    )
    ch_versions = ch_versions.mix(RANK_CHAI.out.versions)

    // 8. Report HTML de la cascada (dos pestañas, NGL, gráficos)
    //    Recolecta los PDB de ESM y los CIF de Chai (NGL lee ambos).
    ch_esm_pdbs  = RUN_ESMFOLD.out.top_ranked_pdb.map { meta, pdb -> pdb }.collect()
    ch_chai_pdbs = RUN_CHAI1.out.top_ranked_cif.map { meta, cif -> cif }.collect()

    CASCADE_REPORT(
        RANK_ESM.out.ranked,
        RANK_CHAI.out.ranked,
        ch_esm_pdbs,
        ch_chai_pdbs
    )

    emit:
    ranked_esm  = RANK_ESM.out.ranked
    ranked_chai = RANK_CHAI.out.ranked
    report      = CASCADE_REPORT.out.report
    metrics     = ch_metrics_reports
    versions    = ch_versions
}
