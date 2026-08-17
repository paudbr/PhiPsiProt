include { RUN_GNINA_LIGAND_FILTER                  } from '../../modules/local/run_gnina_ligand_filter/main'
include { CANDIDATES_TO_BOLTZ                      } from '../../modules/local/candidates_to_boltz/main'
include { CANDIDATES_SMILES_TO_BOLTZ               } from '../../modules/local/candidates_smiles_to_boltz/main'
include { RUN_BOLTZ as RUN_BOLTZ_LIGAND            } from '../../modules/local/run_boltz/main'
include { RUN_BOLTZ as RUN_BOLTZ_SM               } from '../../modules/local/run_boltz/main'
include { BOLTZ_LIGAND_SUMMARY                     } from '../../modules/local/boltz_ligand_summary/main'
include { BOLTZ_LIGAND_SUMMARY as BOLTZ_SM_SUMMARY } from '../../modules/local/boltz_ligand_summary/main'
include { CANDIDATES_TO_COFOLD_SAMPLESHEET         } from '../../modules/local/candidates_to_cofold_samplesheet/main'
include { CANDIDATES_CSV_TO_FASTA                  } from '../../modules/local/candidates_csv_to_fasta/main'
include { ESMFOLD                                  } from '../../workflows/esmfold'
include { RUN_HADDOCK3                             } from '../../modules/local/run_haddock3/main'
include { RUN_ROSETTADOCK                          } from '../../modules/local/run_rosettadock/main'
include { GENERATE_HADDOCK_AIRS                    } from '../../modules/local/generate_haddock_airs/main'
include { FILTER_DEGENERATE_CANDIDATES             } from '../../modules/local/filter_degenerate_candidates/main'

workflow DOCKING_COFOLDING {

    take:
    candidates_csv
    reference_pdb
    receptor_sequence
    ligand_file
    ch_boltz_model
    ch_boltz_ccd
    ch_boltz2_aff
    ch_boltz2_conf
    ch_boltz2_mols
    esmfold_params
    esmfold_num_recycles

    main:
    ch_versions = channel.empty()
    def is_prot_prot = params.docking_tool in ['haddock3', 'rosettadock']

    // ──────────────────────────────────────────────────────────────
    // RAMA PROTEIN-PROTEIN: ESMFold + Boltz2 + HADDOCK3/RosettaDock
    // ──────────────────────────────────────────────────────────────
    if (is_prot_prot) {

        // QC: filtrar candidatos degenerados antes de cualquier folding
        FILTER_DEGENERATE_CANDIDATES(candidates_csv)
        ch_candidates_clean = FILTER_DEGENERATE_CANDIDATES.out.clean

        // Generar YAMLs para Boltz2
        CANDIDATES_TO_BOLTZ(ch_candidates_clean, receptor_sequence)

        ch_boltz_input = CANDIDATES_TO_BOLTZ.out.yaml
            .flatMap { files ->
                def fileList = (files instanceof List) ? files.flatten() : [ files ]
                fileList.collect { yaml ->
                    def runId       = yaml.baseName
                    def candidateId = runId.replaceFirst(/_rep\d+$/, '')
                    [ [ id: runId, candidate_id: candidateId, model: "boltz2_ligand" ], yaml, [] ]
                }
            }

        // Generar FASTAs para ESMFold
        CANDIDATES_CSV_TO_FASTA(ch_candidates_clean)

        ch_esmfold_samplesheet = CANDIDATES_CSV_TO_FASTA.out.fastas
            .flatten()
            .map { fasta -> tuple([id: fasta.baseName], fasta) }

        // ESMFold
        ESMFOLD(
            ch_esmfold_samplesheet,
            channel.empty(),
            esmfold_params,
            esmfold_num_recycles
        )

        // Boltz2 co-folding + affinity
        RUN_BOLTZ_LIGAND(
            ch_boltz_input,
            ch_boltz_model,
            ch_boltz_ccd,
            ch_boltz2_aff,
            ch_boltz2_conf,
            ch_boltz2_mols
        )

        BOLTZ_LIGAND_SUMMARY(
            RUN_BOLTZ_LIGAND.out.confidence
                .map { meta, json -> json }
                .collect()
        )

        CANDIDATES_TO_COFOLD_SAMPLESHEET(
            BOLTZ_LIGAND_SUMMARY.out.summary,
            ch_candidates_clean,
            receptor_sequence
        )

        // AIRs para HADDOCK3
        ch_air_input = ESMFOLD.out.pdb
            .map { meta, pdb ->
                tuple(meta.id, file(params.ligand_reference_pdb), pdb)
            }

        GENERATE_HADDOCK_AIRS(ch_air_input)

        // Join estructuras + restraints para HADDOCK3
        ch_haddock_ready = ESMFOLD.out.pdb
            .map { meta, pdb -> tuple(meta.id, pdb) }
            .join(
                GENERATE_HADDOCK_AIRS.out.airs
                    .map { candidate_id, restraints, active, passive ->
                        tuple(candidate_id, restraints)
                    }
            )
            .map { candidate_id, binder_pdb, restraints ->
                tuple(
                    candidate_id,
                    file(params.ligand_reference_pdb),
                    binder_pdb,
                    restraints
                )
            }

        RUN_HADDOCK3(ch_haddock_ready)

        // Acumular versiones — un solo mix al final
        ch_versions = channel.empty().mix(
            FILTER_DEGENERATE_CANDIDATES.out.versions,
            CANDIDATES_TO_BOLTZ.out.versions,
            CANDIDATES_CSV_TO_FASTA.out.versions,
            ESMFOLD.out.versions,
            RUN_BOLTZ_LIGAND.out.versions,
            BOLTZ_LIGAND_SUMMARY.out.versions,
            CANDIDATES_TO_COFOLD_SAMPLESHEET.out.versions,
            GENERATE_HADDOCK_AIRS.out.versions,
            RUN_HADDOCK3.out.versions
        )

    // ──────────────────────────────────────────────────────────────
    // RAMA SMALL MOLECULE: GNINA + Boltz2
    // ──────────────────────────────────────────────────────────────
    }else if (params.mode != 'structural') {
    // modos de diseño sin ligando: los candidatos son la entrada estructural
    ch_samplesheet = ch_design_candidates
        .splitFasta(record: [id: true, sequence: true])
        .map { record -> tuple([id: record.id], record.sequence)}
    
    } else {
        // Docking clasico con GNINA
        RUN_GNINA_LIGAND_FILTER(
            channel.fromPath(params.ligand_candidates_csv),
            reference_pdb,
            ligand_file
        )

        
        // Leer el CSV (candidate_id,smiles_ligand) y construir un canal por ligando
        ch_ligands = candidates_csv
            .splitCsv(header: true)
            .map { row ->
                [ [ id: row.candidate_id, candidate_id: row.candidate_id, model: "boltz2_ligand" ],
                  row.smiles_ligand ]
            }

        // Generar YAMLs de Boltz2 (receptor fijo + SMILES + affinity)
        CANDIDATES_SMILES_TO_BOLTZ(ch_ligands, receptor_sequence, params.pocket_residues)

        ch_boltz_input = CANDIDATES_SMILES_TO_BOLTZ.out.yaml
            .map { meta, yaml -> [ meta, yaml, [] ] }

        // Boltz2 co-folding + affinity (alias propio para evitar invocacion duplicada)
        RUN_BOLTZ_SM(
            ch_boltz_input,
            ch_boltz_model,
            ch_boltz_ccd,
            ch_boltz2_aff,
            ch_boltz2_conf,
            ch_boltz2_mols
        )

        BOLTZ_SM_SUMMARY(
            RUN_BOLTZ_SM.out.confidence
                .map { meta, json -> json }
                .collect()
        )

        ch_versions = RUN_GNINA_LIGAND_FILTER.out.versions

    }

    // ──────────────────────────────────────────────────────────────
    // Outputs — channel.empty() para la rama inactiva
    // ──────────────────────────────────────────────────────────────
    emit:
            // --- métricas crudas (small molecule: GNINA + Boltz van JUNTOS) ---
            gnina_scores    = is_prot_prot ? channel.empty() : RUN_GNINA_LIGAND_FILTER.out.scores
            gnina_summary   = is_prot_prot ? channel.empty() : RUN_GNINA_LIGAND_FILTER.out.summary
            boltz_scores    = is_prot_prot ? BOLTZ_LIGAND_SUMMARY.out.scores  : BOLTZ_SM_SUMMARY.out.scores
            boltz_summary   = is_prot_prot ? BOLTZ_LIGAND_SUMMARY.out.summary : BOLTZ_SM_SUMMARY.out.summary
            boltz_dirs      = is_prot_prot ? RUN_BOLTZ_LIGAND.out.intermediates.map { meta, dir -> dir }
                                        : RUN_BOLTZ_SM.out.intermediates.map { meta, dir -> dir }
            boltz_pdb       = is_prot_prot ? RUN_BOLTZ_LIGAND.out.pdb : RUN_BOLTZ_SM.out.pdb

            // --- prot-prot ---
            haddock_scores  = params.docking_tool == 'haddock3' ? RUN_HADDOCK3.out.scores : channel.empty()
            samplesheet     = is_prot_prot ? CANDIDATES_TO_COFOLD_SAMPLESHEET.out.samplesheet : channel.empty()
            versions        = ch_versions
}
