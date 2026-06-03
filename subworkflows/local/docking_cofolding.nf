include { RUN_GNINA_LIGAND_FILTER          } from '../../modules/local/run_gnina_ligand_filter/main'
include { CANDIDATES_TO_BOLTZ              } from '../../modules/local/candidates_to_boltz/main'
include { RUN_BOLTZ as RUN_BOLTZ_LIGAND    } from '../../modules/local/run_boltz/main'
include { BOLTZ_LIGAND_SUMMARY             } from '../../modules/local/boltz_ligand_summary/main'
include { CANDIDATES_TO_COFOLD_SAMPLESHEET } from '../../modules/local/candidates_to_cofold_samplesheet/main'
include { CANDIDATES_CSV_TO_FASTA } from '../../modules/local/candidates_csv_to_fasta/main'
include { ESMFOLD                          } from '../../workflows/esmfold'
include { RUN_HADDOCK3                     } from '../../modules/local/run_haddock3/main'
include { RUN_ROSETTADOCK                  } from '../../modules/local/run_rosettadock/main'

include {GENERATE_HADDOCK_AIRS             } from '../../modules/local/generate_haddock_airs/main'


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
    ch_versions  = channel.empty()
    def is_prot_prot = params.docking_tool in ['haddock3', 'rosettadock']

    // ──────────────────────────────────────────────────────────────
    // RAMA PROTEIN-PROTEIN: ESMFold + Boltz2 + HADDOCK3/RosettaDock
    // ──────────────────────────────────────────────────────────────
    if (is_prot_prot) {

        // Generar FASTAs — reutilizamos CANDIDATES_TO_BOLTZ para ambos
        CANDIDATES_TO_BOLTZ(candidates_csv, receptor_sequence)
        ch_versions = ch_versions.mix(CANDIDATES_TO_BOLTZ.out.versions)

        // Canal para Boltz2
        ch_boltz_input = CANDIDATES_TO_BOLTZ.out.fasta
            .flatMap { files ->
                def fileList = (files instanceof List) ? files.flatten() : [ files ]
                fileList.collect { fasta ->
                    def runId       = fasta.baseName
                    def candidateId = runId.replaceFirst(/_rep\d+$/, '')
                    [
                        [ id: runId, candidate_id: candidateId, model: "boltz2_ligand" ],
                        fasta,
                        []
                    ]
                }
            }

        


        //  REEMPLAZA POR ESTO:
        ch_esmfold_samplesheet = CANDIDATES_CSV_TO_FASTA(candidates_csv)
            .map { fasta -> tuple([id: fasta.baseName], fasta) }

        ESMFOLD(
            ch_esmfold_samplesheet,
            ch_versions,
            esmfold_params,
            esmfold_num_recycles
        )
        ch_versions = ch_versions.mix(ESMFOLD.out.versions)

        ch_folded_pdbs =  ESMFOLD.out.pdb
            .map { meta, pdb ->
                def candidate_id = meta.id
                [candidate_id, pdb]
            }

        // Boltz2 co-folding + affinity
        RUN_BOLTZ_LIGAND(
            ch_boltz_input,
            ch_boltz_model,
            ch_boltz_ccd,
            ch_boltz2_aff,
            ch_boltz2_conf,
            ch_boltz2_mols
        )
        ch_versions = ch_versions.mix(RUN_BOLTZ_LIGAND.out.versions)

        BOLTZ_LIGAND_SUMMARY(
            RUN_BOLTZ_LIGAND.out.confidence
                .map { meta, json -> json }
                .collect()
        )
        ch_versions = ch_versions.mix(BOLTZ_LIGAND_SUMMARY.out.versions)

        CANDIDATES_TO_COFOLD_SAMPLESHEET(
            BOLTZ_LIGAND_SUMMARY.out.summary,
            candidates_csv,
            receptor_sequence
        )
        ch_versions = ch_versions.mix(CANDIDATES_TO_COFOLD_SAMPLESHEET.out.versions)

        ch_haddock = ESMFOLD.out.pdb
            .map { meta, pdb ->
                def candidate_id = meta.id
                tuple(candidate_id, pdb, file(params.reference_pdb))
            }

        ch_air_input = ESMFOLD.out.pdb
            .map { meta, pdb ->
                def candidate_id = meta.id
                tuple(candidate_id, file(params.reference_pdb), pdb)
            }
        // Docking — when: en cada módulo controla cuál corre
        GENERATE_HADDOCK_AIRS(
            ch_air_input
        )

        ch_restraints = GENERATE_HADDOCK_AIRS.out
        .map { candidate_id, restraints, active, passive ->
            tuple(candidate_id, restraints)
        }

        ch_haddock_ready = ch_haddock
        .join(ch_restraints)
        .map { cand, pdb, rest ->
            tuple(cand[0], cand[1], cand[2], rest[1])
        }

        RUN_HADDOCK3(
            ch_haddock_ready
        )
        ch_versions = ch_versions.mix(RUN_HADDOCK3.out.versions)

        //RUN_ROSETTADOCK(
            //ch_folded_pdbs,
            //reference_pdb
        //)
        //ch_versions = ch_versions.mix(RUN_ROSETTADOCK.out.versions)

    // ──────────────────────────────────────────────────────────────
    // RAMA SMALL MOLECULE: solo GNINA, sin ESMFold ni Boltz2
    // ──────────────────────────────────────────────────────────────
    } else {

        RUN_GNINA_LIGAND_FILTER(
            candidates_csv,
            reference_pdb,
            ligand_file
        )
        ch_versions = ch_versions.mix(RUN_GNINA_LIGAND_FILTER.out.versions)
    }

    // ──────────────────────────────────────────────────────────────
    // Outputs — channel.empty() para la rama inactiva
    // ──────────────────────────────────────────────────────────────
    emit:
    gnina_scores      = is_prot_prot ? channel.empty() : RUN_GNINA_LIGAND_FILTER.out.scores
    gnina_summary     = is_prot_prot ? channel.empty() : RUN_GNINA_LIGAND_FILTER.out.summary
    boltz_scores      = is_prot_prot ? BOLTZ_LIGAND_SUMMARY.out.scores      : channel.empty()
    boltz_summary     = is_prot_prot ? BOLTZ_LIGAND_SUMMARY.out.summary     : channel.empty()
    boltz_pdb         = is_prot_prot ? RUN_BOLTZ_LIGAND.out.pdb             : channel.empty()
    haddock_scores    = params.docking_tool == 'haddock3'    ? RUN_HADDOCK3.out.scores    : channel.empty()
    //rosetta_scores    = params.docking_tool == 'rosettadock' ? RUN_ROSETTADOCK.out.scores : channel.empty()
    samplesheet       = is_prot_prot ? CANDIDATES_TO_COFOLD_SAMPLESHEET.out.samplesheet  : channel.empty()
    versions          = ch_versions
}