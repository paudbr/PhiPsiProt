#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/proteinfold
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/proteinfold
    Website: https://nf-co.re/proteinfold
    Slack  : https://nfcore.slack.com/channels/proteinfold
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_ALPHAFOLD2_DBS           } from './subworkflows/local/prepare_alphafold2_dbs'
include { PREPARE_ALPHAFOLD3_DBS           } from './subworkflows/local/prepare_alphafold3_dbs'
include { PREPARE_ESMFOLD_DBS              } from './subworkflows/local/prepare_esmfold_dbs'
include { ALPHAFOLD2                       } from './workflows/alphafold2'
include { CHAI1                            } from './workflows/chai1'
include { ALPHAFOLD3                       } from './workflows/alphafold3'
include { ESMFOLD                          } from './workflows/esmfold'

include { PIPELINE_INITIALISATION          } from './subworkflows/local/utils_nfcore_proteinfold_pipeline'
include { PIPELINE_COMPLETION              } from './subworkflows/local/utils_nfcore_proteinfold_pipeline'
include { getColabfoldAlphafold2Params     } from './subworkflows/local/utils_nfcore_proteinfold_pipeline'
include { getColabfoldAlphafold2ParamsPath } from './subworkflows/local/utils_nfcore_proteinfold_pipeline'
include { POST_PROCESSING                  } from './subworkflows/local/post_processing'
include { SATURATION } from './subworkflows/design_modes/saturation'
include { BACKBONE } from './subworkflows/design_modes/backbone'
include { PEPTIDE_DESIGN } from './subworkflows/design_modes/peptide_design'
include { ANTIBODY_DESIGN } from './subworkflows/design_modes/antibody_design'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLABFOLD PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.colabfold_alphafold2_params_link = getColabfoldAlphafold2Params()
params.colabfold_alphafold2_params_path = getColabfoldAlphafold2ParamsPath()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline
//

workflow NFCORE_PROTEINFOLD {

    take:
    samplesheet  // channel: samplesheet read in from --input

    main:
    ch_samplesheet       = samplesheet
    ch_multiqc           = channel.empty()
    ch_versions          = channel.empty()
    ch_report_input      = channel.empty()
    ch_top_ranked_model  = channel.empty()
    requested_modes      = params.structural_tools.toLowerCase().split(",").collect { it.trim() }
    requested_modes_size = requested_modes.size()

    ch_dummy_file = channel.fromPath("$projectDir/assets/NO_FILE")
    ch_dummy_file_pae = channel.fromPath("$projectDir/assets/NO_FILE_PAE")
    
    if (params.mode == 'saturation') {
    SATURATION(file(params.input_pdb), params.target_chain, params.positions)
    }

    if (params.mode == 'backbone') {
    BACKBONE(params.mode)
    }

    if (params.mode == 'peptide_design') {
    PEPTIDE_DESIGN(params.mode)
    }

    if (params.mode == 'antibody_design') {
    ANTIBODY_DESIGN(params.mode)
    }

    else if (params.mode != 'structural') {
        error "Unknown mode: ${params.mode}"
    }
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    NFCORE_PROTEINFOLD(
        PIPELINE_INITIALISATION.out.samplesheet
    )

    //
    // WORKFLOW: Run alphafold2
    //
    if(requested_modes.contains("af2")) {

        //
        // SUBWORKFLOW: Prepare Alphafold2 DBs
        //
        PREPARE_ALPHAFOLD2_DBS (
            params.alphafold2_db,
            params.alphafold2_full_dbs,
            params.alphafold2_bfd_path,
            params.alphafold2_small_bfd_path,
            params.alphafold2_params_path,
            params.alphafold2_mgnify_path,
            params.alphafold2_pdb70_path,
            params.alphafold2_pdb_mmcif_path,
            params.alphafold2_pdb_obsolete_path,
            params.alphafold2_uniref30_path,
            params.alphafold2_uniref90_path,
            params.alphafold2_pdb_seqres_path,
            params.alphafold2_uniprot_path,
            params.alphafold2_bfd_link,
            params.alphafold2_small_bfd_link,
            params.alphafold2_params_link,
            params.alphafold2_mgnify_link,
            params.alphafold2_pdb70_link,
            params.alphafold2_pdb_mmcif_link,
            params.alphafold2_pdb_obsolete_link,
            params.alphafold2_uniref30_link,
            params.alphafold2_uniref90_link,
            params.alphafold2_pdb_seqres_link,
            params.alphafold2_uniprot_sprot_link,
            params.alphafold2_uniprot_trembl_link
        )
        ch_versions = ch_versions.mix(PREPARE_ALPHAFOLD2_DBS.out.versions)

        //
        // WORKFLOW: Run nf-core/alphafold2 workflow
        //
        ALPHAFOLD2 (
            ch_samplesheet,
            ch_versions,
            params.alphafold2_full_dbs,
            params.alphafold2_mode,
            params.alphafold2_model_preset,
            params.uniref30_prefix,
            PREPARE_ALPHAFOLD2_DBS.out.params,
            PREPARE_ALPHAFOLD2_DBS.out.bfd,
            PREPARE_ALPHAFOLD2_DBS.out.small_bfd,
            PREPARE_ALPHAFOLD2_DBS.out.mgnify,
            PREPARE_ALPHAFOLD2_DBS.out.pdb70,
            PREPARE_ALPHAFOLD2_DBS.out.pdb_mmcif,
            PREPARE_ALPHAFOLD2_DBS.out.pdb_obsolete,
            PREPARE_ALPHAFOLD2_DBS.out.uniref30,
            PREPARE_ALPHAFOLD2_DBS.out.uniref90,
            PREPARE_ALPHAFOLD2_DBS.out.pdb_seqres,
            PREPARE_ALPHAFOLD2_DBS.out.uniprot
        )
        ch_multiqc          = ch_multiqc.mix(ALPHAFOLD2.out.multiqc_report.collect())
        ch_versions         = ch_versions.mix(ALPHAFOLD2.out.versions)
        ch_report_input     = ch_report_input
                                .mix(ALPHAFOLD2
                                .out
                                .pdb
                                .map { it ->
                                    [ it[0],
                                        it[1].sort { path ->
                                            def filename = path.name
                                            def matcher = filename =~ /ranked_(\d+)\.pdb/
                                            if (matcher.matches()) {
                                                return matcher[0][1].toInteger()
                                            } else {
                                                return 0  // fallback if no match
                                            }
                                        }.subList(0, Math.min(5, it[1].size() as int))
                                    ]
                                }
                                .join(ALPHAFOLD2.out.msa)
                                .join(ALPHAFOLD2.out.pae)
                            )

        ch_top_ranked_model = ch_top_ranked_model.mix(ALPHAFOLD2.out.top_ranked_pdb)
    }

    //
    // WORKFLOW: Run alphafold3
    //
    if(requested_modes.contains("af3")) {

        //
        // SUBWORKFLOW: Prepare Alphafold3 DBs
        //
        PREPARE_ALPHAFOLD3_DBS (
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

        //
        // WORKFLOW: Run nf-core/alphafold3 workflow
        //
        ALPHAFOLD3 (
            ch_samplesheet,
            ch_versions,
            PREPARE_ALPHAFOLD3_DBS.out.params,
            PREPARE_ALPHAFOLD3_DBS.out.small_bfd,
            PREPARE_ALPHAFOLD3_DBS.out.mgnify,
            PREPARE_ALPHAFOLD3_DBS.out.pdb_mmcif,
            PREPARE_ALPHAFOLD3_DBS.out.uniref90,
            PREPARE_ALPHAFOLD3_DBS.out.pdb_seqres,
            PREPARE_ALPHAFOLD3_DBS.out.uniprot
        )

        ch_multiqc      = ch_multiqc.mix(ALPHAFOLD3.out.multiqc_report)
        ch_versions     = ch_versions.mix(ALPHAFOLD3.out.versions)
        ch_report_input = ch_report_input
                            .mix(
                                ALPHAFOLD3
                                    .out
                                    .pdb
                                    .map { it ->
                                        [
                                            it[0],
                                            it[1].sort { path ->
                                                def filename = path.name
                                                def matcher = filename =~ /.*_ranked_(\d+)\.pdb/
                                                if (matcher.matches()) {
                                                    return matcher[0][1].toInteger()
                                                } else {
                                                    return 0  // fallback if no match
                                                }
                                            }.subList(0, Math.min(5, it[1].size() as int))
                                        ]
                                    }
                                .join(ALPHAFOLD3.out.msa)
                                .join(ALPHAFOLD3.out.pae)
                            )
        ch_top_ranked_model = ch_top_ranked_model.mix(ALPHAFOLD3.out.top_ranked_pdb)
    }


    //
    // WORKFLOW: Run esmfold
    //
    if(requested_modes.contains("esm")) {

        //
        // SUBWORKFLOW: Prepare esmfold DBs
        //
        PREPARE_ESMFOLD_DBS (
            params.esmfold_db,
            params.esmfold_params_path,
            params.esmfold_3B_v1,
            params.esm2_t36_3B_UR50D,
            params.esm2_t36_3B_UR50D_contact_regression
        )
        ch_versions = ch_versions.mix(PREPARE_ESMFOLD_DBS.out.versions)

        //
        // WORKFLOW: Run nf-core/esmfold workflow
        //
        ESMFOLD (
            ch_samplesheet,
            ch_versions,
            PREPARE_ESMFOLD_DBS.out.params,
            params.esmfold_num_recycles
        )

        ch_multiqc      = ch_multiqc.mix(ESMFOLD.out.multiqc_report.collect())
        ch_versions     = ch_versions.mix(ESMFOLD.out.versions)
        ch_report_input = ch_report_input.mix(
            ESMFOLD.out.pdb
                .combine(ch_dummy_file)
                .combine(ch_dummy_file_pae)
        )
        ch_top_ranked_model = ch_top_ranked_model.mix(ESMFOLD.out.pdb)
    }


    //
    // WORKFLOW: Run Chai-1
    //
    if(requested_modes.contains("chai")) {
        //
        // WORKFLOW: Run nf-core/chai1 workflow
        //
        CHAI1 (
            ch_samplesheet,
            ch_versions
        )
        ch_multiqc      = ch_multiqc.mix(CHAI1.out.multiqc_report)
        ch_versions     = ch_versions.mix(CHAI1.out.versions)
        ch_report_input = ch_report_input.mix(
            CHAI1.out.pdb
                .combine(ch_dummy_file)
                .combine(ch_dummy_file_pae)
        )
        ch_top_ranked_model = ch_top_ranked_model.mix(CHAI1.out.top_ranked_pdb)
    }
     
    //
    // POST PROCESSING: generate visualisation reports
    //
    ch_multiqc_config        = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true).first()
    ch_multiqc_custom_config = params.multiqc_config ? channel.fromPath( params.multiqc_config ).first()  : channel.empty()
    ch_multiqc_logo          = params.multiqc_logo   ? channel.fromPath( params.multiqc_logo ).first()    : channel.empty()
    ch_multiqc_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_report_template     = channel.value(file("$projectDir/assets/report_template.html", checkIfExists: true))
    ch_comparison_template = channel.value(file("$projectDir/assets/comparison_template.html", checkIfExists: true))

    POST_PROCESSING(
        params.skip_visualisation,
        requested_modes_size,
        ch_report_input,
        ch_report_template,
        ch_comparison_template,
        params.skip_foldseek,
        params.foldseek_db,
        params.foldseek_db_path,
        params.skip_multiqc,
        params.outdir,
        ch_versions,
        ch_multiqc,
        ch_multiqc_config,
        ch_multiqc_custom_config,
        ch_multiqc_logo,
        ch_multiqc_methods_description,
        ch_top_ranked_model
    )

    emit:
    multiqc_report = ch_multiqc
}
