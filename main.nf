#!/usr/bin/env nextflow
nextflow.enable.dsl=2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/PhiPsiProt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/paudbr/PhiPsiProt
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { PREPARE_ALPHAFOLD2_DBS } from './subworkflows/local/prepare_alphafold2_dbs'
include { PREPARE_ALPHAFOLD3_DBS } from './subworkflows/local/prepare_alphafold3_dbs'
include { PREPARE_ESMFOLD_DBS }    from './subworkflows/local/prepare_esmfold_dbs'

include { ALPHAFOLD2 } from './workflows/alphafold2'
include { ALPHAFOLD3 } from './workflows/alphafold3'
include { ESMFOLD }    from './workflows/esmfold'
include { CHAI1 }      from './workflows/chai1'

include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_proteinfold_pipeline'
include { PIPELINE_COMPLETION }     from './subworkflows/local/utils_nfcore_proteinfold_pipeline'

include { POST_PROCESSING } from './subworkflows/local/post_processing'
include { METRICS_REPORT }  from './subworkflows/local/metrics_report'

workflow NFCORE_PROTEINFOLD {

    take:
    samplesheet

    main:

    ch_samplesheet      = samplesheet
    ch_multiqc          = channel.empty()
    ch_versions         = channel.empty()
    ch_report_input     = channel.empty()
    ch_top_ranked_model = channel.empty()

    def requested_modes = params.structural_tools.toLowerCase().split(",").collect { it.trim() }
    def requested_modes_size = requested_modes.size()

    def metrics_enabled = params.generate_metrics_report

    ch_dummy_file     = channel.fromPath("$projectDir/assets/NO_FILE")
    ch_dummy_file_pae = channel.fromPath("$projectDir/assets/NO_FILE_PAE")

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MODES (design)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (params.mode == 'screening') {
        SATURATION(file(params.input_pdb), params.target_chain, params.positions)
    }

    else if (params.mode == 'backbone') {
        BACKBONE(params.mode)
    }

    else if (params.mode == 'peptide_design') {
        PEPTIDE_DESIGN(params.mode)
    }

    else if (params.mode == 'antibody_design') {
        ANTIBODY_DESIGN(params.mode)
    }

    else if (params.mode != 'structural') {
        error "Unknown mode: ${params.mode}"
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ALPHAFOLD2
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (requested_modes.contains("af2")) {

        PREPARE_ALPHAFOLD2_DBS(
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

        ALPHAFOLD2(
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

        ch_multiqc = ch_multiqc.mix(ALPHAFOLD2.out.multiqc_report.collect())
        ch_versions = ch_versions.mix(ALPHAFOLD2.out.versions)

        ch_top_ranked_model = ch_top_ranked_model.mix(ALPHAFOLD2.out.top_ranked_pdb)
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ALPHAFOLD3
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (requested_modes.contains("af3")) {

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

        ALPHAFOLD3(
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

        ch_multiqc = ch_multiqc.mix(ALPHAFOLD3.out.multiqc_report)
        ch_versions = ch_versions.mix(ALPHAFOLD3.out.versions)

        ch_top_ranked_model = ch_top_ranked_model.mix(ALPHAFOLD3.out.top_ranked_pdb)
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ESMFOLD
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (requested_modes.contains("esm")) {

        PREPARE_ESMFOLD_DBS(
            params.esmfold_db,
            params.esmfold_params_path,
            params.esmfold_3B_v1,
            params.esm2_t36_3B_UR50D,
            params.esm2_t36_3B_UR50D_contact_regression
        )

        ESMFOLD(
            ch_samplesheet,
            ch_versions,
            PREPARE_ESMFOLD_DBS.out.params,
            params.esmfold_num_recycles
        )

        ch_multiqc = ch_multiqc.mix(ESMFOLD.out.multiqc_report.collect())
        ch_versions = ch_versions.mix(ESMFOLD.out.versions)

        ch_top_ranked_model = ch_top_ranked_model.mix(ESMFOLD.out.pdb)
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHAI1
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (requested_modes.contains("chai")) {

        CHAI1(
            ch_samplesheet,
            ch_versions
        )

        ch_multiqc = ch_multiqc.mix(CHAI1.out.multiqc_report)
        ch_versions = ch_versions.mix(CHAI1.out.versions)

        ch_top_ranked_model = ch_top_ranked_model.mix(CHAI1.out.top_ranked_pdb)
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    METRICS REPORT (FIXED)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (metrics_enabled) {

        ch_mr_af2_pdb  = requested_modes.contains("af2") ? ALPHAFOLD2.out.top_ranked_pdb : Channel.empty()

        ch_mr_af3_pdb = requested_modes.contains("af3") ? ALPHAFOLD3.out.pdb : Channel.empty()


        ch_mr_esm_pdb  = requested_modes.contains("esm") ? ESMFOLD.out.pdb : Channel.empty()

        ch_mr_chai_plddt = requested_modes.contains("chai") ? CHAI1.out.plddt : Channel.empty()
        ch_mr_chai_pdb = requested_modes.contains("chai") ? CHAI1.out.pdb : Channel.empty()
        ch_mr_chai_ptm   = requested_modes.contains("chai") ? CHAI1.out.ptms   : Channel.empty()
        ch_mr_chai_iptm  = requested_modes.contains("chai") ? CHAI1.out.iptms  : Channel.empty()

        ch_mr_chai_raw   = requested_modes.contains("chai") ? CHAI1.out.raw  : Channel.empty() 
        ch_mr_af3_plddt = requested_modes.contains("af3") ? ALPHAFOLD3.out.plddt : Channel.empty() 


        METRICS_REPORT(
            ch_mr_af2_pdb,
            ch_mr_af3_pdb,
            ch_mr_esm_pdb,
            ch_mr_chai_pdb,   
            ch_mr_chai_plddt,
            ch_mr_chai_ptm,
            ch_mr_chai_iptm,
            ch_mr_chai_raw,
            ch_mr_af3_plddt
        )
        metrics_report_ch = METRICS_REPORT.out.report

    } else {
        metrics_report_ch = Channel.empty()
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    POSTPROCESSING
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_multiqc_config = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true).first()
    ch_multiqc_custom_config = params.multiqc_config ? channel.fromPath(params.multiqc_config).first() : channel.empty()
    ch_multiqc_logo = params.multiqc_logo ? channel.fromPath(params.multiqc_logo).first() : channel.empty()

    ch_multiqc_methods_description =
        params.multiqc_methods_description ?
            file(params.multiqc_methods_description, checkIfExists: true) :
            file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

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
        multiqc_report  = ch_multiqc
        metrics_report  = metrics_report_ch
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ENTRY WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    PIPELINE_INITIALISATION(
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

    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_PROTEINFOLD.out.multiqc_report
    )
}