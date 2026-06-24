#!/usr/bin/env nextflow
nextflow.enable.dsl=2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/PhiPsiProt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/paudbr/PhiPsiProt
----------------------------------------------------------------------------------------
*/

include { PREPARE_ALPHAFOLD2_DBS } from './subworkflows/local/prepare_alphafold2_dbs'
include { PREPARE_ALPHAFOLD3_DBS } from './subworkflows/local/prepare_alphafold3_dbs'
include { PREPARE_ESMFOLD_DBS }    from './subworkflows/local/prepare_esmfold_dbs'
include { PREPARE_BOLTZ_DBS }      from './subworkflows/local/prepare_boltz_dbs'
include { CANDIDATES_SMILES_TO_AF3_FASTA } from './modules/local/candidates_smiles_to_af3_fasta/main'

include { ALPHAFOLD2 } from './workflows/alphafold2'
include { ALPHAFOLD3 } from './workflows/alphafold3'
include { ESMFOLD }    from './workflows/esmfold'
include { CHAI1 }      from './workflows/chai1'

include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_proteinfold_pipeline'
include { PIPELINE_COMPLETION }     from './subworkflows/local/utils_nfcore_proteinfold_pipeline'

include { POST_PROCESSING } from './subworkflows/local/post_processing'
include { METRICS_REPORT }  from './subworkflows/local/metrics_report'
include { DOCKING_COFOLDING } from './subworkflows/local/docking_cofolding'

include { SATURATION      } from './subworkflows/design_modes/saturation'
include { BACKBONE        } from './subworkflows/design_modes/backbone'
include { PEPTIDE_DESIGN  } from './subworkflows/design_modes/peptide_design'
include { ANTIBODY_DESIGN } from './subworkflows/design_modes/antibody_design'
include { COLLECT_AND_RANK } from './modules/local/collect_and_rank/main'

// Cascada de cribado por métricas estructurales (ESM -> Chai, ranking por YAML)
include { CASCADE } from './subworkflows/design_modes/cascade'

workflow NFCORE_PROTEINFOLD {

    take:
    samplesheet

    main:

    ch_samplesheet      = samplesheet
    ch_multiqc          = channel.empty()
    ch_versions         = channel.empty()
    ch_report_input     = channel.empty()
    ch_top_ranked_model = channel.empty()
    ch_design_candidates = channel.empty()

    def requested_modes = params.structural_tools.toLowerCase().split(",").collect { it.trim() }
    def requested_modes_size = requested_modes.size()

    def metrics_enabled = params.generate_metrics_report

    // ¿Corremos la cascada? Activa con --cascade y modo de diseño backbone.
    // En modo cascada NO corren los bloques paralelos af2/af3/esm/chai ni
    // METRICS_REPORT: la cascada hace ESM+Chai por dentro y saca su propio report.
    def run_cascade = params.cascade && params.mode == 'backbone'

    ch_dummy_file     = channel.value(file("$projectDir/assets/NO_FILE"))
    ch_dummy_file_pae = channel.value(file("$projectDir/assets/NO_FILE_PAE"))

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MODES (design)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (params.mode == 'screening') {
        SATURATION(file(params.input_pdb), params.target_chain, params.positions)
        ch_design_candidates = SATURATION.out.candidates
    }
    else if (params.mode == 'backbone') {
        BACKBONE(params.mode)
        ch_design_candidates = BACKBONE.out.candidates
    }
    else if (params.mode == 'peptide_design') {
        PEPTIDE_DESIGN(params.mode)
        ch_design_candidates = PEPTIDE_DESIGN.out.candidates
    }
    else if (params.mode == 'antibody_design') {
        ANTIBODY_DESIGN(params.mode)
        ANTIBODY_DESIGN.out.iptm_summary
            .subscribe { tsv ->
                log.info "[ANTIBODY_DESIGN] ipTM scores summary: ${tsv}"
            }
    }
    else if (params.mode != 'structural') {
        error "Unknown mode: ${params.mode}"
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CASCADE  (cribado estructural + ranking por YAML, sin AF3)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Resultado: ranking/esm/ranked.csv, ranking/chai/ranked.csv y
    cascade_report/PhiPsiProt_cascade_report.html
    */
    ch_design_candidates = params.input_csv ?
        channel.value(file(params.input_csv, checkIfExists: true)) :
        ch_design_candidates
    if (run_cascade) {
        CASCADE(
            ch_design_candidates,   // CSV de BACKBONE (id,sequence,...)
            ch_versions
        )
        ch_versions = ch_versions.mix(CASCADE.out.versions)
        // La cascada produce su propio ranking y report; no alimenta AF3.
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INTERMEDIATE LIGAND FILTER
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    else if (params.ligand) {
        if (!params.ligand_candidates_csv && params.mode == 'structural') {
            error "Ligand filtering needs --ligand_candidates_csv when --mode structural is used."
        }
        if (params.ligand && !params.receptor_sequence) {
            error "Ligand cofold mode needs --receptor_sequence"
        }
        if (!params.ligand_reference_pdb) {
            error "Ligand filtering needs --ligand_reference_pdb."Otros marcadores
            
        }
        if (!params.pocket_residues) {
            error "Ligand filtering needs --pocket_residues, for example A5,A4,A10,A238,A239,A240,A241,A242,A243,A100."
        }

        PREPARE_BOLTZ_DBS(
            params.boltz_db,
            params.boltz_ccd_path,
            params.boltz_model_path,
            params.boltz2_aff_path,
            params.boltz2_conf_path,
            params.boltz2_mols_path,
            params.boltz_ccd_link,
            params.boltz_model_link,
            params.boltz2_aff_link,
            params.boltz2_conf_link,
            params.boltz2_mols_link
        )
        ch_versions = ch_versions.mix(PREPARE_BOLTZ_DBS.out.versions)

        ch_ligand_candidates = params.ligand_candidates_csv ?
            channel.value(file(params.ligand_candidates_csv, checkIfExists: true)) :
            ch_design_candidates

        ch_ligand_reference = channel.value(file(params.ligand_reference_pdb, checkIfExists: true))
        ch_ligand_file = params.ligand_file ?
            channel.value(file(params.ligand_file, checkIfExists: true)) :
            ch_dummy_file

        def is_prot_prot = params.docking_tool in ['haddock3', 'rosettadock']

        if (is_prot_prot) {
            PREPARE_ESMFOLD_DBS(
                params.esmfold_db,
                params.esmfold_params_path,
                params.esmfold_3B_v1,
                params.esm2_t36_3B_UR50D,
                params.esm2_t36_3B_UR50D_contact_regression
            )
            ch_esmfold_params = PREPARE_ESMFOLD_DBS.out.params
        } else {
            ch_esmfold_params = channel.empty()
        }
        def receptor_seq_str = params.receptor_sequence
            ? "python3 ${projectDir}/bin/parse_fasta.py ${params.receptor_sequence}".execute().text.trim()
            : ''
        ch_receptor_seq = Channel.value(receptor_seq_str)

        DOCKING_COFOLDING(
            ch_ligand_candidates,
            ch_ligand_reference,
            ch_receptor_seq,
            ch_ligand_file,
            PREPARE_BOLTZ_DBS.out.boltz_model,
            PREPARE_BOLTZ_DBS.out.boltz_ccd,
            PREPARE_BOLTZ_DBS.out.boltz2_aff,
            PREPARE_BOLTZ_DBS.out.boltz2_conf,
            PREPARE_BOLTZ_DBS.out.boltz2_mols,
            ch_esmfold_params,
            params.esmfold_num_recycles
        )
        ch_versions = ch_versions.mix(DOCKING_COFOLDING.out.versions)

        DOCKING_COFOLDING.out.gnina_summary.view { "GNINA summary: $it" }
        DOCKING_COFOLDING.out.gnina_scores.view { "GNINA scores: $it" }

        if (is_prot_prot) {
            ch_samplesheet = DOCKING_COFOLDING.out.samplesheet
                .splitCsv(header: true)
                .map { row -> tuple([id: row.id], file(row.fasta)) }
        } else {
            ch_scoring_config = channel.value(
                file(params.scoring_config, checkIfExists: true)
            )

            COLLECT_AND_RANK(
                DOCKING_COFOLDING.out.gnina_scores,
                DOCKING_COFOLDING.out.boltz_dirs.collect(),
                ch_scoring_config,
                ch_ligand_candidates
            )
            ch_versions = ch_versions.mix(COLLECT_AND_RANK.out.versions)

            ch_lig = COLLECT_AND_RANK.out.ranked
                .splitCsv(header: true, sep: '\t')
                .filter { row -> row.in_top_n == 'yes' }
                .map { row -> [ [id: row.candidate_id], row.smiles ] }

            CANDIDATES_SMILES_TO_AF3_FASTA(ch_lig, ch_receptor_seq)

            ch_samplesheet = CANDIDATES_SMILES_TO_AF3_FASTA.out.fasta
            ch_versions = ch_versions.mix(CANDIDATES_SMILES_TO_AF3_FASTA.out.versions)
        }

    } else {
        ch_samplesheet = samplesheet  // el original
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ALPHAFOLD2   (no corre en modo cascada)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (!run_cascade && requested_modes.contains("af2")) {

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

        ch_report_input = ch_report_input.mix(
            ALPHAFOLD2.out.pdb
                .map { meta, files ->
                    def fileList = (files instanceof List) ? files.flatten() : [ files ]
                    [ meta, fileList.sort { a, b -> a.name <=> b.name }.take(5) ]
                }
                .join(ALPHAFOLD2.out.msa)
                .join(ALPHAFOLD2.out.pae)
        )
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ALPHAFOLD3   (no corre en modo cascada)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (!run_cascade && requested_modes.contains("af3") && params.mode == 'structural') {

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

        ch_report_input = ch_report_input.mix(
            ALPHAFOLD3.out.pdb
                .groupTuple()
                .map { meta, files ->
                    [ meta, files.flatten().sort { a, b -> a.name <=> b.name }.take(5) ]
                }
                .join(ALPHAFOLD3.out.msa)
                .join(ALPHAFOLD3.out.pae)
        )
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ESMFOLD   (no corre en modo cascada; la cascada lo hace por dentro)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (!run_cascade && requested_modes.contains("esm")) {

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

        ch_report_input = ch_report_input.mix(
            ESMFOLD.out.pdb
                .combine(ch_dummy_file)
                .combine(ch_dummy_file_pae)
        )
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHAI1   (no corre en modo cascada; la cascada lo hace por dentro)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (!run_cascade && requested_modes.contains("chai")) {

        CHAI1(
            ch_samplesheet,
            ch_versions
        )

        ch_multiqc = ch_multiqc.mix(CHAI1.out.multiqc_report)
        ch_versions = ch_versions.mix(CHAI1.out.versions)
        ch_top_ranked_model = ch_top_ranked_model.mix(CHAI1.out.top_ranked_pdb)

        ch_report_input = ch_report_input.mix(
            CHAI1.out.pdb
                .combine(ch_dummy_file)
                .combine(ch_dummy_file_pae)
        )
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    METRICS REPORT   (no corre en modo cascada; la cascada saca su propio report)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (metrics_enabled && !run_cascade && params.mode == 'structural') {

        ch_mr_af2_pdb  = requested_modes.contains("af2") ? ALPHAFOLD2.out.top_ranked_pdb : Channel.empty()
        ch_mr_af3_pdb  = requested_modes.contains("af3") ? ALPHAFOLD3.out.pdb : Channel.empty()
        ch_mr_esm_pdb  = requested_modes.contains("esm") ? ESMFOLD.out.pdb : Channel.empty()

        ch_mr_chai_plddt = requested_modes.contains("chai") ? CHAI1.out.plddt : Channel.empty()
        ch_mr_chai_pdb   = requested_modes.contains("chai") ? CHAI1.out.pdb : Channel.empty()
        ch_mr_chai_ptm   = requested_modes.contains("chai") ? CHAI1.out.ptms   : Channel.empty()
        ch_mr_chai_iptm  = requested_modes.contains("chai") ? CHAI1.out.iptms  : Channel.empty()
        ch_mr_chai_raw   = requested_modes.contains("chai") ? CHAI1.out.raw  : Channel.empty()

        ch_mr_af3_plddt = requested_modes.contains("af3") ? ALPHAFOLD3.out.plddt : Channel.empty()
        ch_mr_af3_ptm   = requested_modes.contains("af3") ? ALPHAFOLD3.out.ptms   : Channel.empty()
        ch_mr_af3_iptm  = requested_modes.contains("af3") ? ALPHAFOLD3.out.iptms  : Channel.empty()

        METRICS_REPORT(
            ch_mr_af2_pdb,
            ch_mr_af3_pdb,
            ch_mr_esm_pdb,
            ch_mr_chai_plddt,
            ch_mr_chai_ptm,
            ch_mr_chai_iptm,
            ch_mr_chai_pdb,
            ch_mr_chai_raw,
            ch_mr_af3_plddt,
            ch_mr_af3_ptm,
            ch_mr_af3_iptm
        )
        metrics_report_ch = METRICS_REPORT.out.report

    } else {
        metrics_report_ch = Channel.empty()
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    POSTPROCESSING   (no corre en modo cascada)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (!run_cascade) {
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
    }

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
