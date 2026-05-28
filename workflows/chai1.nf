/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { RUN_CHAI1                         } from '../modules/local/run_chai1'
include { MMCIF2PDB as MMCIF2PDB_TOP_RANKED } from '../modules/local/mmcif2pdb/main.nf'
include { MMCIF2PDB as MMCIF2PDB_MODELS     } from '../modules/local/mmcif2pdb/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELPER PROCESS: Download Chai-1 weights once and cache them via storeDir
    storeDir makes Nextflow skip this process if the output already exists on disk
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CHAI1 {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_versions    // channel: [ path(versions.yml) ]

    main:
    ch_pdb_final      = channel.empty()
    ch_top_ranked_pdb = channel.empty()
    ch_multiqc_report = channel.empty()
    ch_ptms_final     = channel.empty()

    //
    // Download weights once — storeDir skips this automatically on reruns
    // workflow
    ch_weights = Channel.value(file(params.chai1_weights_path, checkIfExists: true))
    //

    //
    // MODULE: Run Chai-1 — no FASTA→JSON conversion needed, takes FASTA directly
    //
    RUN_CHAI1 (
        ch_samplesheet,
        ch_weights
    )
    ch_versions = ch_versions.mix(RUN_CHAI1.out.versions)

    // Convert all ranked mmcifs to pdb
    MMCIF2PDB_MODELS (
        RUN_CHAI1
            .out
            .cif
            .groupTuple()
            .map { meta, files ->
                [ meta, files.flatten() ]
            }
    )
    ch_versions = ch_versions.mix(MMCIF2PDB_MODELS.out.versions)

    MMCIF2PDB_MODELS
        .out
        .pdb
        .map { it ->
            def meta   = it[0].clone()
            meta.model = "chai1"
            def files  = (it[1] instanceof List) ? it[1] : [ it[1] ]
            [ meta, files ]
        }
        .set { ch_pdb_final }

    // Convert top-ranked mmcif to pdb
    MMCIF2PDB_TOP_RANKED (
        RUN_CHAI1.out.top_ranked_cif
    )
    ch_versions = ch_versions.mix(MMCIF2PDB_TOP_RANKED.out.versions)

    MMCIF2PDB_TOP_RANKED
        .out
        .pdb
        .map { it ->
            def meta   = it[0].clone()
            meta.model = "chai1"
            [ meta, it[1] ]
        }
        .set { ch_top_ranked_pdb }


    RUN_CHAI1
        .out
        .ptms
        .map { it ->
            def meta   = it[0].clone()
            meta.model = "chai1"
            [ meta, it[1] ]
        }
        .set { ch_ptms_final }

    // Prepare multiqc report input
    RUN_CHAI1
        .out
        .multiqc
        .map { it -> it[1] }
        .toSortedList()
        .map { it ->
            [ [ "model": "chai1" ], it.flatten() ]
        }
        .set { ch_multiqc_report }

    emit:
    top_ranked_pdb = ch_top_ranked_pdb // channel: [ meta, /path/to/*.pdb ]
    pdb            = ch_pdb_final      // channel: [ meta, /path/to/*.pdb, ..., /path/to/*.pdb ]
    ptms           = ch_ptms_final     // channel: [ meta, /path/to/*_ptm.tsv ] — replaces msa + pae
    multiqc_report = ch_multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions       // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/