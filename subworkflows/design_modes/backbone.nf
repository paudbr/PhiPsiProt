/*
===============================================================================
PhiPsiProt - Mode 2: Backbone Stabilization
===============================================================================

Workflow:

STEP 1 Backbone design
STEP 2 Backbone summary
STEP 3 ProteinMPNN sequence design
STEP 4 Result aggregation
STEP 5 Biophysical annotation

Design strategy:
- BACKBONE_DESIGN generates candidate backbone structures using RFdiffusion.
- ProteinMPNN designs one or more sequences compatible with each backbone.
- BIOPHYSICAL_ANNOTATION adds developability descriptors.
- No candidates are removed at this stage.

Output:
- designed backbone PDB files
- backbone summary
- ProteinMPNN-designed sequences
- merged backbone results
- biophysically annotated candidates
- final_backbone_results

Notes:
- Mode 2 focuses on structural stabilization and controlled backbone redesign.
- RFdiffusion is used in partial diffusion mode.
- Ligand-aware design is handled in Mode 3 using LigandMPNN.
===============================================================================
*/

include { BACKBONE_DESIGN } from '../../modules/local/backbone_design/main'
include { BACKBONE_SUMMARY } from '../../modules/local/backbone_summary/main'
include { PROTEINMPNN_DESIGN } from '../../modules/local/proteinmpnn_design/main'
include { MERGE_BACKBONE_RESULTS } from '../../modules/local/merge_backbone_results/main'
include { BIOPHYSICAL_ANNOTATION } from '../../modules/local/biophysical_annotation/main'
include { BACKBONE_REPORT } from '../../modules/local/backbone_report/main'


workflow BACKBONE {

    take:
    input_pdb

    main:

    /*
     * STEP 1
     * Generate candidate backbone structures from the input PDB.
     */

    BACKBONE_DESIGN(input_pdb)

    /*
     * Flatten the backbone design output so that each designed PDB can be
     * processed independently by ProteinMPNN.
     */

    designed_pdb_ch = BACKBONE_DESIGN.out[0].flatten()

    /*
     * STEP 2
     * Summarize the generated backbone candidates.
     */

    BACKBONE_SUMMARY(BACKBONE_DESIGN.out[0])

    /*
     * STEP 3
     * Design sequences for each generated backbone using ProteinMPNN.
     */

    PROTEINMPNN_DESIGN(designed_pdb_ch)

    /*
     * STEP 4
     * Merge backbone metadata with ProteinMPNN sequence-design results.
     */

    MERGE_BACKBONE_RESULTS(
        BACKBONE_SUMMARY.out[0],
        PROTEINMPNN_DESIGN.out.collect()
    )

    /*
     * STEP 5
     * Annotate candidates with sequence-based biophysical descriptors.
     * No candidates are removed here.
     */

    BIOPHYSICAL_ANNOTATION(
        MERGE_BACKBONE_RESULTS.out
    )

    

    BACKBONE_REPORT(
        BIOPHYSICAL_ANNOTATION.out,
        input_pdb,
        BACKBONE_DESIGN.out[0]
    )

    final_results_ch = BIOPHYSICAL_ANNOTATION.out

    /*
     * Expose the final Mode 2 output.
     */

    emit:
    final_backbone_results = final_results_ch
    backbone_report = BACKBONE_REPORT.out
}