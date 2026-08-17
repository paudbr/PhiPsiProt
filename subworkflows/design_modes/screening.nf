/*
===============================================================================
PhiPsiProt - Mode 1: Mutational Screening
===============================================================================

Workflow:

STEP 1 Residue selection
STEP 2 PyRosetta screening
STEP 3 DDG annotation
STEP 4 Biophysical annotation
STEP 5 Candidate ranking
STEP 6 Plot generation
STEP 7 Result aggregation
STEP 8 HTML report generation

Selection modes:
- manual
- pocket
- auto_pocket
- interface
- all

Output:
- candidates.csv
- candidates_ddg_annotated.csv
- candidates_biophysical_annotated.csv
- ranked_candidates.csv
- plots/
- screening_report.html

===============================================================================
*/

include { PYROSETTA_SCREENING } from '../../modules/local/pyrosetta_screening/main'
include { DDG_ANNOTATION } from '../../modules/local/ddg_annotation/main'
include { RESIDUE_SELECTION } from '../../modules/local/residue_selection/main'
include { AUTO_POCKET_DETECTION } from '../../modules/local/auto_pocket_detection/main'
include { RANK_SCREENING } from '../../modules/local/ranking/main'
include { BIOPHYSICAL_ANNOTATION } from '../../modules/local/biophysical_annotation/main'
include { MERGE_SCREENING_RESULTS } from '../../modules/local/merge_screening_results/main'
include { SCREENING_PLOTS } from '../../modules/local/screening_plots/main'
include { SCREENING_REPORT } from '../../modules/local/screening_report/main'
include { PYMOL_SCREENING_IMAGES } from '../../modules/local/pymol_screening_images/main'

workflow SCREENING {

    take:
    input_pdb
    target_chain
    positions

    main:

    /*
     * Define local workflow variables from the inputs received via `take`.
     * This avoids mixing global params with subworkflow inputs.
     */
    def chain = target_chain ?: 'A'
    def pos = positions ?: 'ALL'
    def mode = params.selection_mode ?: 'manual'
    def cutoff = params.distance_cutoff ?: 6.0

    /*
     * STEP 1
     * Select mutable residues according to:
     * - manual positions
     * - ligand-defined pocket
     * - automatically predicted pocket using fpocket
     * - protein-protein interface
     * - all residues in the target chain
     */

    if (mode == 'auto_pocket') {

        AUTO_POCKET_DETECTION(
            input_pdb,
            chain,
            params.pocket_rank ?: 1
        )

        selected_positions_ch = AUTO_POCKET_DETECTION.out[0]
        selection_summary_ch = AUTO_POCKET_DETECTION.out[1]

    } else {

        RESIDUE_SELECTION(
            input_pdb,
            chain,
            pos,
            mode,
            params.ligand_resname ?: '',
            params.interface_chain ?: '',
            cutoff
        )

        selected_positions_ch = RESIDUE_SELECTION.out[0]
        selection_summary_ch = RESIDUE_SELECTION.out[1]
    }

    /*
     * STEP 2
     * Generate all single amino-acid substitutions
     * and calculate ΔΔG using PyRosetta.
     */

    PYROSETTA_SCREENING(
        input_pdb,
        chain,
        selected_positions_ch
    )

    /*
     * STEP 3
     * Annotate candidates with ddG threshold information.
     * No candidates are removed at this stage.
     */

    DDG_ANNOTATION(
        PYROSETTA_SCREENING.out[0]
    )

    /*
     * STEP 4
     * Annotate candidates with sequence-based developability descriptors.
     * No candidates are removed at this stage.
     */

    BIOPHYSICAL_ANNOTATION(
        DDG_ANNOTATION.out
    )

    /*
     * STEP 5
     * Rank candidates using a combined score based on PyRosetta ddG
     * and biophysical developability annotations.
     */

    RANK_SCREENING(
        BIOPHYSICAL_ANNOTATION.out
    )

    /*
     * STEP 6
     * Generate visual summaries:
     * heatmap, ranking plot and top mutation plots.
     */

    SCREENING_PLOTS(
        RANK_SCREENING.out
    )

    /*
     * STEP 7
     * Merge ranked candidates with residue selection metadata.
     */

    MERGE_SCREENING_RESULTS(
        RANK_SCREENING.out,
        selection_summary_ch
    )

    PYMOL_SCREENING_IMAGES(
        input_pdb,
        MERGE_SCREENING_RESULTS.out
    )

    /*
     * STEP 8
     * Generate final HTML report.
     */

    SCREENING_REPORT(
        MERGE_SCREENING_RESULTS.out,
        SCREENING_PLOTS.out.collect(),
        PYMOL_SCREENING_IMAGES.out.images.collect()
    )

    /*
     * Expose intermediate and final outputs so that the main workflow
     * can publish, test, or reuse them in downstream design modes.
     */

    emit:
    selected_positions = selected_positions_ch
    candidates_csv = PYROSETTA_SCREENING.out[0]
    candidates_fasta = PYROSETTA_SCREENING.out[1]
    ddg_annotated_candidates = DDG_ANNOTATION.out
    biophysical_annotated_candidates = BIOPHYSICAL_ANNOTATION.out
    selection_summary = selection_summary_ch
    ranked_candidates = RANK_SCREENING.out
    final_screening_results = MERGE_SCREENING_RESULTS.out
    pymol_images = PYMOL_SCREENING_IMAGES.out.images
}