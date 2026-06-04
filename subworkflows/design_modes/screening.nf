/*
===============================================================================
PhiPsiProt - Mode 1: Mutational Screening
===============================================================================

Workflow:

1. Residue selection
2. PyRosetta mutational scan
3. ddG filtering
4. Candidate ranking
5. Plot generation
6. Result aggregation
7. HTML report generation

Selection modes:
- manual
- pocket
- auto_pocket
- interface
- all

Output:
- candidates.csv
- ranked_candidates.csv
- plots/
- screening_report.html

===============================================================================
*/

include { PYROSETTA_SCREENING } from '../../modules/local/pyrosetta_screening/main'
include { DDG_FILTER } from '../../modules/local/ddg_filter/main'
include { RESIDUE_SELECTION } from '../../modules/local/residue_selection/main'
include { AUTO_POCKET_DETECTION } from '../../modules/local/auto_pocket_detection/main'
include { RANK_SCREENING } from '../../modules/local/ranking/main'
include { BIOPHYSICAL_FILTER } from '../../modules/local/biophysical_filter/main'
include { MERGE_SCREENING_RESULTS } from '../../modules/local/merge_screening_results/main'
include { SCREENING_PLOTS } from '../../modules/local/screening_plots/main'
include { SCREENING_REPORT } from '../../modules/local/screening_report/main'

workflow SCREENING {

    take:
    input_pdb
    target_chain
    positions

    main:

    /*
     * STEP 1
     * Select mutable residues according to:
     * - manual positions
     * - ligand-defined pocket
     * - automatically predicted pocket using fpocket
     * - protein-protein interface
     * - all residues in the target chain
     */

    if ((params.selection_mode ?: 'manual') == 'auto_pocket') {

        AUTO_POCKET_DETECTION(
            input_pdb,
            params.target_chain ?: 'A',
            params.pocket_rank ?: 1
        )

        selected_positions_ch = AUTO_POCKET_DETECTION.out[0]
        selection_summary_ch = AUTO_POCKET_DETECTION.out[1]

    } else {

        RESIDUE_SELECTION(
            input_pdb,
            params.target_chain ?: 'A',
            params.positions ?: 'ALL',
            params.selection_mode ?: 'manual',
            params.ligand_resname ?: '',
            params.interface_chain ?: '',
            params.distance_cutoff ?: 6.0
        )

        selected_positions_ch = RESIDUE_SELECTION.out[0]
        selection_summary_ch = RESIDUE_SELECTION.out[1]
    }

    /*
     * STEP 2
     * Generate all single amino acid substitutions
     * and calculate ΔΔG using PyRosetta
     */

    PYROSETTA_SCREENING(
        input_pdb,
        params.target_chain ?: 'A',
        selected_positions_ch
    )

    /*
     * STEP 3
     * Keep stabilizing mutations according
     * to the configured ddG threshold
     */

    DDG_FILTER(
        PYROSETTA_SCREENING.out[0]
    )

    BIOPHYSICAL_FILTER(
        PYROSETTA_SCREENING.out[0]
    )


    /*
     * STEP 4
     * Rank candidates by ddG score
     */

    RANK_SCREENING(
        BIOPHYSICAL_FILTER.out
    )

    /*
     * STEP 5
     * Generate visual summaries
     * (heatmap, ranking plot, top mutations)
     */

    SCREENING_PLOTS(
        RANK_SCREENING.out
    )

    /*
     * STEP 6
     * Merge ranked candidates with selection metadata
     */

    MERGE_SCREENING_RESULTS(
        RANK_SCREENING.out,
        selection_summary_ch
    )

    /*
     * STEP 7
     * Generate final HTML report
     */

    SCREENING_REPORT(
        MERGE_SCREENING_RESULTS.out,
        SCREENING_PLOTS.out.collect()
    )

    emit:
    selected_positions = selected_positions_ch
    candidates_csv = PYROSETTA_SCREENING.out[0]
    candidates_fasta = PYROSETTA_SCREENING.out[1]
    filtered_candidates_csv = DDG_FILTER.out[0]
    filtered_candidates_fasta = DDG_FILTER.out[1]
    screening_samplesheet = DDG_FILTER.out[2]
    selection_summary = selection_summary_ch
    ranked_candidates = RANK_SCREENING.out
    final_screening_results = MERGE_SCREENING_RESULTS.out
}