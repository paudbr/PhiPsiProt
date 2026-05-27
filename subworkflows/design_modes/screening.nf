include { PYROSETTA_SCREENING } from '../../modules/local/pyrosetta_screening/main'
include { DDG_FILTER } from '../../modules/local/ddg_filter/main'
include { RESIDUE_SELECTION } from '../../modules/local/residue_selection/main'
include { RANK_SCREENING } from '../../modules/local/ranking/main'

workflow SCREENING {

    take:
    input_pdb
    target_chain
    positions

    main:

    RESIDUE_SELECTION(
        input_pdb,
        target_chain,
        positions,
        params.selection_mode ?: 'manual',
        params.ligand_resname ?: '',
        params.interface_chain ?: '',
        params.distance_cutoff ?: 6.0
    )

    PYROSETTA_SCREENING(
        input_pdb,
        target_chain,
        RESIDUE_SELECTION.out[0]
    )

    DDG_FILTER(
        PYROSETTA_SCREENING.out[0]
    )

    RANK_SCREENING(
    DDG_FILTER.out[0]
    )

    emit:
    selected_positions = RESIDUE_SELECTION.out[0]
    candidates_csv = PYROSETTA_SCREENING.out[0]
    candidates_fasta = PYROSETTA_SCREENING.out[1]
    filtered_candidates_csv = DDG_FILTER.out[0]
    filtered_candidates_fasta = DDG_FILTER.out[1]
    screening_samplesheet = DDG_FILTER.out[2]
    selection_summary = RESIDUE_SELECTION.out[1]
    ranked_candidates = RANK_SCREENING.out
}
