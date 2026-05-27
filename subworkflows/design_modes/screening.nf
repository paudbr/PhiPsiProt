include { PYROSETTA_SCREENING } from '../../modules/local/pyrosetta_screening/main'
include { DDG_FILTER } from '../../modules/local/ddg_filter/main'

workflow SCREENING {

    take:
    input_pdb
    target_chain
    positions

    main:
    PYROSETTA_SCREENING(input_pdb, target_chain, positions)
    DDG_FILTER(PYROSETTA_SCREENING.out[0])

    emit:
    candidates_csv = PYROSETTA_SCREENING.out[0]
    candidates_fasta = PYROSETTA_SCREENING.out[1]
    filtered_candidates = DDG_FILTER.out
}
