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
    filtered_candidates_csv = DDG_FILTER.out[0]
    filtered_candidates_fasta = DDG_FILTER.out[1]
    screening_samplesheet = DDG_FILTER.out[2]
}
