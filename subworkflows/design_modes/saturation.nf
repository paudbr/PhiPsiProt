include { PYROSETTA_SATURATION } from '../../modules/local/pyrosetta_saturation/main'

workflow SATURATION {

    take:
    input_pdb
    target_chain
    positions

    main:
    PYROSETTA_SATURATION(input_pdb, target_chain, positions)

    emit:
    candidates = PYROSETTA_SATURATION.out
}
