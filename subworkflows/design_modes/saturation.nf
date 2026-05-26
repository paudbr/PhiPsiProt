include { PYROSETTA_SATURATION } from '../../modules/local/pyrosetta_saturation/main'

workflow SATURATION {

    take:
    input_pdb

    main:
    PYROSETTA_SATURATION(input_pdb)
}
