include { DESIGN_DUMMY } from '../../modules/local/design_dummy_ab/main'

workflow ANTIBODY_DESIGN {

    take:
    mode

    main:
    DESIGN_DUMMY_AB(mode)
}
