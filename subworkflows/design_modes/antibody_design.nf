include { DESIGN_DUMMY } from '../../modules/local/design_dummy/main'

workflow ANTIBODY_DESIGN {

    take:
    mode

    main:
    DESIGN_DUMMY(mode)
}
