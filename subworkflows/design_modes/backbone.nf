include { DESIGN_DUMMY } from '../../modules/local/design_dummy/main'

workflow BACKBONE {

    take:
    mode

    main:
    DESIGN_DUMMY(mode)
}
