include { DESIGN_DUMMY } from '../../modules/local/design_dummy/main'

workflow SATURATION {

    take:
    mode

    main:
    DESIGN_DUMMY(mode)
}
