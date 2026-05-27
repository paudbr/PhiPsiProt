include { BACKBONE_DESIGN } from '../../modules/local/backbone_design/main'

workflow BACKBONE {

    take:
    input_pdb

    main:
    BACKBONE_DESIGN(input_pdb)

    emit:
    designed_pdbs = BACKBONE_DESIGN.out[0]
    design_metadata = BACKBONE_DESIGN.out[1]
}