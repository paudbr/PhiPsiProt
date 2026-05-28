include { BACKBONE_DESIGN } from '../../modules/local/backbone_design/main'
include { BACKBONE_SUMMARY } from '../../modules/local/backbone_summary/main'
include { PROTEINMPNN_DESIGN } from '../../modules/local/proteinmpnn_design/main'
include { MERGE_BACKBONE_RESULTS } from '../../modules/local/merge_backbone_results/main'

workflow BACKBONE {

    take:
    input_pdb

    main:

    BACKBONE_DESIGN(input_pdb)

    BACKBONE_SUMMARY(BACKBONE_DESIGN.out[0])

    PROTEINMPNN_DESIGN(BACKBONE_DESIGN.out[0])

    MERGE_BACKBONE_RESULTS(
        BACKBONE_SUMMARY.out[0],
        PROTEINMPNN_DESIGN.out
    )

    emit:
    final_backbone_results = MERGE_BACKBONE_RESULTS.out
}
