include { RFANTIBODY_DESIGN } from '../../modules/local/rfantibody_design/main'
include { RFANTIBODY_PROTEINMPNN } from '../../modules/local/rfantibody_proteinmpnn/main'
include { RFANTIBODY_RF2 } from '../../modules/local/rfantibody_rf2/main'
include { RFANTIBODY_SCORES } from '../../modules/local/rfantibody_scores/main'
include { FINALIZE_ANTIBODY_RESULTS } from '../../modules/local/finalize_antibody_results/main'

workflow ANTIBODY_DESIGN {

    take:
    target_pdb
    framework_pdb

    main:

    RFANTIBODY_DESIGN(target_pdb, framework_pdb)

    RFANTIBODY_PROTEINMPNN(RFANTIBODY_DESIGN.out)

    RFANTIBODY_RF2(RFANTIBODY_PROTEINMPNN.out)

    RFANTIBODY_SCORES(RFANTIBODY_RF2.out)

    FINALIZE_ANTIBODY_RESULTS(RFANTIBODY_SCORES.out)

    emit:
    final_antibody_candidates = FINALIZE_ANTIBODY_RESULTS.out
}
