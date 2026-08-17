/*
===============================================================================
UNCONDITIONAL DESIGN WORKFLOW

RFdiffusion unconditional backbone generation
optionally followed by ProteinMPNN sequence design.
===============================================================================
*/

include {
    RFDIFFUSION_UNCONDITIONAL
} from '../../modules/local/rfdiffusion_unconditional/main'

include {
    PROTEINMPNN_UNCONDITIONAL
} from '../../modules/local/proteinmpnn_unconditional/main'

workflow UNCONDITIONAL_DESIGN {

    take:
    run_mpnn

    main:

    RFDIFFUSION_UNCONDITIONAL()

    if (run_mpnn) {

        designed_pdb_ch = RFDIFFUSION_UNCONDITIONAL.out.pdbs.flatten()

        PROTEINMPNN_UNCONDITIONAL(designed_pdb_ch)
    }

    emit:
    pdbs = RFDIFFUSION_UNCONDITIONAL.out.pdbs
}
