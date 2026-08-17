/*
===============================================================================
PARTIAL DIFFUSION WORKFLOW

Partial diffusion from an existing structure,
optionally followed by ProteinMPNN.
===============================================================================
*/

include {
    RFDIFFUSION_PARTIAL
} from '../../modules/local/rfdiffusion_partial/main'

include {
    PROTEINMPNN_UNCONDITIONAL
} from '../../modules/local/proteinmpnn_unconditional/main'

workflow PARTIAL_DESIGN {

    take:
    input_pdb_ch
    contig
    num_designs
    partial_t
    run_mpnn

    main:

    RFDIFFUSION_PARTIAL(
        input_pdb_ch,
        contig,
        num_designs,
        partial_t
    )

    if (run_mpnn) {

        partial_pdb_ch = RFDIFFUSION_PARTIAL.out.pdbs.flatten()

        PROTEINMPNN_UNCONDITIONAL(partial_pdb_ch)
    }

    emit:
    pdbs = RFDIFFUSION_PARTIAL.out.pdbs
}
