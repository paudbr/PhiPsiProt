/*
===============================================================================
Test workflow:
RFdiffusion unconditional backbone generation followed by ProteinMPNN
sequence design.
===============================================================================
*/

include { RFDIFFUSION_UNCONDITIONAL } from '../../modules/local/rfdiffusion_unconditional/main'
include { PROTEINMPNN_UNCONDITIONAL } from '../../modules/local/proteinmpnn_unconditional/main'

workflow {

    RFDIFFUSION_UNCONDITIONAL()

    designed_pdb_ch = RFDIFFUSION_UNCONDITIONAL.out.pdbs.flatten()

    PROTEINMPNN_UNCONDITIONAL(designed_pdb_ch)
}
