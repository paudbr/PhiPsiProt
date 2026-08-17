/*
===============================================================================
TEST:
RFdiffusion partial diffusion followed by ProteinMPNN sequence design.

Workflow:
1. Read an existing input PDB.
2. Generate structural variants using partial diffusion.
3. Design complete sequences for every generated backbone using ProteinMPNN.
===============================================================================
*/

include { RFDIFFUSION_PARTIAL } from '../../modules/local/rfdiffusion_partial/main'

include {
    PROTEINMPNN_UNCONDITIONAL
} from '../../modules/local/proteinmpnn_unconditional/main'

workflow {

    /*
    ---------------------------------------------------------------------------
    Input structure
    ---------------------------------------------------------------------------
    */

    input_pdb_ch = Channel.fromPath(
        params.input_pdb,
        checkIfExists: true
    )

    /*
    ---------------------------------------------------------------------------
    RFdiffusion partial diffusion
    ---------------------------------------------------------------------------
    */

    RFDIFFUSION_PARTIAL(
        input_pdb_ch,
        params.rf_contig,
        params.rf_num_designs,
        params.partial_T
    )

    /*
    Each generated PDB becomes an independent ProteinMPNN input.
    */

    partial_pdb_ch = RFDIFFUSION_PARTIAL.out.pdbs.flatten()

    /*
    ---------------------------------------------------------------------------
    ProteinMPNN complete sequence design
    ---------------------------------------------------------------------------
    */

    PROTEINMPNN_UNCONDITIONAL(partial_pdb_ch)
}
