/*
===============================================================================
TEST 2: RFdiffusion partial diffusion

Reproduce el ejemplo oficial design_partialdiffusion.sh.
===============================================================================
*/

include { RFDIFFUSION_PARTIAL } from '../../modules/local/rfdiffusion_partial/main'

workflow {

    input_pdb_ch = Channel.fromPath(
        params.input_pdb,
        checkIfExists: true
    )

    RFDIFFUSION_PARTIAL(
        input_pdb_ch,
        params.rf_contig,
        params.rf_num_designs,
        params.partial_T
    )
}
