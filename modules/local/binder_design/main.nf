process BINDER_DESIGN {

    tag "$target_pdb"

    executor 'local'

    publishDir "${params.outdir}/designs/binders", mode: 'copy'

    input:
    path target_pdb

    output:
    path "binder_design_*.pdb"
    path "binder_design_*.trb"

    script:
    """
    mkdir -p rfdiffusion_output

    INPUT_PDB=\$(readlink -f $target_pdb)

    docker run --rm --gpus all \
        -v \$(dirname \$INPUT_PDB):/input \
        -v \$PWD/rfdiffusion_output:/output \
        docker.io/rosettacommons/rfdiffusion:latest \
        inference.output_prefix=/output/binder_design \
        inference.input_pdb=/input/\$(basename \$INPUT_PDB) \
        inference.num_designs=${params.num_designs ?: 1} \
        'contigmap.contigs=[${params.target_contig ?: "A1-100"}/0 ${params.binder_length ?: 40}-${params.binder_length ?: 40}]' \
        'ppi.hotspot_res=[${params.hotspot_res}]' \
        diffuser.T=${params.diffusion_T ?: 50} \
        denoiser.noise_scale_ca=${params.noise_scale_ca ?: 0.5} \
        denoiser.noise_scale_frame=${params.noise_scale_frame ?: 0.5}

    sudo chown -R \$(id -u):\$(id -g) rfdiffusion_output

    cp rfdiffusion_output/binder_design_*.pdb .
    cp rfdiffusion_output/binder_design_*.trb .
    """
}
