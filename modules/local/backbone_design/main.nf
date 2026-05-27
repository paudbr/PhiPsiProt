process BACKBONE_DESIGN {

    tag "$input_pdb"

    executor 'local'

    publishDir "${params.outdir}/backbone", mode: 'copy'

    input:
    path input_pdb

    output:
    path "backbone_design_0.pdb"
    path "backbone_design_0.trb"

    script:
    """
    mkdir -p rfdiffusion_output

    INPUT_PDB=\$(readlink -f $input_pdb)

    docker run --rm --gpus all \
        -v \$(dirname \$INPUT_PDB):/input \
        -v \$PWD/rfdiffusion_output:/output \
        docker.io/rosettacommons/rfdiffusion:latest \
        inference.output_prefix=/output/backbone_design \
        inference.input_pdb=/input/\$(basename \$INPUT_PDB) \
        inference.num_designs=${params.num_designs ?: 1} \
        'contigmap.contigs=[${params.contig ?: "4-4"}]' \
        diffuser.partial_T=${params.partial_T ?: 1}

    cp rfdiffusion_output/backbone_design_0.pdb .
    cp rfdiffusion_output/backbone_design_0.trb .
    """
}
