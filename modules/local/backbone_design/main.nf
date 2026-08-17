/*
===============================================================================
BACKBONE_DESIGN

Purpose:
Generate candidate backbone structures from an input PDB using RFdiffusion.

Inputs:
- input_pdb

Outputs:
- backbone_design_*.pdb
- backbone_design_*.trb

Notes:
- Runs RFdiffusion inside the official Docker image.
- Can run in partial diffusion mode when params.partial_T > 0.
- If params.partial_T is null or 0, partial diffusion is not applied.
- The contig parameter controls which regions are redesigned or kept fixed.
- PDB files are used downstream by ProteinMPNN.
- TRB files contain RFdiffusion metadata.
===============================================================================
*/

process BACKBONE_DESIGN {

    tag "$input_pdb"

    executor 'local'

    publishDir "${params.outdir}/designs/backbones", mode: 'copy'

    input:
    path input_pdb

    output:
    path "backbone_design_*.pdb"
    path "backbone_design_*.trb"

    script:
    def partialArg = ''
    if (
        params.partial_T != null &&
        params.partial_T.toString() != 'null' &&
        params.partial_T.toString().isInteger() &&
        params.partial_T.toInteger() > 0
    ) {
        partialArg = "diffuser.partial_T=${params.partial_T}"
    }

    """
    mkdir -p rfdiffusion_output

    INPUT_PDB=\$(readlink -f $input_pdb)

    python3 $projectDir/bin/build_backbone_contig.py \
        --pdb $input_pdb \
        --chain ${params.chain ?: "A"} \
        --design_strategy ${params.design_strategy ?: "local_redesign"} \
        --redesign_region "${params.redesign_region ?: "121-126"}" \
        --output contig.txt

    CONTIG=\$(cat contig.txt)

    echo "Using contig: \$CONTIG"

    docker run --rm --gpus all \
        -v \$(dirname \$INPUT_PDB):/input \
        -v \$PWD/rfdiffusion_output:/output \
        docker.io/rosettacommons/rfdiffusion:latest \
        inference.output_prefix=/output/backbone_design \
        inference.input_pdb=/input/\$(basename \$INPUT_PDB) \
        inference.num_designs=${params.num_designs ?: 1} \
        'contigmap.contigs=['"\$CONTIG"']' \
        diffuser.T=${params.diffusion_T ?: 50} \
        ${partialArg} \
        denoiser.noise_scale_ca=${params.noise_scale_ca ?: 1.0} \
        denoiser.noise_scale_frame=${params.noise_scale_frame ?: 1.0} \
        ${params.final_step ? "inference.final_step=${params.final_step}" : ""}

    sudo chown -R \$(id -u):\$(id -g) rfdiffusion_output

    cp rfdiffusion_output/backbone_design_*.pdb .
    cp rfdiffusion_output/backbone_design_*.trb .
    """
}