/*
===============================================================================
RFDIFFUSION_PARTIAL

Genera variantes estructurales alrededor de un PDB existente utilizando
partial diffusion.

Ejemplo oficial reproducido:
- PDB: 2KL8.pdb
- longitud: 79 residuos
- diseños: 10
- partial_T: 10
===============================================================================
*/

process RFDIFFUSION_PARTIAL {

    tag "${input_pdb.simpleName}_partialT${partial_t}"

    executor 'local'

    publishDir "${params.outdir}/rfdiffusion_partial", mode: 'copy'

    input:
    path input_pdb
    val contig
    val num_designs
    val partial_t

    output:
    path "partial_*.pdb", emit: pdbs
    path "partial_*.trb", emit: trbs

    script:
    """
    mkdir -p output

    docker run --rm --gpus all \
        -v "\$PWD/${input_pdb}:/input/${input_pdb}:ro" \
        -v "\$PWD/output:/output" \
        rosettacommons/rfdiffusion:latest \
        inference.input_pdb=/input/${input_pdb} \
        inference.output_prefix=/output/partial \
        'contigmap.contigs=[${contig}]' \
        inference.num_designs=${num_designs} \
        diffuser.partial_T=${partial_t}

    cp output/partial_*.pdb .
    cp output/partial_*.trb .
    """
}
