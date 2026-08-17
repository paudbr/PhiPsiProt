process RFDIFFUSION_UNCONDITIONAL {

    tag "unconditional_150aa"

    executor 'local'

    publishDir "${params.outdir}/rfdiffusion_unconditional", mode: 'copy'

    output:
    path "test_*.pdb", emit: pdbs
    path "test_*.trb", emit: trbs

    script:
    """
    mkdir -p output

    docker run --rm --gpus all \
        -v \$PWD/output:/output \
        rosettacommons/rfdiffusion:latest \
        'contigmap.contigs=[150-150]' \
        inference.output_prefix=/output/test \
        inference.num_designs=10

    cp output/test_*.pdb .
    cp output/test_*.trb .
    """
}
