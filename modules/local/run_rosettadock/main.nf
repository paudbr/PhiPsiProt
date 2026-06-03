process RUN_ROSETTADOCK {
    when:
    params.docking_tool == 'rosettadock'


    tag "$id"

        input:
        tuple val(id), path(receptor), path(ligand)

        output:
        tuple val(id), path("*.score"), emit: scores

        script:
        """
        rosetta_scripts \
            -s ${receptor} \
            -l ${ligand} \
            -out:file:scorefile ${id}.score
    """
}