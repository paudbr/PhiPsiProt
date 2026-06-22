process CANDIDATES_SMILES_TO_BOLTZ {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.11 conda-forge::pyyaml=6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    tuple val(meta), val(smiles)
    val   receptor_sequence
    val   pocket_residues

    output:
    tuple val(meta), path("*.yaml"), emit: yaml
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def ligand_id  = "L"
    def protein_id = "A"
    // si pocket_residues viene vacio/null, no se pasa el flag
    def pocket_arg = pocket_residues ? "--pocket-residues \"${pocket_residues}\"" : ""
    def chains_arg = params.protein_chains ? "--protein-chains \"${params.protein_chains}\"" : ""
    """
    candidates_smiles_to_boltz.py \\
        --candidate-id "${meta.id}" \\
        --smiles "${smiles}" \\
        --receptor-sequence "${receptor_sequence}" \\
        --protein-id "${protein_id}" \\
        --ligand-id "${ligand_id}" \\
        --output "${prefix}.yaml" \\
        ${chains_arg} \\
        ${pocket_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pyyaml: \$(python -c "import yaml; print(yaml.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.yaml
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pyyaml: stub
END_VERSIONS
    """
}