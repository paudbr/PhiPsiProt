process CANDIDATES_SMILES_TO_AF3_FASTA {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c7/c7dabd3f132a613fb11ee27c66e9517eb7649eee64f4e4f63747841105883b40/data' :
        'community.wave.seqera.io/library/biopython_python:06582b7b722f3db3' }"

    input:
    tuple val(meta), val(smiles)
    val   receptor_sequence

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // cadenas de proteina: "A,B" -> homodimero. Si no, una sola cadena.
    def chain_ids = (params.protein_chains ?: 'A')
                        .split(',')
                        .collect { it.trim() }
                        .findAll { it }
    // construir las entradas de proteina (misma secuencia, una por cadena)
    def protein_blocks = chain_ids
                            .collect { cid -> ">protein|${cid}\n${receptor_sequence}" }
                            .join('\n')
    """
    cat > ${prefix}.fasta <<'FASTA'
${protein_blocks}
>smiles|${meta.id}
${smiles}
FASTA

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | sed 's/.*version //; s/ .*//')
    END_VERSIONS
    """
}