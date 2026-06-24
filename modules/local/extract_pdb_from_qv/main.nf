// Módulo nuevo: EXTRACT_PDB_FROM_QV
// usa rfantibody_local:1.0.0 que SÍ tiene qvextractspecific
process EXTRACT_PDB_FROM_QV {
    tag "${meta.id}"
    container 'quay.io/rfantibody_local:1.0.0'
    containerOptions = '-v /mnt/alphafold3_db:/mnt/alphafold3_db --entrypoint ""'

    input:
    tuple val(meta), path(quiver), val(design_id)

    output:
    tuple val(meta), path("${design_id}.pdb"), emit: pdb

    script:
    """
    qvextractspecific ${quiver} ${design_id} -o ./
    """

    stub:
    """
    touch ${design_id}.pdb
    """
}