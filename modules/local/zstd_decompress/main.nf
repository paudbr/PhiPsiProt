process ZSTD_DECOMPRESS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0a/0a27033ae5d8add5059f44c62a6004bfcd061d33020edee095fbb204e6f32fee/data' :
        'community.wave.seqera.io/library/zstd:b5faa75d5b75be7f' }"

    input:
    tuple val(meta), path(archive, stageAs: 'input.zst')

    output:
    tuple val(meta), path("*.decompressed"), emit: decompressed
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    # nombre base seguro (NO depende de output scope)
    OUTNAME="${meta.id ?: "output"}"

    # asegurar archivo real (Wave-safe)
    cp -L input.zst input_real.zst

    # descompresión
    zstd -d input_real.zst

    # detectar archivo resultante sin ls/grep (robusto)
    file_out=\$(find . -maxdepth 1 -type f ! -name "*.zst" | head -n 1)

    # renombrar salida de forma estable
    mv "\$file_out" "\${OUTNAME}.decompressed"

    # versión
    cat <<EOF > versions.yml
${task.process}:
    zstd: \$(zstd --version | grep -o 'v[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+')
EOF
    """
}