/*
 * COLLECT_AND_RANK_STRUCTURAL — junta los *_metrics.tsv de todos los candidatos,
 * aplica el scoring ponderado del YAML y produce ranked.csv.
 *
 * Se usa DOS veces en la cascada:
 *   - tras ESMFold: apply_gates=true,  select_top_n=true   (filtra y recorta)
 *   - tras Chai-1:  apply_gates=false, select_top_n=false  (solo ordena)
 *
 * NO necesita PyRosetta: solo Python + PyYAML. Por eso usa una imagen python
 * ligera de biocontainers (igual que COLLECT_AND_RANK del modo ligando),
 * no la imagen pesada de pyrosetta.
 */
process COLLECT_AND_RANK_STRUCTURAL {
    tag "${stage}"
    label 'process_single'
    conda "conda-forge::python=3.11 conda-forge::pyyaml=6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"
    publishDir "${params.outdir}/ranking/${stage}", mode: 'copy'

    input:
    path  metrics_tsvs            // muchos *_metrics.tsv (collect)
    path  scoring_config
    val   apply_gates             // true/false
    val   select_top_n            // true/false
    val   stage                   // 'esm' | 'chai'

    output:
    path "ranked.csv"  , emit: ranked
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def gates_flag = apply_gates   ? "--apply-gates"   : "--no-apply-gates"
    def topn_flag  = select_top_n  ? "--select-top-n"  : "--no-select-top-n"
    """
    rank_structural.py \\
        --config ${scoring_config} \\
        --metrics-glob '*_metrics.tsv' \\
        ${gates_flag} \\
        ${topn_flag} \\
        --out ranked.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pyyaml: \$(python -c "import yaml; print(yaml.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
