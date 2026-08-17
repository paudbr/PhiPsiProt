/*
 * REFERENCE_METRICS — calcula ΔG y SAP del PDB original UNA sola vez
 * (relax restringido para el ΔG), para usarlos como referencia de los
 * deltas (rosetta_ddg y delta_sap) en todos los candidatos.
 */
process REFERENCE_METRICS {
    tag "reference"
    label 'process_medium'
    container "phipsiprot/pyrosetta:dev-metrics"

    input:
    path reference_pdb
    val  relax_mode

    output:
    path "reference_ref.tsv", emit: ref
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    structural_metrics.py \\
        --id reference \\
        --model ${reference_pdb} \\
        --relax ${relax_mode} \\
        --dg-only \\
        --out reference_ref.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyrosetta: \$(python3 -c "import pyrosetta; print(getattr(pyrosetta,'__version__','unknown'))" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}

/*
 * STRUCTURAL_METRICS — métricas estructurales sobre el PDB de un modelo.
 * Emite TSV largo: id  metric  value. Columnas:
 *   rosetta_dg, rosetta_ddg, tm_score, mean_plddt, n_cysteines,
 *   net_charge, net_charge_abs, sap_score, delta_sap
 *
 * ref_dg / ref_sap: valores del original (de REFERENCE_METRICS) para los deltas.
 */
process STRUCTURAL_METRICS {
    tag "$meta.id"
    label 'process_medium'
    container "phipsiprot/pyrosetta:dev-metrics"

    input:
    tuple val(meta), path(model_pdb)
    path  reference_pdb
    val   ph
    val   relax_mode
    val   ref_dg            // ΔG del original; '' o 'NA' si no aplica
    val   ref_sap           // SAP del original; '' o 'NA' si no aplica

    output:
    tuple val(meta), path("${meta.id}_metrics.tsv"), emit: metrics
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ref_arg    = (reference_pdb.name in ["NO_FILE", "NO_FILE_PAE"]) ? "" : "--ref ${reference_pdb}"
    def refdg_arg  = (ref_dg  == null || ref_dg.toString()  in ["", "NA"]) ? "" : "--ref-dg ${ref_dg}"
    def refsap_arg = (ref_sap == null || ref_sap.toString() in ["", "NA"]) ? "" : "--ref-sap ${ref_sap}"
    """
    structural_metrics.py \\
        --id ${meta.id} \\
        --model ${model_pdb} \\
        ${ref_arg} \\
        ${refdg_arg} \\
        ${refsap_arg} \\
        --ph ${ph} \\
        --relax ${relax_mode} \\
        --out ${meta.id}_metrics.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyrosetta: \$(python3 -c "import pyrosetta; print(getattr(pyrosetta,'__version__','unknown'))" 2>/dev/null || echo "unknown")
        propka: \$(python3 -c "import propka; print(getattr(propka,'__version__','unknown'))" 2>/dev/null || echo "unknown")
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)" 2>/dev/null || echo "unknown")
        freesasa: \$(python3 -c "import freesasa; print(getattr(freesasa,'__version__','unknown'))" 2>/dev/null || echo "unknown")
        tmtools: \$(python3 -c "import tmtools; print(getattr(tmtools,'__version__','unknown'))" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
