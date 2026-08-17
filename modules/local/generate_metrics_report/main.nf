process GENERATE_METRICS_REPORT {
    label 'process_single'
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    path(plddt_files, stageAs: "plddt/*")
    path(scores_files, stageAs: "scores/*")
    path(pdb_files,    stageAs: "pdbs/*")

    output:
    path "PhiPsiProt_metrics_report.html", emit: report

    script:
    """
    python3 ${moduleDir}/generate_report.py
    """
}
