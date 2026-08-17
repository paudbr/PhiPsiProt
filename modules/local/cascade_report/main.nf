/*
 * CASCADE_REPORT — report HTML de la cascada (dos pestañas ESM/Chai, NGL por
 * pLDDT, tabla de métricas, gráficos score y ΔG-vs-TM).
 *
 * Stagea las entradas en las subcarpetas que espera generate_cascade_report.py:
 *   metrics_esm/  metrics_chai/  pdbs_esm/  pdbs_chai/  + los dos ranked.csv
 */
process CASCADE_REPORT {
    label 'process_single'
    container "nf-core/proteinfold_esmfold:2.0.0"   // trae python+biopython si hiciera falta
    publishDir "${params.outdir}/cascade_report", mode: 'copy'

    input:
    path(ranked_esm,  stageAs: "ranked_esm.csv")   // ranked.csv ronda ESM
    path(ranked_chai, stageAs: "ranked_chai.csv")  // ranked.csv ronda Chai
    path(esm_pdbs,  stageAs: "pdbs_esm/*")         // PDBs de ESMFold
    path(chai_pdbs, stageAs: "pdbs_chai/*")        // PDBs de Chai (convertidos a pdb)

    output:
    path "PhiPsiProt_cascade_report.html", emit: report

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    generate_cascade_report.py
    """
}