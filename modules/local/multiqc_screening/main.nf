process MULTIQC_SCREENING {

    tag "multiqc_screening"

    container 'quay.io/phipsiprot/reporting:dev'

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path final_results
    path plots

    output:
    path "multiqc_report.html"

    script:
    """
    mkdir -p multiqc_input

    cp ${final_results} multiqc_input/final_screening_results.csv
    cp ${plots} multiqc_input/ || true

    python3 - <<'PY'
    import csv
    from pathlib import Path

    csv_path = Path("multiqc_input/final_screening_results.csv")
    rows = list(csv.DictReader(open(csv_path)))

    ddgs = []
    for r in rows:
        try:
            ddgs.append(float(r["ddg"]))
        except Exception:
            pass

    best = min(ddgs) if ddgs else "NA"
    worst = max(ddgs) if ddgs else "NA"
    mean = sum(ddgs) / len(ddgs) if ddgs else "NA"

    with open("multiqc_input/screening_summary_mqc.yaml", "w") as f:
        f.write("id: 'screening_summary'\\n")
        f.write("section_name: 'PhiPsiProt Screening Summary'\\n")
        f.write("description: 'Summary of PyRosetta mutational screening results.'\\n")
        f.write("plot_type: 'table'\\n")
        f.write("data:\\n")
        f.write("  Screening:\\n")
        f.write(f"    candidates: {len(rows)}\\n")
        f.write(f"    best_ddg: {best}\\n")
        f.write(f"    worst_ddg: {worst}\\n")
        f.write(f"    mean_ddg: {mean}\\n")
    PY

    cat > multiqc_config.yaml <<'EOF'
    custom_data:
      screening_summary:
        file_format: 'yaml'
    sp:
      screening_summary:
        fn: 'screening_summary_mqc.yaml'
    EOF

    multiqc multiqc_input \
        --config multiqc_config.yaml \
        --filename multiqc_report.html \
        --outdir .
    """
}
