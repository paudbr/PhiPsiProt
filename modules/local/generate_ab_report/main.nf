/*
 * GENERATE_AB_REPORT
 * HTML report con las tres etapas del pipeline de anticuerpos.
 * Llama a bin/generate_ab_report.py.
 */
process GENERATE_AB_REPORT {
    tag "antibody_report"
    label 'process_single'

    container 'docker.io/python:3.11-slim'

    input:
    path ranked_esm
    path ranked_af3
    path ranked_haddock
    path pdbs_esm
    path pdbs_af3

    output:
    path "PhiPsiProt_antibody_report.html", emit: report
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def haddock_arg = ranked_haddock.name != 'NO_FILE' ? "--ranked-haddock ${ranked_haddock}" : ""
    def esm_dir_arg = pdbs_esm.name      != 'NO_FILE' ? "--pdbs-esm ${pdbs_esm}"             : ""
    def af3_dir_arg = pdbs_af3.name      != 'NO_FILE' ? "--pdbs-af3 ${pdbs_af3}"             : ""
    """
    generate_ab_report.py \\
        --ranked-esm     ${ranked_esm} \\
        --ranked-af3     ${ranked_af3} \\
        ${haddock_arg} \\
        ${esm_dir_arg} \\
        ${af3_dir_arg} \\
        --out PhiPsiProt_antibody_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo "<html><body>stub report</body></html>" > PhiPsiProt_antibody_report.html
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
