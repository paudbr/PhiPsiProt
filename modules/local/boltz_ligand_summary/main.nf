process BOLTZ_LIGAND_SUMMARY {
    tag "boltz_ligand_summary"
    label 'process_single'

    input:
    path confidence_jsons, stageAs: "confidence_jsons/*"

    output:
    path "boltz2_ligand_scores.tsv"  , emit: scores
    path "boltz2_ligand_summary.tsv" , emit: summary
    path "versions.yml"              , emit: versions

    script:
    """
    boltz_ligand_summary.py confidence_jsons/*.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
