process DESIGN_DUMMY {
    tag "$mode"

    input:
    val mode

    output:
    path "candidates.csv"

    script:
    """
    echo "id,fasta,pdb,mode" > candidates.csv
    echo "dummy_${mode},dummy.fasta,dummy.pdb,${mode}" >> candidates.csv
    """
}
