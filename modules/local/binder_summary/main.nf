process BINDER_SUMMARY {

    tag "binder_summary"

    publishDir "${params.outdir}/intermediate", mode: 'copy'

    input:
    path binder_pdbs

    output:
    path "binder_candidates.csv"

    script:
    """
    echo "candidate_id,mode,method,pdb,binder_length,status" > binder_candidates.csv

    for pdb in ${binder_pdbs}; do
        id=\$(basename "\$pdb" .pdb)
        length=\$(grep '^ATOM' "\$pdb" | awk '{print \$5,\$6}' | sort -u | wc -l)

        echo "\$id,binder,rfdiffusion,\$(basename \$pdb),\$length,ready_for_sequence_design" >> binder_candidates.csv
    done
    """
}
