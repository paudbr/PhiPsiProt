process PYMOL_VISUALIZATION {

    tag "$input_pdb"

    executor 'local'

    publishDir "${params.outdir}/visualization", mode: 'copy'

    input:
    path input_pdb
    path selected_positions
    val target_chain

    output:
    path "selected_residues.pml"
    path "selected_residues.png"
    path "selected_residues.pse"
    path "selected_residues.html"

    script:
    """
    python $projectDir/bin/generate_pymol_report.py \
        --input_pdb $input_pdb \
        --positions $selected_positions \
        --target_chain $target_chain \
        --output_prefix selected_residues

    docker run --rm \
        -v \$PWD:/data \
        quay.io/phipsiprot_pymol_visualization:dev \
        xvfb-run -a pymol -cq /data/selected_residues.pml

    sudo chown -R \$(id -u):\$(id -g) selected_residues.*
    """
}
