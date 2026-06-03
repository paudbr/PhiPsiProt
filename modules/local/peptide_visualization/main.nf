process PEPTIDE_VISUALIZATION {

    tag "peptide_visualization"

    executor 'local'

    publishDir "${params.outdir}/visualization", mode: 'copy'

    input:
    path target_pdb
    path binder_pdbs

    output:
    path "*.png"

    script:
    """
    for pdb in $binder_pdbs; do
        name=\$(basename \$pdb .pdb)

        cat > view_\${name}.pml <<EOF
load $target_pdb, target
load \$pdb, binder

hide everything
show cartoon, target
show cartoon, binder and chain A

color cyan, target
color green, binder and chain A

show sticks, target and resi 103
color red, target and resi 103

show spheres, target and resn ZN
color yellow, target and resn ZN

zoom target or binder, 10
png \${name}_overview.png, width=1600, height=1200, ray=1

zoom target and resi 103, 12
png \${name}_GLN103_zoom.png, width=1600, height=1200, ray=1

quit
EOF

        docker run --rm \
            -v \$PWD:/work \
            phipsiprot/pymol_visualization:dev \
            pymol -cq /work/view_\${name}.pml
    done

    sudo chown -R \$(id -u):\$(id -g) .
    """
}
