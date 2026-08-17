/*
===============================================================================
RESIDUE_SELECTION

Purpose:
Select mutable residues for mutational screening.

Modes handled directly by this module:
- manual: user-defined positions
- pocket: residues near a ligand
- interface: residues near another chain
- all: all residues in target chain

Note:
- auto_pocket is handled by AUTO_POCKET_DETECTION upstream.

Outputs:
- selected_positions.txt
- selection_summary.csv
===============================================================================
*/

process RESIDUE_SELECTION {

    tag "$input_pdb"

    container 'quay.io/phipsiprot_pyrosetta:dev'

    publishDir "${params.outdir}/residue_selection", mode: 'copy'

    input:
    path input_pdb
    val target_chain
    val positions
    val selection_mode
    val ligand_resname
    val interface_chain
    val distance_cutoff

    output:
    path "selected_positions.txt"
    path "selection_summary.csv"

    script:
    def mode = selection_mode ?: 'manual'
    def chain = target_chain ?: 'A'
    def cutoff = distance_cutoff ?: 6.0
    def pos = positions ?: 'ALL'
    def ligand = ligand_resname ?: ''
    def iface = interface_chain ?: ''

    def optional_args = ""
    if (mode == 'manual') {
        optional_args += " --positions '${pos}'"
    }
    if (mode == 'pocket') {
        optional_args += " --ligand_resname '${ligand_resname}'"
    }
    if (mode == 'interface') {
        optional_args += " --interface_chain '${interface_chain}'"
    }

    """
    python3 $projectDir/bin/select_mutable_residues.py \
        --input_pdb $input_pdb \
        --target_chain '${chain}' \
        --selection_mode '${mode}' \
        --distance_cutoff '${cutoff}' \
        ${optional_args} \
        --output selected_positions.txt

    selected_positions=\$(cat selected_positions.txt)

    python3 - << EOF
import csv

row = {
    "selection_mode": "${mode}",
    "target_chain": "${chain}",
    "positions": "${pos}",
    "ligand_resname": "${ligand}",
    "interface_chain": "${iface}",
    "distance_cutoff": "${cutoff}",
    "selected_positions": "\${selected_positions}",
}

with open("selection_summary.csv", "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=row.keys())
    writer.writeheader()
    writer.writerow(row)
EOF
    """
}