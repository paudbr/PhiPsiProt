/*
 * Generate HADDOCK3 AIRs for antibody-antigen docking.
 * Reads CDR residues from HLT REMARK annotations automatically.
 * Hotspots in chain T format: "T305,T456"
 */
process GENERATE_AB_HADDOCK_AIRS {
    tag "${meta.id}"
    label 'process_single'

    container 'docker.io/python:3.11-slim'

    input:
    tuple val(meta),
          path(antibody_pdb),
          path(antigen_pdb)
    val  hotspot_residues   // "T305,T456"

    output:
    tuple val(meta),
          path(antibody_pdb),
          path(antigen_pdb),
          path("${meta.id}_airs.tbl"), emit: docking_input
    path "versions.yml",              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 <<'PYSCRIPT'
sample_id     = "${meta.id}"
hotspot_raw   = "${hotspot_residues}"

hotspot_resnums = [int(h[1:]) for h in hotspot_raw.split(",")
                   if h.strip().startswith("T")]

# Read CDR residues from HLT REMARK PDBinfo-LABEL annotations
cdr_res = []
with open("${antibody_pdb}") as fh:
    for line in fh:
        if "PDBinfo-LABEL" in line:
            parts = line.split()
            try:
                cdr_res.append(int(parts[-2]))
            except (ValueError, IndexError):
                pass

lines = [
    "! HADDOCK3 AIRs — de novo antibody vs antigen",
    f"! Epitope (chain T → segid B): {hotspot_raw}",
    f"! CDR residues (chain H/L → segid A): {cdr_res}", "!"
]
if cdr_res and hotspot_resnums:
    cdr_sel = " or resid ".join(str(r) for r in cdr_res)
    hot_sel = " or resid ".join(str(r) for r in hotspot_resnums)
    for r in hotspot_resnums:
        lines.append(f"assign ( resid {r} and segid B ) ( ( resid {cdr_sel} ) and segid A ) 2.0 2.0 0.0")
    for r in cdr_res:
        lines.append(f"assign ( resid {r} and segid A ) ( ( resid {hot_sel} ) and segid B ) 2.0 2.0 0.0")

with open(f"{sample_id}_airs.tbl", "w") as fh:
    fh.write("\\n".join(lines) + "\\n")
print(f"Written {len(lines)-4} AIRs for {sample_id}")
PYSCRIPT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo "! stub AIRs" > ${meta.id}_airs.tbl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
