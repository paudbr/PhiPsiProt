process COLLECT_METRICS {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(pdb_or_cif), val(tool), path(scores_json)

    output:
    tuple val(meta), path("${meta.id}_${tool}_plddt.tsv"),  emit: plddt
    tuple val(meta), path("${meta.id}_${tool}_scores.tsv"), emit: scores

    script:
    """
    python3 - << 'EOF'
import json, os

TAB = chr(9)
NL  = chr(10)

tool   = "${tool}"
prefix = "${meta.id}_${tool}"
pdb    = "${pdb_or_cif}"
scores = "${scores_json}"

plddt_vals = []
ext = os.path.splitext(pdb)[1].lower()

if ext == ".pdb":
    with open(pdb) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    chain = line[21]
                    res   = int(line[22:26])
                    bfac  = float(line[60:66])
                    plddt_vals.append((chain, res, bfac))
                except ValueError:
                    pass

elif ext == ".cif":
    with open(pdb) as f:
        for line in f:
            parts = line.split()
            if len(parts) > 14 and parts[0] in ("ATOM", "HETATM"):
                try:
                    plddt_vals.append((parts[6], int(parts[8]), float(parts[14])))
                except (ValueError, IndexError):
                    pass

seen = {}
for chain, res, val in plddt_vals:
    key = (chain, res)
    if key not in seen:
        seen[key] = val

with open(f"{prefix}_plddt.tsv", "w") as fh:
    fh.write("chain" + TAB + "residue" + TAB + "plddt" + NL)
    for (chain, res), val in sorted(seen.items()):
        fh.write(chain + TAB + str(res) + TAB + f"{val:.4f}" + NL)

ptm, iptm = None, None
if os.path.isfile(scores) and os.path.getsize(scores) > 0:
    with open(scores) as f:
        d = json.load(f)
    ptm  = d.get("ptm")  or d.get("iptm_ptm", {}).get("ptm")
    iptm = d.get("iptm") or d.get("iptm_ptm", {}).get("iptm")

with open(f"{prefix}_scores.tsv", "w") as fh:
    fh.write("metric" + TAB + "value" + NL)
    if ptm  is not None: fh.write("ptm"  + TAB + f"{ptm:.4f}"  + NL)
    if iptm is not None: fh.write("iptm" + TAB + f"{iptm:.4f}" + NL)

EOF
    """
}
