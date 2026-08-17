/*
 * Run Chai-1
 */
process RUN_CHAI1 {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'
    container 'quay.io/nf-core/proteinfold_chai1:1.3.0'
  
    publishDir "${params.outdir}/run/chai1", mode: 'copy'

    shell '/bin/bash', '-euo', 'pipefail'

    input:
    tuple val(meta), path(fasta)
    path  weights_dir

    output:
    tuple val(meta), path("raw/**"), emit: raw
    tuple val(meta), path("${meta.id}_chai1.cif"), emit: top_ranked_cif
    tuple val(meta), path("raw/*ranked_*.cif")   , emit: cif
    tuple val(meta), path("${meta.id}_plddt.tsv"), emit: multiqc
    tuple val(meta), path("${meta.id}_ptm.tsv")  , emit: ptms
    tuple val(meta), path("${meta.id}_iptm.tsv") , optional: true, emit: iptms
    path  "versions.yml"                          , emit: versions

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    set -euo pipefail

    export CHAI_DOWNLOADS_DIR="${weights_dir}"
    export PATH="/chai1_venv/bin:\$PATH"

    awk '
    /^>/ {
        header = substr(\$0, 2)
        n = split(header, parts, "|")
        entity = parts[1]
        if (entity == "protein" || entity == "dna" || entity == "rna" ||
            entity == "ligand"  || entity == "glycan") {
            print \$0
        } else if (n >= 3) {
            split(parts[3], tmp, " ")
            print ">protein|name=" tmp[1]
        } else if (n == 2) {
            split(parts[2], tmp, " ")
            print ">protein|name=" tmp[1]
        } else {
            gsub(/ /, "_", header)
            print ">protein|name=" header
        }
        next
    }
    { print }
    ' ${fasta} > chai_ready.fasta

    /chai1_venv/bin/chai-lab fold \\
        chai_ready.fasta \\
        ./chai1_output \\
        ${args}

    mkdir -p raw

    cat > process_scores.py << 'PYEOF'
import sys, numpy as np, glob, os, shutil

prefix     = sys.argv[1]
output_dir = sys.argv[2]

score_files = sorted(glob.glob(os.path.join(output_dir, "scores.model_idx_*.npz")))
if not score_files:
    raise FileNotFoundError("No Chai-1 score files found in: " + output_dir)

# Rankings por aggregate_score
records = []
for sf in score_files:
    idx  = int(sf.split("model_idx_")[1].replace(".npz", ""))
    data = np.load(sf)
    agg  = float(data["aggregate_score"].flat[0])
    records.append((idx, agg))

records.sort(key=lambda x: x[1], reverse=True)

for rank, (idx, _score) in enumerate(records):
    src = os.path.join(output_dir, f"pred.model_idx_{idx}.cif")
    shutil.copy(src, f"raw/model_idx_{idx}_ranked_{rank}.cif")
    if rank == 0:
        shutil.copy(src, f"{prefix}_chai1.cif")

# pLDDT desde el CIF del modelo top
# B_iso_or_equiv contiene el pLDDT en los CIFs de Chai-1
top_idx = records[0][0]
top_cif = os.path.join(output_dir, f"pred.model_idx_{top_idx}.cif")

def parse_plddt_from_cif(cif_path):
    # Extrae pLDDT por residuo del bloque _atom_site del mmCIF.
    # Devuelve lista de (chain_id, seq_id, plddt) - un valor por residuo (CA o primer atomo).
    col_names = []
    in_loop   = False
    in_atom   = False
    col_b     = None
    col_chain = None
    col_seq   = None
    col_atom  = None
    seen      = {}

    with open(cif_path) as f:
        for line in f:
            line = line.rstrip()

            if line.startswith("loop_"):
                col_names = []
                in_loop   = True
                in_atom   = False
                continue

            if in_loop and line.startswith("_atom_site."):
                col_names.append(line.strip())
                if line.strip() == "_atom_site.B_iso_or_equiv":
                    col_b = len(col_names) - 1
                if line.strip() == "_atom_site.auth_asym_id":
                    col_chain = len(col_names) - 1
                if line.strip() == "_atom_site.auth_seq_id":
                    col_seq = len(col_names) - 1
                if line.strip() == "_atom_site.label_atom_id":
                    col_atom = len(col_names) - 1
                if col_b is not None:
                    in_atom = True
                continue

            if in_atom and (line.startswith("ATOM") or line.startswith("HETATM")):
                parts = line.split()
                indices = [i for i in [col_b, col_chain, col_seq, col_atom] if i is not None]
                if not indices or len(parts) <= max(indices):
                    continue
                atom_name = parts[col_atom]  if col_atom  is not None else "CA"
                chain     = parts[col_chain] if col_chain is not None else "A"
                seq_id    = parts[col_seq]   if col_seq   is not None else "0"
                b_val     = parts[col_b]     if col_b     is not None else "0"
                key = (chain, seq_id)
                # Queda con CA preferentemente, o el primer atomo del residuo
                if key not in seen or atom_name == "CA":
                    try:
                        seen[key] = float(b_val)
                    except ValueError:
                        pass

            elif in_atom and line.startswith("#"):
                in_atom = False

    return [(chain, int(seq), val) for (chain, seq), val in sorted(
        seen.items(), key=lambda x: (x[0][0], int(x[0][1]))
    )]

plddt_data = parse_plddt_from_cif(top_cif)

with open(f"{prefix}_plddt.tsv", "w") as fh:
    fh.write("chain\\tresidue\\tplddt\\n")
    for chain, res, val in plddt_data:
        fh.write(f"{chain}\\t{res}\\t{val:.4f}\\n")

# pTM
with open(f"{prefix}_ptm.tsv", "w") as fh:
    fh.write("model\\tptm\\taggregate_score\\n")
    for rank, (idx, agg) in enumerate(records):
        d   = np.load(os.path.join(output_dir, f"scores.model_idx_{idx}.npz"))
        ptm = float(d["ptm"].flat[0])
        fh.write(f"ranked_{rank}\\t{ptm:.4f}\\t{agg:.4f}\\n")

# ipTM
try:
    with open(f"{prefix}_iptm.tsv", "w") as fh:
        fh.write("model\\tiptm\\n")
        for rank, (idx, _) in enumerate(records):
            d    = np.load(os.path.join(output_dir, f"scores.model_idx_{idx}.npz"))
            iptm = float(d["iptm"].flat[0])
            fh.write(f"ranked_{rank}\\t{iptm:.4f}\\n")
except Exception as e:
    print(f"[WARN] ipTM not available: {e}", file=sys.stderr)

for f in os.listdir(output_dir):
    shutil.move(os.path.join(output_dir, f), os.path.join("raw", f))

print("Done.")
PYEOF

    python3 process_scores.py "${prefix}" "./chai1_output"

    cat > versions.yml << EOF
${task.process}:
  python: \$(python3 --version 2>/dev/null | sed 's/Python //g')
  chai_lab: \$(python3 -c "import chai_lab; print(chai_lab.__version__)" 2>/dev/null || echo "unknown")
  numpy: \$(python3 -c "import numpy; print(numpy.__version__)" 2>/dev/null || echo "unknown")
  torch: \$(python3 -c "import torch; print(torch.__version__)" 2>/dev/null || echo "unknown")
EOF
    """
}
