/*
 * Run Chai-1
 */
process RUN_CHAI1 {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'
    container 'quay.io/nf-core/proteinfold_chai1:1.3.0'

    shell '/bin/bash', '-euo', 'pipefail'

    input:
    tuple val(meta), path(fasta)
    path  weights_dir

    output:
    path  ("raw/**")                              , emit: raw
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

records = []
for sf in score_files:
    idx  = int(sf.split("model_idx_")[1].replace(".npz", ""))
    data = np.load(sf)
    records.append((idx, float(data["aggregate_score"])))

records.sort(key=lambda x: x[1], reverse=True)

for rank, (idx, _score) in enumerate(records):
    src = os.path.join(output_dir, f"pred.model_idx_{idx}.cif")
    shutil.copy(src, f"raw/model_idx_{idx}_ranked_{rank}.cif")
    if rank == 0:
        shutil.copy(src, f"{prefix}_chai1.cif")

top_idx  = records[0][0]
top_data = np.load(os.path.join(output_dir, f"scores.model_idx_{top_idx}.npz"))
plddt    = top_data["per_atom_plddt"]

with open(f"{prefix}_plddt.tsv", "w") as fh:
    fh.write("atom_index\tplddt\n")
    for i, v in enumerate(plddt):
        fh.write(f"{i+1}\t{v:.4f}\n")

with open(f"{prefix}_ptm.tsv", "w") as fh:
    fh.write("model\tptm\taggregate_score\n")
    for rank, (idx, agg) in enumerate(records):
        d   = np.load(os.path.join(output_dir, f"scores.model_idx_{idx}.npz"))
        ptm = float(d["ptm"])
        fh.write(f"ranked_{rank}\t{ptm:.4f}\t{agg:.4f}\n")

try:
    with open(f"{prefix}_iptm.tsv", "w") as fh:
        fh.write("model\tiptm\n")
        for rank, (idx, _) in enumerate(records):
            d    = np.load(os.path.join(output_dir, f"scores.model_idx_{idx}.npz"))
            iptm = float(d["iptm"])
            fh.write(f"ranked_{rank}\t{iptm:.4f}\n")
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