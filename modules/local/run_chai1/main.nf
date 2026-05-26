/*
 * Run Chai-1
 */
process RUN_CHAI1 {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'
    container "nf-core/proteinfold_chai1:1.0.0"

    input:
    tuple val(meta), path(fasta)
    path weights_dir

    output:
    path ("raw/**")                                     , emit: raw
    tuple val(meta), path ("${meta.id}_chai1.cif")      , emit: top_ranked_cif
    tuple val(meta), path ("raw/*ranked_*.cif")         , emit: cif
    tuple val(meta), path ("${meta.id}_plddt.tsv")      , emit: multiqc
    tuple val(meta), path ("${meta.id}_ptm.tsv")        , emit: ptms
    tuple val(meta), path ("${meta.id}_iptm.tsv")       , optional: true, emit: iptms
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Local RUN_CHAI1 module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Point Chai-1 to the locally staged weights (avoids runtime download)
    export CHAI_DOWNLOADS_DIR=\$PWD/weights

    # Run Chai-1 — no MSA databases needed
    chai fold \\
        ${fasta} \\
        ./chai1_output \\
        $args

    mkdir -p raw

    # Parse scores, rank models, emit TSVs
    python3 << 'PYEOF'
import numpy as np
import glob, os, shutil

prefix     = "${prefix}"
output_dir = "./chai1_output"

score_files = sorted(glob.glob(f"{output_dir}/scores.model_idx_*.npz"))
if not score_files:
    raise FileNotFoundError("No Chai-1 score files found in output directory")

# Load aggregate_score for each seed
records = []
for sf in score_files:
    idx  = int(sf.split("model_idx_")[1].replace(".npz", ""))
    data = np.load(sf)
    records.append((idx, float(data["aggregate_score"]), data))

# Sort by aggregate_score descending → rank 0 is best
records.sort(key=lambda x: x[1], reverse=True)

# ── Ranked CIF files ──────────────────────────────────────────────────────────
for rank, (idx, score, _) in enumerate(records):
    src = f"{output_dir}/pred.model_idx_{idx}.cif"
    shutil.copy(src, f"raw/model_idx_{idx}_ranked_{rank}.cif")
    if rank == 0:
        shutil.copy(src, f"{prefix}_chai1.cif")

# ── pLDDT TSV (per-atom, top-ranked model) ───────────────────────────────────
top_idx  = records[0][0]
top_data = np.load(f"{output_dir}/scores.model_idx_{top_idx}.npz")
plddt    = top_data["per_atom_plddt"]
with open(f"{prefix}_plddt.tsv", "w") as fh:
    fh.write("atom_index\tplddt\n")
    for i, v in enumerate(plddt):
        fh.write(f"{i + 1}\t{v:.4f}\n")

# ── PTM TSV (all models, ranked) ─────────────────────────────────────────────
with open(f"{prefix}_ptm.tsv", "w") as fh:
    fh.write("model\tptm\taggregate_score\n")
    for rank, (idx, agg, _) in enumerate(records):
        ptm = float(np.load(f"{output_dir}/scores.model_idx_{idx}.npz")["ptm"])
        fh.write(f"ranked_{rank}\t{ptm:.4f}\t{agg:.4f}\n")

# ── iPTM TSV (optional — only meaningful for multimers) ──────────────────────
try:
    with open(f"{prefix}_iptm.tsv", "w") as fh:
        fh.write("model\tiptm\n")
        for rank, (idx, _, __) in enumerate(records):
            iptm = float(np.load(f"{output_dir}/scores.model_idx_{idx}.npz")["iptm"])
            fh.write(f"ranked_{rank}\t{iptm:.4f}\n")
except Exception:
    pass  # monomers may not expose iptm; output is optional

# ── Move raw Chai-1 outputs into raw/ for save_intermediates ─────────────────
for f in os.listdir(output_dir):
    shutil.move(f"{output_dir}/{f}", f"raw/{f}")
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        chai_lab: \$(python3 -c "import chai_lab; print(chai_lab.__version__)" 2>/dev/null || echo "unknown")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)" 2>/dev/null || echo "unknown")
        torch: \$(python3 -c "import torch; print(torch.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p raw
    touch ${prefix}_chai1.cif
    touch raw/${prefix}_ranked_0.cif
    touch raw/${prefix}_ranked_1.cif
    touch raw/${prefix}_ranked_2.cif
    touch raw/${prefix}_ranked_3.cif
    touch raw/${prefix}_ranked_4.cif
    touch ${prefix}_plddt.tsv
    touch ${prefix}_ptm.tsv
    touch ${prefix}_iptm.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>/dev/null | sed 's/Python //g' || echo "unknown")
        chai_lab: \$(python3 -c "import chai_lab; print(chai_lab.__version__)" 2>/dev/null || echo "unknown")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)" 2>/dev/null || echo "unknown")
        torch: \$(python3 -c "import torch; print(torch.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}