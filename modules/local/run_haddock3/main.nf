process RUN_HADDOCK3 {
    tag "run_haddock3"
    label 'process_high'
    container "quay.io/nf-core/haddock:3.0.0"
    input:
    tuple val(candidate_id),
        path(receptor_pdb),
        path(ligand_pdb),
        path(restraints)
    output:
    tuple val(candidate_id), path("${candidate_id}_haddock.tsv"),         emit: scores
    tuple val(candidate_id), path("${candidate_id}_haddock_summary.tsv"), emit: summary
    path "versions.yml", emit: versions
    script:
    def task_process = task.process
    """
    mkdir -p run

    awk '/^ATOM|^HETATM/{sub(/^(.{21})./, substr(\$0,1,21)"A"); print} !/^ATOM|^HETATM/{print}' \\
        ${receptor_pdb} > run/receptor.pdb

    awk '/^ATOM|^HETATM/{sub(/^(.{21})./, substr(\$0,1,21)"B"); print} !/^ATOM|^HETATM/{print}' \\
        ${ligand_pdb} > run/ligand.pdb

    ln -s \$(realpath ${restraints}) run/restraints.tbl

    cat << EOF > run/haddock3.cfg
run_dir = "haddock_out"
molecules = [
    "receptor.pdb",
    "ligand.pdb"
]
[topoaa]
autohis = true
ncores = ${task.cpus}
[rigidbody]
sampling = 500
ambig_fname = "restraints.tbl"
ncores = ${task.cpus}
[caprieval]
ncores = ${task.cpus}
EOF

    cd run
    haddock3 haddock3.cfg

    python3 << PY
import glob, pandas as pd
f = glob.glob("haddock_out/analysis/**/capri_ss.tsv", recursive=True)[0]
df = pd.read_csv(f, sep="\\t")
df.insert(0, "candidate_id", "${candidate_id}")
score_col = [c for c in df.columns if "score" in c.lower()][0]
df_sorted = df.sort_values(score_col)

df_sorted.to_csv("../${candidate_id}_haddock.tsv", sep="\\t", index=False)

summary = pd.DataFrame([{
    "candidate_id":     "${candidate_id}",
    "best_score":       df_sorted[score_col].iloc[0],
    "top10_mean_score": df_sorted[score_col].head(10).mean(),
    "median_score":     df_sorted[score_col].median(),
    "n_poses":          len(df_sorted)
}])
summary.to_csv("../${candidate_id}_haddock_summary.tsv", sep="\\t", index=False)
PY

    cd ..

    cat << END_VERSIONS > versions.yml
"${task_process}":
    haddock3: \$(haddock3 --version 2>&1 | sed 's/HADDOCK3 version: //')
    python: \$(python3 --version | sed 's/Python //')
END_VERSIONS
    """
}