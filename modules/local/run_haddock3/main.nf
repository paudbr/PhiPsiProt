process RUN_HADDOCK3 {

    tag "run_haddock3"
    label 'process_high'

    container "quay.io/nf-core/haddock:3.0.0"

    input:
    tuple val(candidate_id),
          path(ligand_pdb),
          path(receptor_pdb),
          path(restraints)

    output:
    tuple val(candidate_id), path("${candidate_id}_haddock.tsv"), emit: scores
    path "versions.yml", emit: versions

    script:
    """
    mkdir run
    cd run

    ln -s ${receptor_pdb} receptor.pdb
    ln -s ${ligand_pdb} ligand.pdb
    ln -s ${restraints} restraints.tbl

    cat <<EOF > haddock3.cfg
run_dir = "run"

molecules = [
    "receptor.pdb",
    "ligand.pdb"
]

restraints = "restraints.tbl"

[topoaa]
autohis = true

[rigidbody]
sampling = 500

[caprieval]
EOF

    haddock3 haddock3.cfg

    python <<PY
import glob, pandas as pd

f = glob.glob("run/analysis/**/capri_ss.tsv", recursive=True)[0]
df = pd.read_csv(f, sep="\\t")

score_col = [c for c in df.columns if "score" in c.lower()][0]
best = df.sort_values(score_col).iloc[0]

pd.DataFrame([{
    "candidate_id": "${candidate_id}",
    "score": best[score_col]
}]).to_csv("../${candidate_id}_haddock.tsv", sep="\\t", index=False)
PY

    cd ..
    """
}