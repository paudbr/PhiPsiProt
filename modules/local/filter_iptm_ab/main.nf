/*
 * FILTER_IPTM_AB
 * ──────────────────────────────────────────────────────────────────────────
 * Reads the iptm.tsv emitted by the existing RUN_ALPHAFOLD3 process and
 * applies the ipTM threshold from Watson et al. (Nature 2026):
 *
 *   VHH:  ipTM > 0.60  (AUC = 0.86, Extended Data Fig. 10a,b)
 *   scFv: ipTM > 0.85  (5/6 confirmed binders pass, Extended Data Fig. 10d)
 *
 * Takes the MAXIMUM ipTM over all seeds (AF3 JSON was built with n_seeds).
 *
 * Inputs come directly from RUN_ALPHAFOLD3 outputs — no extra AF3 run needed.
 */
process FILTER_IPTM_AB {
    tag "${meta.id}"
    label 'process_single'

    container 'docker.io/python:3.11-slim'

    input:
    tuple val(meta),
          path(iptm_tsv),      // *_iptm.tsv from RUN_ALPHAFOLD3
          path(top_cif)        // *_alphafold3.cif  (top ranked model)
    val   iptm_threshold       // 0.85 for scFv, 0.60 for VHH

    output:
    tuple val(meta), path("${meta.id}_af3_top.cif"), optional: true, emit: passed_cif
    tuple val(meta), path("${meta.id}_iptm_result.tsv"),              emit: result
    path  "versions.yml",                                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def threshold = iptm_threshold ?: params.ab_iptm_threshold ?: 0.85
    """
    python3 <<'PYSCRIPT'
import csv, shutil, sys

threshold = float("${threshold}")
sample_id = "${meta.id}"
iptm_tsv  = "${iptm_tsv}"
top_cif   = "${top_cif}"

# The iptm.tsv from RUN_ALPHAFOLD3 contains one row per seed
# Columns: model  ptm  iptm  aggregate_score   (ranked format)
# OR simple:  metric  value
iptm_vals = []

with open(iptm_tsv) as fh:
    content = fh.read().strip()

# Try ranked format first (model  ptm  iptm  aggregate_score)
lines = content.splitlines()
if lines:
    header = lines[0].split("\\t")
    if "iptm" in [h.lower() for h in header]:
        iptm_col = [h.lower() for h in header].index("iptm")
        for line in lines[1:]:
            parts = line.split("\\t")
            if len(parts) > iptm_col:
                try:
                    iptm_vals.append(float(parts[iptm_col]))
                except ValueError:
                    pass
    else:
        # Simple metric/value format
        for line in lines[1:]:
            parts = line.split("\\t")
            if len(parts) >= 2 and parts[0].lower() == "iptm":
                try:
                    iptm_vals.append(float(parts[1]))
                except ValueError:
                    pass

if not iptm_vals:
    print(f"WARNING: could not parse ipTM from {iptm_tsv}", file=sys.stderr)
    max_iptm = 0.0
    passed   = False
else:
    max_iptm = max(iptm_vals)
    passed   = max_iptm >= threshold

print(f"{sample_id}: max_iptm={max_iptm:.4f} n_seeds={len(iptm_vals)} "
      f"threshold={threshold} passed={passed}", file=sys.stderr)

# Write result TSV
with open(f"{sample_id}_iptm_result.tsv", "w") as fh:
    fh.write("id\\tmax_iptm\\tn_seeds\\tthreshold\\tpassed\\n")
    fh.write(f"{sample_id}\\t{max_iptm:.4f}\\t{len(iptm_vals)}\\t{threshold}\\t{passed}\\n")

if passed:
    shutil.copy(top_cif, f"{sample_id}_af3_top.cif")
    print(f"PASS → copied {top_cif} as {sample_id}_af3_top.cif", file=sys.stderr)
else:
    print(f"FAIL → {sample_id} rejected (ipTM {max_iptm:.4f} < {threshold})", file=sys.stderr)
PYSCRIPT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    cp ${top_cif} ${meta.id}_af3_top.cif
    printf "id\\tmax_iptm\\tn_seeds\\tthreshold\\tpassed\\n" > ${meta.id}_iptm_result.tsv
    printf "${meta.id}\\t0.9100\\t10\\t${iptm_threshold}\\ttrue\\n" >> ${meta.id}_iptm_result.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
