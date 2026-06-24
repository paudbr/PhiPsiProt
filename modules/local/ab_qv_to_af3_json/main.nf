/*
 * AB_QV_TO_AF3_JSON
 * ──────────────────────────────────────────────────────────────────────────
 * Bridges rfantibody Quiver files → AlphaFold3 JSON input format.
 *
 * Reads the Quiver file from ProteinMPNN, extracts H[L] sequences for each
 * design, and pairs them with the antigen (chain T) sequence to produce
 * one AF3 JSON per design.
 *
 * AF3 JSON strategy (Watson et al. 2026):
 *   - Antibody H[L] chains: template-only, NO MSA (CDRs are de novo)
 *   - Antigen chain T:       full MSA + templates (standard AF3 behaviour)
 *   - Seeds: ab_af3_seeds per design, take max ipTM across seeds
 *
 * The existing RUN_ALPHAFOLD3 process is used as-is — no modifications needed.
 */
process AB_QV_TO_AF3_JSON {
    tag "${meta.id}"
    label 'process_single'

    container 'rfantibody_local:1.0.0'

    input:
    tuple val(meta), path(mpnn_quiver)   // Quiver from RUN_PROTEINMPNN_RFAB
    path  target_hlt_pdb                 // antigen PDB in HLT format (chain T)
    val   ab_type                        // 'VHH' | 'scFv'
    val   n_seeds                        // AF3 seeds per design (default 10)

    output:
    tuple val(meta), path("af3_inputs/*.json"), emit: jsons   // one JSON per design
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p af3_inputs extracted_pdbs

    # ── Step 1: Extract all PDBs from Quiver ──────────────────────────────
    qvextract ${mpnn_quiver} -o extracted_pdbs/

    # ── Step 2: Build one AF3 JSON per design ─────────────────────────────
    python3 <<'PYSCRIPT'
import json, os, glob, sys, string

aa3to1 = dict(ALA="A",ARG="R",ASN="N",ASP="D",CYS="C",GLN="Q",
              GLU="E",GLY="G",HIS="H",ILE="I",LEU="L",LYS="K",
              MET="M",PHE="F",PRO="P",SER="S",THR="T",TRP="W",
              TYR="Y",VAL="V")

def extract_chain_seq(pdb_path, chain_id):
    seq = {}
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[13:15].strip() == "CA":
                ch  = line[21]
                if ch != chain_id:
                    continue
                rno = int(line[22:26])
                aa  = aa3to1.get(line[17:20].strip(), "X")
                seq[rno] = aa
    return "".join(seq[k] for k in sorted(seq))

ab_type      = "${ab_type}"
n_seeds      = int("${n_seeds}")
target_pdb   = "${target_hlt_pdb}"
pdbs         = sorted(glob.glob("extracted_pdbs/*.pdb"))

chains_found = set()
with open(target_pdb) as fh:
    for line in fh:
        if line.startswith("ATOM") and line[13:15].strip() == "CA":
            chains_found.add(line[21])

antigen_chain = "T" if "T" in chains_found else list(chains_found)[0]
print(f"Antigen chain detected: {antigen_chain}", file=sys.stderr)
antigen_seq = extract_chain_seq(target_pdb, antigen_chain)
if not antigen_seq:
    print("ERROR: no antigen chain found in target PDB", file=sys.stderr)
    sys.exit(1)

print(f"Antigen sequence length: {len(antigen_seq)}", file=sys.stderr)
print(f"Processing {len(pdbs)} designs from Quiver...", file=sys.stderr)

for pdb_path in pdbs:
    design_id = os.path.splitext(os.path.basename(pdb_path))[0]
    # Sanitise for AF3 filename requirements (lowercase, no special chars)
    safe_id = design_id.lower().replace(" ", "_")
    allowed = set(string.ascii_lowercase + string.digits + "_-.")
    safe_id = "".join(c for c in safe_id if c in allowed)

    heavy_seq = extract_chain_seq(pdb_path, "H")
    light_seq = extract_chain_seq(pdb_path, "L") if ab_type == "scFv" else ""

    if not heavy_seq:
        print(f"WARNING: no chain H in {pdb_path}, skipping", file=sys.stderr)
        continue

    # Build sequence list
    # Antibody chains: unpairedMsa="" + pairedMsa="" (no MSA, CDRs are de novo)
    # Antigen chain: AF3 will use its DB MSA automatically (no override)
    sequences = [
        {
            "protein": {
                "id": "H",
                "sequence": heavy_seq,
                "unpairedMsa": "",
                "pairedMsa":   "",
                "templates":   []
            }
        }
    ]
    if light_seq:
        sequences.append({
            "protein": {
                "id": "L",
                "sequence": light_seq,
                "unpairedMsa": "",
                "pairedMsa":   "",
                "templates":   []
            }
        })
    sequences.append({
        "protein": {
            "id": "T",
            "sequence": antigen_seq
            # No MSA override → AF3 uses full DB search for antigen
        }
    })

    af3_json = {
        "name":       safe_id,
        "sequences":  sequences,
        "modelSeeds": list(range(n_seeds)),
        "dialect":    "alphafold3",
        "version":    1
    }

    out_path = f"af3_inputs/{safe_id}.json"
    with open(out_path, "w") as fh:
        json.dump(af3_json, fh, indent=2)

print(f"Written {len(list(glob.glob('af3_inputs/*.json')))} AF3 JSON files", file=sys.stderr)
PYSCRIPT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        rfantibody: \$(qvextract --version 2>/dev/null || echo "rfantibody_local:1.0.0")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p af3_inputs
    echo '{"name":"stub_design_0001","sequences":[{"protein":{"id":"H","sequence":"QVQLVESGG"}}],"modelSeeds":[0],"dialect":"alphafold3","version":1}' \
        > af3_inputs/stub_design_0001.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
