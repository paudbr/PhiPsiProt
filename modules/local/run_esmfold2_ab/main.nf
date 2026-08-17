/*
 * RUN_ESMFOLD2_AB  —  batch mode
 * ──────────────────────────────────────────────────────────────────────────
 * Recibe TODOS los PDBs de diseño de una vez, carga el modelo una sola vez
 * y procesa en batch. Evita el OOM de cargar ESMC-6B 40 veces.
 *
 * Logic in bin/run_esmfold2_ab.py
 */
process RUN_ESMFOLD2_AB {
    tag "esmfold2_batch"
    label 'process_high_memory'

    container 'quay.io/esmfold2_local:1.0.0'

    input:
    tuple val(meta),
          path(designs_dir),     // directorio con todos los PDBs de diseño
          path(target_pdb),      // antígeno truncado
          path(framework_pdb)    // framework HLT PDB (CDR REMARKs)

    output:
    tuple val(meta), path("esm2_all_metrics.tsv"), emit: metrics
    tuple val(meta), path("esm2_cifs/"),           emit: cifs
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def antigen_chain = params.ab_antigen_chain ?: 'B'
    def ab_type       = params.ab_type          ?: 'VHH'
    def hf_home       = "/mnt/alphafold3_db/ESMFold2/hub"
    def bin_dir       = "${projectDir}/bin"
    def num_loops     = params.ab_esm_num_loops  ?: 1
    def num_steps     = params.ab_esm_num_steps  ?: 10
    """
    export HF_HOME=${hf_home}
    export HUGGINGFACE_HUB_CACHE=${hf_home}
    export TRANSFORMERS_OFFLINE=1
    export HF_HUB_OFFLINE=1
    export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

    mkdir -p esm2_cifs

    run_esmfold2_ab.py \\
        --designs-dir   \$(readlink -f ${designs_dir}) \\
        --target        \$(readlink -f ${target_pdb}) \\
        --framework     \$(readlink -f ${framework_pdb}) \\
        --antigen-chain ${antigen_chain} \\
        --ab-type       ${ab_type} \\
        --out-tsv       esm2_all_metrics.tsv \\
        --out-cif-dir   esm2_cifs/ \\
        --num-loops     ${num_loops} \\
        --num-steps     ${num_steps} \\
        --seed          42

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esmfold2: "biohub/ESMFold2-Fast"
        python: \$(python3.12 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p esm2_cifs
    printf "id\tiptm_ht\tiptm_global\tpae_global\tpae_interface\tplddt_global\tplddt_cdr\trmsd_global\trmsd_cdr\n" \
        > esm2_all_metrics.tsv
    printf "stub_design_0\t0.45\t0.40\t8.5\t10.2\t0.75\t0.72\t1.8\t1.5\n" \
        >> esm2_all_metrics.tsv
    touch esm2_cifs/stub_design_0_esm2.cif
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esmfold2: stub
    END_VERSIONS
    """
}
