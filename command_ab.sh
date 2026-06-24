#!/usr/bin/env bash
# run_antibody_design.sh — RFdiffusion → ProteinMPNN → AF3 ipTM → HADDOCK3

# ── VHH example (influenza HA, ipTM > 0.60) ──────────────────────────────
nextflow run . \
    -profile docker \
    --mode              antibody_design \
    --outdir            /mnt/alphafold3_db/test_ab/vhh_ha \
    --input             samplesheet_papain.csv \
    --ab_type               VHH \
    --ab_target_pdb         /mnt/alphafold3_db/test_ab/flu_HA_truncated.pdb \
    --ab_framework_pdb      /mnt/alphafold3_db/test_ab/h-NbBCII10.pdb \
    --ab_hotspot_residues   "B131,B132,B133,B156,B157,B158" \
    --ab_design_loops       "H1:7,H2:6,H3:5-13" \
    \
    --ab_num_designs        10\
    --ab_seqs_per_struct    4 \
    --ab_mpnn_temperature   0.2 \
    \
    --ab_iptm_threshold     0.60 \
    --ab_af3_seeds          10 \
    --ab_run_haddock        true \
    --use_gpu               true \
    -resume
