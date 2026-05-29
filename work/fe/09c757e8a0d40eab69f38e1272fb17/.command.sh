#!/usr/bin/env bash -C -e -u -o pipefail
mkdir -p proteinmpnn_output

INPUT_PDB=$(readlink -f backbone_design_2.pdb)

docker run --rm --gpus all         -u $(id -u):$(id -g)         -v $(dirname $INPUT_PDB):/input         -v $PWD/proteinmpnn_output:/output         docker.io/rosettacommons/proteinmpnn:latest         --pdb_path /input/$(basename $INPUT_PDB)         --out_folder /output         --num_seq_per_target 1         --sampling_temp 0.1         --model_name v_48_020

cp -r proteinmpnn_output/seqs .
