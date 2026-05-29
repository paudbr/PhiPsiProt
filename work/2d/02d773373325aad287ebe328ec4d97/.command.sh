#!/usr/bin/env bash -C -e -u -o pipefail
mkdir -p ligandmpnn_output

INPUT_PDB=$(readlink -f binder_design_1.pdb)

docker run --rm --gpus all         -v $(dirname $INPUT_PDB):/input         -v $PWD/ligandmpnn_output:/output         docker.io/rosettacommons/ligandmpnn:latest         --model_type ligand_mpnn         --pdb_path /input/$(basename $INPUT_PDB)         --out_folder /output         --seed 111         --temperature 0.1

sudo chown -R $(id -u):$(id -g) ligandmpnn_output

cp -r ligandmpnn_output/seqs .
cp -r ligandmpnn_output/backbones .
