#!/usr/bin/env bash -C -e -u -o pipefail
mkdir -p rfdiffusion_output

INPUT_PDB=$(readlink -f test_interface.pdb)

docker run --rm --gpus all         -v $(dirname $INPUT_PDB):/input         -v $PWD/rfdiffusion_output:/output         docker.io/rosettacommons/rfdiffusion:latest         inference.output_prefix=/output/backbone_design         inference.input_pdb=/input/$(basename $INPUT_PDB)         inference.num_designs=1         'contigmap.contigs=[4-4]'         diffuser.partial_T=1

cp rfdiffusion_output/backbone_design_0.pdb .
cp rfdiffusion_output/backbone_design_0.trb .
