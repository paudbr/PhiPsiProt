#!/usr/bin/env bash -C -e -u -o pipefail
mkdir -p rfdiffusion_output

INPUT_PDB=$(readlink -f 1PPN.pdb)

docker run --rm --gpus all         -v $(dirname $INPUT_PDB):/input         -v $PWD/rfdiffusion_output:/output         docker.io/rosettacommons/rfdiffusion:latest         inference.output_prefix=/output/binder_design         inference.input_pdb=/input/$(basename $INPUT_PDB)         inference.num_designs=3         'contigmap.contigs=[A1-212/0 30-30]'         'ppi.hotspot_res=[A25,A26,A27]'         diffuser.T=50         denoiser.noise_scale_ca=1.0         denoiser.noise_scale_frame=1.0

sudo chown -R $(id -u):$(id -g) rfdiffusion_output

cp rfdiffusion_output/binder_design_*.pdb .
cp rfdiffusion_output/binder_design_*.trb .
