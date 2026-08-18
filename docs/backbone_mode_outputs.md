# Backbone Mode Outputs

## Mode

`--mode backbone`

## Purpose

Generate stabilized/diversified protein backbones using RFdiffusion and design sequences using ProteinMPNN.

## Workflow

```text
input PDB
→ RFdiffusion partial diffusion
→ generated backbone PDB
→ ProteinMPNN sequence design
→ final backbone candidates
