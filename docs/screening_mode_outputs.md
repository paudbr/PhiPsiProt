# Screening Mode Outputs

## Mode

`--mode screening`

## Main purpose

Generate saturation mutagenesis candidates from an input protein structure.

## Residue selection modes

- `manual`: mutate user-provided positions
- `all`: mutate all residues in the target chain
- `pocket`: mutate residues close to a ligand
- `interface`: mutate residues close to another protein chain

## Important outputs

### `residue_selection/selected_positions.txt`

Comma-separated list of selected mutable positions.

### `residue_selection/selection_summary.csv`

Metadata describing how mutable residues were selected.

### `screening/candidates.csv`

All generated mutation candidates.

### `screening/candidates.fasta`

FASTA sequences for all generated mutants.

### `filtering/filtered_candidates.csv`

Candidates after ΔΔG filtering.

Currently ΔΔG is a placeholder until PyRosetta/FoldX scoring is integrated.

### `filtering/filtered_candidates.fasta`

FASTA sequences for filtered candidates.

### `filtering/screening_samplesheet.csv`

Samplesheet for structural validation.

Format:

```csv
id,fasta
screening_000001|AA1C,filtered_candidates.fasta
