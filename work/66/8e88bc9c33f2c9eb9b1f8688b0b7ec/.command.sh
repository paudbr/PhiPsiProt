#!/usr/bin/env bash -C -e -u -o pipefail
echo "candidate_id,mode,method,pdb,binder_length,status" > binder_candidates.csv

for pdb in binder_design_0.pdb; do
    id=$(basename "$pdb" .pdb)
    length=$(grep '^ATOM' "$pdb" | awk '{print $5,$6}' | sort -u | wc -l)

    echo "$id,binder,rfdiffusion,$(basename $pdb),$length,ready_for_sequence_design" >> binder_candidates.csv
done
