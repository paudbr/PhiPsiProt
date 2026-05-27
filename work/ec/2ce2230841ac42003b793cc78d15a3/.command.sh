#!/usr/bin/env bash -C -e -u -o pipefail
filter_ddg.py         --input candidates.csv         --output filtered_candidates.csv         --threshold 2.0
