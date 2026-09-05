#!/usr/bin/env bash
python3 subsample_grid.py --in full_grid_design_v1.csv --n "${1:-100}" --out "draft_grid_${1:-100}.csv"
