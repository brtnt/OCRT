@echo off
set OCRT_PY_GPU=0
python produce_grid.py --grid "%~1" --out "%~2" --data data --n-mu-water 24 --nt-atm 400 --m-max 16 --max-it-water 500 --chunk 32
