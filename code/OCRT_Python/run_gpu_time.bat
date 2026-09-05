@echo off
REM Measure real production per-case time at FULL params (GPU only, one band).
REM Full params = nt-atm 400, m-max 16, max-it-water 500 (gpu_gate defaults).
REM Usage: run_gpu_time.bat <grid.csv> <chunk>
REM   e.g. run_gpu_time.bat grid100.csv 24
python gpu_gate.py --grid %1 --data data --bands 555 --chunk %2 --gpu-only
