@echo off
REM Scaling test: does larger B fill the GPU?
REM Runs the same grid at chunk 6/24/96 (GPU only) and prints ms/case for each.
REM If ms/case drops as chunk grows, GPU utilization improved.
REM Usage: run_gpu_scaling.bat            (default grid100.csv, needs >=96 rows)
REM        run_gpu_scaling.bat mygrid.csv
set GRID=%1
if "%GRID%"=="" set GRID=grid100.csv
echo ==================== chunk=6 ====================
python gpu_gate.py --grid %GRID% --data data --bands 555 --chunk 6 --gpu-only --nt-atm 60 --m-max 4 --max-it-water 150
echo ==================== chunk=24 ====================
python gpu_gate.py --grid %GRID% --data data --bands 555 --chunk 24 --gpu-only --nt-atm 60 --m-max 4 --max-it-water 150
echo ==================== chunk=96 ====================
python gpu_gate.py --grid %GRID% --data data --bands 555 --chunk 96 --gpu-only --nt-atm 60 --m-max 4 --max-it-water 150
echo ==================== done: compare ms/case above ====================
