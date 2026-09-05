@echo off
python gpu_gate.py --grid "%~1" --data data --bands 443,555,865
