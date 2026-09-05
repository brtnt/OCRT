@echo off
chcp 65001 > nul
setlocal

rem ==========================================================================
rem  OCRT polarization sensitivity run - part 2 of 3
rem
rem  Put this file and run_star.py in the OCRT package root or its build\
rem  folder. The package root is the folder holding both "inputs" and "src".
rem  No path editing is needed.
rem
rem  Run the three parts on three machines, or one after another on one
rem  machine. Each part writes its own index_part2of3.csv, so they never
rem  overwrite each other.
rem
rem  Usage (PowerShell needs the leading .\ ):
rem    .\run_part2.bat --verify     check the executable first (do this once)
rem    .\run_part2.bat --dry-run    show the plan only
rem    .\run_part2.bat              run this part
rem
rem  Extra options: --workers N   --out PATH   --max-orders N
rem ==========================================================================

set "SCRIPT=%~dp0run_star.py"
if not exist "%SCRIPT%" (
    echo [ERROR] run_star.py not found next to this batch file.
    echo         Expected: %SCRIPT%
    pause
    exit /b 1
)

set "PY=python"
%PY% --version >nul 2>nul
if errorlevel 1 set "PY=py -3"
%PY% --version >nul 2>nul
if errorlevel 1 goto :nopython

%PY% "%SCRIPT%" --part 2/3 --workers 24 %*
goto :done

:nopython
echo [ERROR] Python was not found. Install Python 3.8 or newer,
echo         make sure "python" works in this window, then retry.

:done
echo.
pause
