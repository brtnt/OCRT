@echo off
chcp 65001 > nul
setlocal

rem ==========================================================================
rem  OCRT polarization sensitivity - PART 4-2  (addition to part 4)
rem
rem  This does NOT redo anything. It only ADDS two aerosol optical depth
rem  levels to the S branch: 0.02 and 0.60, on top of the existing
rem  0.05 / 0.10 / 0.20 / 0.30 / 0.40.
rem
rem  Why: the S metric compares aerosol models at matched TOA intensity.
rem  With five AOD levels the models only share 8.5 percent of the intensity
rem  range, so the comparison often cannot be made at all. Adding 0.02 and
rem  0.60 widens that shared range to 30.9 percent. Measured at 555 nm,
rem  solar zenith 50 degrees.
rem
rem  New runs: 240 only. About 40 minutes on 24 cores.
rem
rem  Run this AFTER run_part4.bat has finished. Point --out at the SAME
rem  folder as before so everything already computed is reused.
rem
rem  Usage (PowerShell needs the leading .\ ):
rem    .\run_part4-2.bat --dry-run
rem    .\run_part4-2.bat
rem ==========================================================================

set "SCRIPT=%~dp0run_star.py"
if not exist "%SCRIPT%" (
    echo [ERROR] run_star.py not found next to this batch file.
    pause
    exit /b 1
)

set "PY=python"
%PY% --version >nul 2>nul
if errorlevel 1 set "PY=py -3"
%PY% --version >nul 2>nul
if errorlevel 1 goto :nopython

echo ==== stage p4c : S-branch AOD 0.02 and 0.60 added ====
%PY% "%SCRIPT%" --stage p4c --workers 24 %*
goto :done

:nopython
echo [ERROR] Python was not found. Install Python 3.8 or newer,
echo         make sure "python" works in this window, then retry.

:done
echo.
pause
