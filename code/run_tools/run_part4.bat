@echo off
chcp 65001 > nul
setlocal

rem ==========================================================================
rem  OCRT polarization sensitivity - PART 4
rem  Expansion that does NOT touch wind speed or solar zenith angle.
rem
rem  Three sub-stages, run in order. Each one reuses everything already
rem  computed, so nothing is repeated.
rem
rem    p4a  aerosol models  3 -> 8            +375 runs
rem    p4b  water cases    27 -> 64           +555 runs
rem    p4   N-branch AOD  0.1 -> 0.05/0.1/0.3 +1890 runs
rem
rem  Total new: 2820 runs on top of the 615 already done.
rem
rem  Point --out at the SAME folder used for the first run so the existing
rem  615 results are found and skipped. Default is <package root>\star_run.
rem
rem  Usage (PowerShell needs the leading .\ ):
rem    .\run_part4.bat --dry-run     show the plan only
rem    .\run_part4.bat               run all three sub-stages in order
rem    .\run_part4.bat p4a           run one sub-stage only
rem
rem  Extra options after the stage name: --workers N  --out PATH  --part K/N
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

rem If the first argument names a sub-stage, run only that one.
if /I "%~1"=="p4a" goto :one
if /I "%~1"=="p4b" goto :one
if /I "%~1"=="p4"  goto :one

echo ==== stage p4a : aerosol 3 -^> 8 models ====
%PY% "%SCRIPT%" --stage p4a --workers 24 %*
if errorlevel 1 goto :stop

echo.
echo ==== stage p4b : water 27 -^> 64 cases ====
%PY% "%SCRIPT%" --stage p4b --workers 24 %*
if errorlevel 1 goto :stop

echo.
echo ==== stage p4 : N-branch AOD 0.05 / 0.1 / 0.3 ====
%PY% "%SCRIPT%" --stage p4 --workers 24 %*
goto :done

:one
set "ST=%~1"
shift
echo ==== stage %ST% only ====
%PY% "%SCRIPT%" --stage %ST% --workers 24 %1 %2 %3 %4 %5 %6 %7 %8 %9
goto :done

:stop
echo.
echo [STOP] a sub-stage reported a problem. Later sub-stages were not started.
goto :done

:nopython
echo [ERROR] Python was not found. Install Python 3.8 or newer,
echo         make sure "python" works in this window, then retry.

:done
echo.
pause
