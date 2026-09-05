@echo off
chcp 65001 > nul
setlocal

rem ==========================================================================
rem  OCRT polarization sensitivity - PART 6
rem  Aerosol degeneracy grid, and the water-free reference atmosphere.
rem
rem  Why this replaces part 5:
rem  Single-scattering albedo and particle size both shift the spectral slope
rem  of intensity reflectance, so intensity alone cannot separate them. The
rem  earlier eight aerosol models were picked to spread albedo evenly, which
rem  puts them on a diagonal in the (size, albedo) plane and never creates the
rem  degenerate case. This part adds a 4 x 4 grid where fine-mode fraction
rem  sets size and relative humidity sets absorption, so the two vary
rem  independently.
rem
rem  Two sub-stages, in order:
rem
rem    p6   16 new aerosol models, coupled ocean       +1680 runs   ~4.3 h
rem    p6a  water-free reference, all 24 models        +2520 runs   ~0.4 h
rem
rem  The eight earlier models are kept, so nothing already computed is lost.
rem  Run part 4 and part 4-2 first. Point --out at the SAME folder as before.
rem
rem  Usage (PowerShell needs the leading .\ ):
rem    .\run_part6.bat --dry-run     show the plan only
rem    .\run_part6.bat               run both sub-stages in order
rem    .\run_part6.bat p6            run one sub-stage only
rem
rem  Extra options: --workers N   --out PATH   --part K/N
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

if /I "%~1"=="p6"  goto :one
if /I "%~1"=="p6a" goto :one

echo ==== stage p6 : aerosol degeneracy grid, coupled ocean ====
%PY% "%SCRIPT%" --stage p6 --workers 24 %*
if errorlevel 1 goto :stop

echo.
echo ==== stage p6a : water-free reference atmosphere ====
%PY% "%SCRIPT%" --stage p6a --workers 24 %*
goto :done

:one
set "ST=%~1"
shift
echo ==== stage %ST% only ====
%PY% "%SCRIPT%" --stage %ST% --workers 24 %1 %2 %3 %4 %5 %6 %7 %8 %9
goto :done

:stop
echo.
echo [STOP] the first sub-stage reported a problem. The second was not started.
goto :done

:nopython
echo [ERROR] Python was not found. Install Python 3.8 or newer,
echo         make sure "python" works in this window, then retry.

:done
echo.
pause
