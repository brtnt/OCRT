@echo off
setlocal enabledelayedexpansion
rem ==========================================================================
rem  Rebuild OCRT for Windows.
rem  Put this file and ocrt_win_compat.c in the OCRT package root
rem  (the folder that contains "src"), then run it.
rem  Requires MSYS2 MinGW64 gcc on PATH.
rem
rem  Usage:
rem    build_win.bat                 use -march=x86-64-v3  (AVX2, safe default)
rem    build_win.bat cascadelake     use -march=cascadelake (AVX-512, faster
rem                                  but only on CPUs that support it)
rem    build_win.bat native          let gcc detect this machine
rem ==========================================================================
cd /d "%~dp0"

set "ARCH=%~1"
if "%ARCH%"=="" set "ARCH=x86-64-v3"

if not exist "src" (
    echo [ERROR] Run this from the OCRT package root ^(the folder with "src"^).
    pause
    exit /b 1
)
if not exist "ocrt_win_compat.c" (
    echo [ERROR] ocrt_win_compat.c not found next to this batch file.
    pause
    exit /b 1
)
if not exist "build" mkdir build

echo Building with -march=%ARCH% ...

set "SRC="
for /r "%~dp0src" %%f in (*.c) do set "SRC=!SRC! "%%f""

gcc -std=c11 -O3 -march=%ARCH% -ffp-contract=fast -fassociative-math ^
    -fno-signed-zeros -fno-trapping-math -DOCRT_FAST_KERNELS -fopenmp ^
    -Isrc !SRC! "%~dp0ocrt_win_compat.c" -o "build\ocrt_%ARCH%.exe" -lm -static

if errorlevel 1 (
    echo [ERROR] build failed
) else (
    echo [OK] built build\ocrt_%ARCH%.exe
    echo      Verify it before use:
    echo        python build\run_star.py --exe build\ocrt_%ARCH%.exe --verify
)
pause
