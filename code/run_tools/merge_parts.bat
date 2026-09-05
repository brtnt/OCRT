@echo off
chcp 65001 > nul
setlocal

rem ==========================================================================
rem  Merge the three parts into one dataset.
rem  Edit the three folder paths below to match where each part was run,
rem  then run this file. Use --copy to also gather the result CSV files.
rem ==========================================================================

set "PART1=D:\temp\star_run"
set "PART2=D:\temp\star_run"
set "PART3=D:\temp\star_run"
set "MERGED=D:\temp\star_merged"

set "SCRIPT=%~dp0merge_parts.py"
if not exist "%SCRIPT%" (
    echo [ERROR] merge_parts.py not found next to this batch file.
    pause
    exit /b 1
)

set "PY=python"
%PY% --version >nul 2>nul
if errorlevel 1 set "PY=py -3"

%PY% "%SCRIPT%" "%PART1%" "%PART2%" "%PART3%" --out "%MERGED%" %*

echo.
pause
