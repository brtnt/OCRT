$ErrorActionPreference = 'Stop'
$PkgWin = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$PkgWsl = (wsl wslpath -a "$PkgWin").Trim()
wsl bash -lc "cd '$PkgWsl' && python3 runner/run_campaign.py --runtime-root '$PkgWsl/runtime/MIGRATION_PKG_2026-08-19' --matrix '$PkgWsl/matrices/PURE_RAYLEIGH_15RUN_MATRIX.csv' --output-dir '$PkgWsl/results/pure_rayleigh_5band_pssa' --workers 16"
if ($LASTEXITCODE -ne 0) { throw "WSL campaign failed: $LASTEXITCODE" }
Write-Host "Results: $PkgWin\results\pure_rayleigh_5band_pssa"
