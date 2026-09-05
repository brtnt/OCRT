$ErrorActionPreference = 'Stop'
$PkgWin = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$PartsWin = (Resolve-Path (Join-Path $PkgWin '..')).Path
$PkgWsl = (wsl wslpath -a "$PkgWin").Trim()
$PartsWsl = (wsl wslpath -a "$PartsWin").Trim()
wsl bash -lc "cd '$PkgWsl' && python3 runner/setup_runtime.py --parts-dir '$PartsWsl' --runtime-dir '$PkgWsl/runtime' --force"
if ($LASTEXITCODE -ne 0) { throw "WSL setup failed: $LASTEXITCODE" }
Write-Host "Setup complete: $PkgWin\runtime"
