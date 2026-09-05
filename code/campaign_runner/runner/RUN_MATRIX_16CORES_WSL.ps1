param(
    [Parameter(Mandatory=$true)][string]$Matrix,
    [Parameter(Mandatory=$true)][string]$OutputDir
)
$ErrorActionPreference = 'Stop'
$PkgWin = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$MatrixWin = (Resolve-Path $Matrix).Path
$PkgWsl = (wsl wslpath -a "$PkgWin").Trim()
$MatrixWsl = (wsl wslpath -a "$MatrixWin").Trim()
# Output directory may not exist yet; convert its parent and append the leaf.
$OutFull = [System.IO.Path]::GetFullPath($OutputDir)
$OutParent = Split-Path $OutFull -Parent
$OutLeaf = Split-Path $OutFull -Leaf
New-Item -ItemType Directory -Force -Path $OutParent | Out-Null
$OutParentWsl = (wsl wslpath -a "$OutParent").Trim()
$OutWsl = "$OutParentWsl/$OutLeaf"
wsl bash -lc "cd '$PkgWsl' && python3 runner/run_campaign.py --runtime-root '$PkgWsl/runtime/MIGRATION_PKG_2026-08-19' --matrix '$MatrixWsl' --output-dir '$OutWsl' --workers 16"
if ($LASTEXITCODE -ne 0) { throw "WSL campaign failed: $LASTEXITCODE" }
