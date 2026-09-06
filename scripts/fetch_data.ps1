# ----------------------------------------------------------------------------
# fetch_data.ps1 — download the OCRT Mie tables (GitHub Release "data-v1") and
# install them into
#     code\OCRT_C\inputs\      (C reference implementation)
#     code\OCRT_Python\data\   (Python batch package; same 181 files)
#
# Usage (PowerShell 5.1 or 7):   powershell -ExecutionPolicy Bypass -File scripts\fetch_data.ps1 [-NoVerify]
# About 1.3 GB download, 8.6 GB installed.  Existing archives in package\release_data-v1 are reused.
# ----------------------------------------------------------------------------
param([switch]$NoVerify)
$ErrorActionPreference = "Stop"
$Root  = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$Tag   = "data-v1"
$Base  = "https://github.com/brtnt/OCRT/releases/download/$Tag"
$Assets = @("ocrt_mie_opac16_v1.zip", "ocrt_mie_hydrosol_v1.zip",
            "ocrt_mie_ahmad2010_paper_v1.zip", "ocrt_mie_ahmad2010_accurt_v1.zip")
$Sums  = Join-Path $Root "scripts\fetch_data.sha256"
$CIn   = Join-Path $Root "code\OCRT_C\inputs"
$PyIn  = Join-Path $Root "code\OCRT_Python\data"
$Cache = Join-Path $Root "package\release_$Tag"
New-Item -ItemType Directory -Force -Path $Cache, $CIn, $PyIn | Out-Null

[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
foreach ($z in $Assets) {
    $dst = Join-Path $Cache $z
    if (-not (Test-Path $dst) -or (Get-Item $dst).Length -eq 0) {
        Write-Host "downloading $z"
        Invoke-WebRequest -Uri "$Base/$z" -OutFile "$dst.part" -UseBasicParsing
        Move-Item -Force "$dst.part" $dst
    } else { Write-Host "reusing $dst" }
}

if (-not $NoVerify) {
    Write-Host "verifying archives"
    $expected = @{}
    Get-Content $Sums | ForEach-Object { $p = $_ -split '\s+', 2; if ($p.Count -eq 2) { $expected[$p[1].TrimStart('*')] = $p[0].ToLower() } }
    foreach ($z in $Assets) {
        $h = (Get-FileHash -Algorithm SHA256 (Join-Path $Cache $z)).Hash.ToLower()
        if ($h -ne $expected[$z]) { throw "SHA-256 mismatch for $z" }
        Write-Host "  $z OK"
    }
}

foreach ($z in $Assets) {
    Write-Host "unpacking $z -> $CIn"
    Expand-Archive -Path (Join-Path $Cache $z) -DestinationPath $CIn -Force
}
$n = (Get-ChildItem -Path $CIn -Recurse -Filter *.mie -File).Count
Write-Host "installed $n Mie tables in $CIn (expected 181)"

Write-Host "populating $PyIn"
Get-ChildItem -Path $CIn -Recurse -Filter *.mie -File | ForEach-Object {
    $rel = $_.FullName.Substring($CIn.Length).TrimStart('\')
    $dst = Join-Path $PyIn $rel
    New-Item -ItemType Directory -Force -Path (Split-Path -Parent $dst) | Out-Null
    Copy-Item -Force $_.FullName $dst
}
Write-Host "done."
