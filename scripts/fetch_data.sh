#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# fetch_data.sh — download the OCRT Mie tables (GitHub Release "data-v1") and
# install them into
#     code/OCRT_C/inputs/      (C reference implementation)
#     code/OCRT_Python/data/   (Python batch package; same 181 files)
#
# Usage:  bash scripts/fetch_data.sh [--no-verify] [--link]
#   --no-verify   skip the SHA-256 check of the downloaded archives
#   --link        hard-link the Python copy instead of copying (same filesystem)
#
# Requires: curl, unzip, sha256sum.  About 1.3 GB download, 8.6 GB installed
# (4.3 GB per tree).  Re-running is safe: existing archives are reused.
# ----------------------------------------------------------------------------
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TAG="data-v1"
BASE="https://github.com/brtnt/OCRT/releases/download/${TAG}"
ASSETS=(ocrt_mie_opac16_v1.zip ocrt_mie_hydrosol_v1.zip
        ocrt_mie_ahmad2010_paper_v1.zip ocrt_mie_ahmad2010_accurt_v1.zip)
SUMS="${ROOT}/scripts/fetch_data.sha256"
C_IN="${ROOT}/code/OCRT_C/inputs"
PY_IN="${ROOT}/code/OCRT_Python/data"
CACHE="${ROOT}/package/release_${TAG}"
VERIFY=1; LINK=0
for a in "$@"; do
  case "$a" in
    --no-verify) VERIFY=0 ;;
    --link) LINK=1 ;;
    *) echo "unknown option: $a" >&2; exit 2 ;;
  esac
done

mkdir -p "$CACHE" "$C_IN" "$PY_IN"
for z in "${ASSETS[@]}"; do
  if [ ! -s "$CACHE/$z" ]; then
    echo "downloading $z"
    curl -L --fail --retry 5 --retry-delay 10 -o "$CACHE/$z.part" "$BASE/$z"
    mv "$CACHE/$z.part" "$CACHE/$z"
  else
    echo "reusing $CACHE/$z"
  fi
done

if [ "$VERIFY" = 1 ]; then
  echo "verifying archives"
  (cd "$CACHE" && sha256sum -c "$SUMS")
fi

for z in "${ASSETS[@]}"; do
  echo "unpacking $z -> $C_IN"
  unzip -q -o "$CACHE/$z" -d "$C_IN"
done

n=$(find "$C_IN" -name '*.mie' -type f | wc -l)
echo "installed $n Mie tables in $C_IN (expected 181)"

echo "populating $PY_IN"
(cd "$C_IN" && find . -name '*.mie' -type f -print0) | while IFS= read -r -d '' f; do
  mkdir -p "$PY_IN/$(dirname "$f")"
  if [ "$LINK" = 1 ]; then ln -f "$C_IN/$f" "$PY_IN/$f"; else cp -f "$C_IN/$f" "$PY_IN/$f"; fi
done
echo "done. Optional per-file check: (cd code/OCRT_C/inputs && sha256sum -c ../../../scripts/MIE_SHA256SUMS_${TAG}.txt)"
