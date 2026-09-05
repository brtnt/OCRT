#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
for f in Phytoplankton.mie TSM_mineral.mie Calcareous_sand_from_OSOAA_phase_443_550.mie; do
  test ! -e "$ROOT/inputs/$f" || { echo "FAIL: legacy Mie remains: $f" >&2; exit 1; }
done
echo "PASS: legacy narrow Mie files absent"
