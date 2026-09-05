#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
exec python3 scripts/test_stage2_coupling_clamp_fix.py --root "$ROOT" --binary build/ocrt
