#!/usr/bin/env bash
set -euo pipefail
PKG_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
PARTS_DIR=${1:-$(dirname "$PKG_ROOT")}
RUNTIME_DIR=${2:-$PKG_ROOT/runtime}
python3 "$PKG_ROOT/runner/setup_runtime.py" --parts-dir "$PARTS_DIR" --runtime-dir "$RUNTIME_DIR" --force
