#!/usr/bin/env bash
# Compatibility entry point. SIMPLE was removed in OCRT v1.18; run the
# canonical OCRT/CCRR/IOP water-branch contract smoke instead.
set -euo pipefail
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
exec "$HERE/smoke_cli_water_branch_v118.sh" "$@"
