# OCRT C safe kernel/invariant update — 2026-07-29

Version: `OCRT-v1.2-2026-07-29-KST-safe-kernel-invariant-update`

Applied exact, default-output-preserving items from `MIGRATION_WORKORDER_2026-07-29_KO.md`:

- Skip zero-quadrature-weight incoming columns in scalar/vector SOS source contractions.
  Disable only for regression with `OCRT_SOS_ZERO_COL_SKIP_OFF=1`.
- Hoist RAA-invariant Beer factors out of the native coupled angular-LUT cell loop.

Not adopted in the default code:

- water SOS geometric-tail acceleration: approximate and source details/guards were not supplied;
- `--precision standard` as the default: it intentionally changes optionless outputs;
- withdrawn multi-slot moment cache and Rww zero-column skip.
