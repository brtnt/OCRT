# EAP phytoplankton option gate — 2026-07-26

- Default remains `micro`.
- `micro` and `eap_diatoms_centric` are the same regenerated D=6 um phase file.
- Omitting `--phyto` / `--phyto-group` preserves the previous usage.
- `pico`, `nano`, and every other non-default `eap_*` selection require
  `OCRT_ADVANCED=1`.
- The gate is applied to `ocrt_solve.py`, `ocrt_fullgrid.py`, and
  `produce_grid.py`.
- Programmatic low-level APIs retain access to the full 17-species catalog.

Linux example:

```bash
OCRT_ADVANCED=1 python ocrt_solve.py single ... \
  --phyto-group eap_synechococcus
```

Windows CMD:

```bat
set OCRT_ADVANCED=1
python ocrt_solve.py single ... --phyto-group eap_synechococcus
```
