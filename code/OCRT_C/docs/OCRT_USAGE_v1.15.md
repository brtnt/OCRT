# OCRT v1.15 TSM usage

## Basic TSM run

```bash
./build/v2_solver_vk \
  --surface ocean --wind-speed 0 \
  --simple-mode \
  --simple-chl 0 \
  --simple-min 0.01 \
  --simple-adom440 0 \
  --tsm-species red_clay \
  --sza 30 --vza 0 --raa 90 \
  --wavelength 443 --pressure 0
```

`--simple-min` is dry-weight TSM in `g m^-3`. If `--tsm-species` is omitted,
`red_clay` is used. The available alternatives are `brown_earth`,
`yellow_clay`, and `calcareous_sand`.

The selected Ahn `.mie` phase is loaded automatically. Do not combine a
positive TSM input with `--water-mie-phase`, `--ccrr-phase-moments`, or
`--ccrr-particle-phase-lut`.

## Installed-data override

```bash
OCRT_TSM_DIR=/absolute/path/to/tsm_ahn ./build/v2_solver_vk ...
```

The directory must contain all four `*_AHN.mie` files. The `astarmin_*` and
`bstarmin_*` audit tables are optional because OCRT can read `a*` and `b*`
directly from the Mie bulk spectral blocks.
