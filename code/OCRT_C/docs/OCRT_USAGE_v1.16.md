# OCRT v1.16 TSM usage

## Build

Use a deterministic source-file order so a rebuild with the same compiler and
flags reproduces the packaged executable byte-for-byte.

```bash
cd ocrt
gcc -std=c11 -O3 -march=native -ffp-contract=fast -fassociative-math \
  -fno-signed-zeros -fno-trapping-math -DOCRT_FAST_KERNELS -fopenmp -Isrc \
  $(find src -name '*.c' -print | LC_ALL=C sort) \
  -o build/v2_solver_vk -lm
```

Run the executable from the `ocrt/` directory. Water-IOP, gas, atmosphere and
TSM data are resolved relative to that working directory. `OMP_NUM_THREADS=1`
is the reference setting for bitwise regression and performance measurement.

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

## Phase-cache diagnostic

Normal runs require no cache option. For development verification only:

```bash
OCRT_DUMP_SOS_PHASE_CACHE=1 ./build/v2_solver_vk ...
```

One `SOS_PHASE_CACHE_BUILD` record is emitted for each finalized Fourier mode.
The record count must depend on the number of Fourier modes, not on the number
of SOS scattering orders. `scripts/regression_tsm_phase_cache_v116.sh` checks
this invariant with water SOS caps 2 and 100 and also compares the prepared and
legacy fallback source operators bit-for-bit on synthetic I/Q/U fields.
