# OCRT black Fresnel ocean surface naming validation

Date: 2026-07-21

## Decision

The previous high-level name `coxmunk` was not an accurate name for the
surface boundary condition. It conflated:

1. the boundary condition — a rough Fresnel air-water interface with a
   perfectly black ocean below it, so no water-leaving radiance is present; and
2. the slope-variance law used inside the rough-surface microfacet kernel.

The canonical name is now `black_fresnel_ocean` in both C and Python.

## Surface taxonomy

| Surface | Meaning |
|---|---|
| `black` | Perfectly absorbing lower boundary; no Fresnel reflection |
| `flat` | Flat Fresnel interface over a black ocean |
| `black_fresnel_ocean` | Rough Fresnel interface over a black ocean; no water-leaving radiance |
| `ocean` | Coupled atmosphere-ocean RT with water-leaving radiance |

## Slope variance

The slope model is independent of the surface boundary name.

- Default `ocrt-floor`:
  `sigma^2 = 0.003 + 0.00512 * max(0.01, W)`
- Alternative `nakajima-tanaka`:
  `sigma = 0.0731 * sqrt(W)`

The default is CM54-derived but includes an explicit OCRT low-wind floor. It is
therefore not described as an unmodified Cox-Munk model.

## C changes

- Canonical enum: `RT_SURFACE_BLACK_FRESNEL_OCEAN`
- Canonical CLI: `--surface black_fresnel_ocean`
- Also accepts `black-fresnel-ocean`
- Legacy `--surface coxmunk` remains accepted with a deprecation warning
- New slope string selector:
  - `--sigma-model ocrt-floor`
  - `--sigma-model nakajima-tanaka`
- Legacy numeric `--sigma-type` remains supported
- Version:
  `OCRT-v1.2-2026-07-21-KST-black-fresnel-surface-naming`

## Python changes

- Canonical high-level functions:
  - `solve_atm_black_fresnel_ocean`
  - `solve_atm_aerosol_black_fresnel_ocean`
  - `solve_aerosol_black_fresnel_ocean_value`
- Canonical mode string: `black_fresnel_ocean`
- Legacy `coxmunk` strings and functions remain deprecated aliases
- Canonical low-level rough-Fresnel aliases added in `surface.py`
- Slope constants:
  - `SLOPE_MODEL_NAKAJIMA_TANAKA = 0`
  - `SLOPE_MODEL_OCRT_FLOOR = 1`
- Version:
  `pyOCRT-v1.2-2026-07-21-black-fresnel-surface-naming`

## Numerical regression

### C

Reduced Rayleigh + rough Fresnel test:

- SZA/VZA/RAA = 30/30/90 deg
- wavelength = 555 nm
- wind = 3 m/s
- gas absorption columns = 0
- n_mu = 8, n_layers = 20, m_max = 2, max_orders = 10
- direct sunglint decoupled

The following stdout files have the same SHA-256:

`c7b761a226b7325672fa40b97c0b769d800a4f4b208058fe95cad75fcb734316`

- original binary, `--surface coxmunk`
- new binary, `--surface black_fresnel_ocean`
- new binary, deprecated `--surface coxmunk`
- new binary, default slope model
- new binary, `--sigma-model ocrt-floor`
- new binary, `--sigma-type 1`

Thus the naming change produces byte-identical stdout.

The named and numeric Nakajima-Tanaka selections also match each other exactly:

`087d1bd7f4e4dd5e32573eaae1474de009e6d1981a58dd4d8927cab620309736`

### Python

Canonical and deprecated high-level functions returned identical tuples:

```text
canonical=(0.03880440998334026, 0.09354851181485427, [6, 6, 6])
legacy   =(0.03880440998334026, 0.09354851181485427, [6, 6, 6])
```

Slope values:

```text
OCRT floor, W=0: 0.0030512
OCRT floor, W=3: 0.01836
Nakajima-Tanaka, W=3: 0.01603083
```

## Regression tests

- C build: PASS
- C Mie P11/P12/P33 wavelength PCHIP: PASS
- C external bottom mode-bound and shape guards: PASS
- Python compileall: PASS
- Python full-grid aerosol-object tests: 3/3 PASS

## Numerical impact

None intended and none observed. This change separates physical naming from the
slope-model naming and preserves legacy inputs for migration.
