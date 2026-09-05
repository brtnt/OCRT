> **역사 문서:** 이 문서는 v1.15 통합 당시 옵션명을 기록한다. OCRT v1.18의
> canonical 입력은 `--water-model ocrt --ocrt-chl ... --ocrt-tsm ...
> --ocrt-adom440 ... [--ocrt-tsm-species ...]`이다. 최신 계약은
> `WATER_INPUT_INTERFACE_v1.18_2026-07-18.md`를 따른다.

# OCRT v1.15 — Ahn four-species TSM integration

## Scope

This change integrates dry-weight TSM concentration into the existing OCRT
water radiative-transfer path. It does not replace the OCRT RT solver and does
not implement the broader water-optics interface redesign.

## Public input

```text
--simple-min C
--ccrr-min-g-m3 C          # deprecated alias
--tsm-species NAME
```

`C` is dry-weight TSM concentration in `g m^-3`. Supported species are:

| CLI name | Phase/optical file | Default |
|---|---|---:|
| `red_clay` | `Red_clay_AHN.mie` | yes |
| `brown_earth` | `Brown_earth_AHN.mie` | no |
| `yellow_clay` | `Yellow_clay_AHN.mie` | no |
| `calcareous_sand` | `Calcareous_sand_AHN.mie` | no |

The user does not supply `--water-mie-phase`. The selected TSM species maps to
its matching vector `.mie` automatically. Supplying both is an error.

## Optical conversion

For each wavelength:

```text
a_min(lambda)  = a*(lambda) C
b_min(lambda)  = b*(lambda) C
bb_min(lambda) = b_min(lambda) [bb/b](lambda)
```

The delivered `.mie` bulk spectral block is used as follows:

```text
a*(lambda) = Extinct_Co(lambda) - Scatter_Co(lambda)
b*(lambda) = Scatter_Co(lambda)
```

The audit tables in `inputs/tsm_ahn/astarmin_*` and `bstarmin_*` contain the
same 757 values without resampling. If those tables are absent, the adapter
recovers the values directly from the `.mie` file. `[bb/b](lambda)` is
integrated from the matching P11 phase and interpolated in wavelength.

## OCRT coupling

The TSM coefficients are added to the existing pure-water, CDOM and pigment
components. The selected mineral `.mie` is passed to the existing OCRT vector
water-phase moment path. It supplies P11/P12/P33 and preserves the existing
Mie-moment and forward-peak processing of the OCRT solver.

When both Chl and TSM are nonzero, Chl-derived absorption/scattering
coefficients remain in the bulk IOP. In the current v1.14 constituent branch,
the selected TSM `.mie` is the particulate vector phase source. A separately
weighted organic-particle/mineral vector-phase mixture is intentionally not
introduced in this narrow integration; it belongs to the pending full
water-optics interface design.

## Data replacement

The former package files below were deleted:

```text
aux/phase_mie/Red_clay.mie
aux/phase_mie/Brown_earth.mie
aux/phase_mie/Yellow_clay.mie
aux/phase_mie/Calcareous_sand.mie
```

The delivered replacements are canonical under `inputs/tsm_ahn/`.

## 443 nm verification at C = 2 g m^-3

| Species | a_min [m^-1] | b_min [m^-1] | bb/b |
|---|---:|---:|---:|
| red_clay | 0.15594 | 1.72062 | 0.0224741733515 |
| brown_earth | 0.20748 | 1.56368 | 0.0201694512447 |
| yellow_clay | 0.07834 | 1.77760 | 0.0201259908891 |
| calcareous_sand | 0.05120 | 2.07360 | 0.0228875249566 |

## Runtime data resolution

OCRT resolves `inputs/tsm_ahn` from normal package working directories and
from the executable location. `OCRT_TSM_DIR` is the authoritative override.
A missing or malformed override is a fatal input error.
