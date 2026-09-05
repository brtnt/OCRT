# CHANGELOG — RT-derived transmittance output fix (2026-07-21)

## Fixed

- Removed the invalid public approximation `T_diff_up = T_diff_dn_hemi`.
- Removed the derived `T_total_up_view = Beer_up + T_diff_dn_hemi` output.
- Added rigorous ocean-source-derived `T_total_up_view`, `T_diff_up_view`, and `TOA_water_signal_I`.
- Added `T_up_rt_valid` so undefined upward transmission is explicit.
- Added `T_total_dn_hemi` as the exact sum of direct and diffuse downward components.
- Removed `Tup(vza)=Tdn(vza)` reconstruction from `scripts/ocrt_ac_lut.py`.

## Also integrated

- Shape-preserving PCHIP wavelength interpolation for Mie `P11/P12/P33`.
- Compatibility update to the PSSA numerical smoke test for intentional invalid-output NaNs.

## Preserved

- Water `m_max` external-bottom-source bounds and shape guard.
- TOA I/Q/U, Rrs/rrs, Lu/Ed/Eu, IOP, order, and convergence values under the output-only correction, relative to the same PCHIP baseline.
