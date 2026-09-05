# OCRT public-output RT derivation audit — 2026-07-21

- Issue: `OCRT-OUTPUT-TRANS-RT-001`
- Priority: P0, scientific-output correctness
- Scope: public stdout, ocean full-grid CSV, atmospheric-correction LUT scripts
- Constraint: solver state, TOA reflectance and `Rrs/rrs` must not change

## 1. Confirmed defect

The previous public output used

```text
T_diff_up_hemi  = T_diff_dn_hemi
T_total_up_view = exp(-tau_atm/mu_view) + T_diff_up_hemi
```

This is a reciprocity/hemispheric approximation.  It is not the directional
transmission of the actual ocean water-leaving BRDF and therefore must not be
published as an OCRT result.

The current `scripts/ocrt_ac_lut.py` also reconstructed `Tup_total` by evaluating
the downward irradiance-transmittance table at `theta=vza`.  That assignment is
removed.  Black and black-Fresnel surface runs have no prescribed upward BOA
source, so their upward transmission is undefined rather than reciprocal.

## 2. Correct RT definitions

### Downward

The atmospheric SOS exports the BOA downward field.  OCRT integrates its `m=0`
diffuse radiance over the downward hemisphere:

```text
T_diff_dn_hemi = E_diff_dn(BOA) / E_solar_TOA
T_total_dn_hemi = T_dir_dn + T_diff_dn_hemi
```

`T_dir_dn` is the exact unscattered attenuation for the selected plane-parallel
or PSSA direct path.  `T_diff_dn_dir` is reconstructed from the solved Fourier
field at the requested direction.

### Upward

For an ocean-coupled solve, OCRT already propagates the actual water-leaving
bottom source through a second atmospheric SOS pass.  The exact effective
intensity transmission is therefore

```text
TOA_water_signal_I = I_TOA_total - I_atm_path
T_total_up_view = TOA_water_signal_I / Lu(0+)
T_diff_up_view = T_total_up_view - T_dir_up_view
```

This retains water BRDF, rough-surface transmission, multiple scattering and
view azimuth.  The value is source dependent; it is not a medium-only scalar.
When the rigorous bottom-source solve is absent or the denominator is zero,
OCRT returns `NaN` and `T_up_rt_valid=0` rather than a fallback approximation.

## 3. Full public-output audit

| Output or path | Derivation | Status after fix |
|---|---|---|
| `TOA_rho_I/Q/U` | Atmospheric and ocean SOS assembly | RT result; unchanged |
| `Rrs0plus`, `rrs0minus` | Solved boundary radiance divided by solved irradiance | RT result; unchanged |
| `Lu/Qu/Uu`, `Ed/Eu` at 0+/0- | Solved angular fields and quadrature | RT result; unchanged |
| `T_dir_dn` | Direct-beam attenuation; PSSA when enabled | Exact direct sub-solution |
| `T_diff_dn_hemi` | BOA diffuse radiance hemispheric integration | Exact RT diagnostic |
| `T_total_dn_hemi` | Direct plus solved diffuse | Exact RT diagnostic |
| `T_diff_dn_dir` | Solved BOA Fourier field at view | Exact RT diagnostic |
| old `T_diff_up_hemi=T_diff_dn_hemi` | Reciprocity approximation | Removed |
| `T_total_up_view` | `(I_TOA-I_atm_path)/Lu0plus` | Exact effective RT diagnostic |
| `T_diff_up_view` | Exact total minus exact direct | RT-derived diagnostic |
| old `T_sg_up_dir` | `rho_glint*mu_sun` | Kept only under explicit `diag_*` label |
| old `T_total_up_dir` | `rho_TOA*mu_sun` | Kept only under explicit `diag_*` label |
| `Kd0minus`, `Ku0minus` | One-layer finite difference of solved `Ed/Eu` | RT-field numerical diagnostic; not an analytic transport replacement |
| gas absorption output | Layer absorption is integrated in SOS in current v1.2 path | RT result; legacy standalone correction helpers are unused |
| `ocrt_ac_lut.py:Tup_total` for black surfaces | Former downward-table reciprocity | Removed; now NaN/invalid |
| `ocrt_ac_lut_ORIGINAL_pre3flag.py` | Historical reciprocity implementation | Marked historical/non-production |
| ocean full-grid with aerosol request | current driver calls the ocean solver with `aer=NULL` | Independent wiring defect: aerosol is not represented; not modified by this output-only patch |


## 4. Independent full-grid aerosol wiring issue

`run_ocean_rrs_full_grid_csv()` currently invokes `rt_solve_case_ocean(..., aer=NULL, ...)`.
Therefore an ocean `--output-full-grid` request cannot be used as an aerosol-coupled
transmittance product even when the base command supplied a Mie model and AOD.  This is
not an analytic-transmittance error and is outside the requested output-only patch,
because correcting it would change full-grid TOA/Rrs values.  It is recorded as a
separate high-priority wiring issue.  The exact transmittance columns added here are
valid for the atmosphere actually solved; the validation grid uses no aerosol.

## 5. Important independent fallback

The ocean assembly retains a legacy first-cut TOA fallback only when the
rigorous atmospheric bottom-source pass fails.  This patch does not alter that
path because the requested invariant is that TOA and `Rrs/rrs` remain
unchanged.  Crucially, the fallback is no longer exposed as a valid
transmittance: `T_up_rt_valid=0` and the upward fields are `NaN`.

A later solver-policy change may replace the fallback with fail-fast behavior,
but that is a separate core-output availability decision and must receive its
own regression review.

## 6. Acceptance requirements

1. All pre-existing physical outputs except transmission/diagnostic labels are
   numerically identical before and after this patch.
2. For valid ocean runs:
   `T_total_up_view * Lu0plus == TOA_water_signal_I` within floating precision.
3. `T_total_dn_hemi == T_dir_dn + T_diff_dn_hemi`.
4. Black-surface atmospheric runs output `NaN`, `T_up_rt_valid=0` for upward T.
5. No reciprocity-derived `Tup_total` remains in the production LUT script.
6. Full-grid ocean CSV includes exact downward and upward RT diagnostics.
