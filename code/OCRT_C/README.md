# Current canonical note (2026-08-26)

The current performance-patched C canonical tree is platform-neutral. Native-Windows-only gates and launch workflows are withdrawn. See `README_FINAL_PLATFORM_NEUTRAL_2026-08-26_KO.md` and `docs/NATIVE_WINDOWS_ROLLBACK_PLATFORM_NEUTRAL_FINAL_2026-08-26_KO.md`.

# OCRT

Ocean-color coupled vector radiative-transfer software.

## Current validated release

`OCRT_C_FINAL_PLATFORM_NEUTRAL_PERF_2026-08-26`

Base lineage: `OCRT-v1.2-2026-08-16-KST-mie-fr631-direct-truncation`

### Mie FR631 direct-particle release (2026-08-16)

- The active aerosol catalog (176 models), EAP catalog (17 models), and AHN TSM catalog (4 models) use the fixed 631-node Forward-Robust scattering-angle grid where source data are available or reproducibly remapped.
- Particle P11/P12/P33 are evaluated with one common piecewise-linear weight in scattering angle. Cubic angle interpolation and particle L=200 reconstruction are not used by the production direct kernel.
- The direct vector Fourier kernel is cached in worker-private RAM. One case performs one native all-view solve and reconstructs every requested VZA×RAA cell; production per-cell replay is prohibited.
- The OSOAA-style hydrosol forward-cap option uses `mu1=0.85`, `mu2=0.92`, and `A_TRONCA>=0.1`. The residual phase remains a direct theta-linear LUT and transport scattering is scaled as `b_eff=b*(1-A/2)`.
- When the truncation criterion is not met, the truncation operation is an exact no-op: phase arrays, metadata, scalar IOPs, and RT output remain unchanged.
- A pressure-zero, no-atmosphere case has a native one-water-solve all-view path; it does not fall back to VZA×RAA cell replay.
- Independent batch rows reset case-local atmospheric all-view state and retain no mutable dependency on preceding or following rows.

This source includes:

- worker-private caching of complete parsed Mie models and complete direct-vector Fourier kernels;
- mandatory native all-view execution: one case produces every requested VZA×RAA cell without per-cell solver replay;
- independent batch-row state, including explicit reset of atmospheric S7/S7b all-view caches at every case boundary;
- fixed-grid θ-linear direct-particle kernel reuse with cached Fourier trigonometric tables;
- exact SOS zero-quadrature-column skipping with a byte-regression disable switch;
- RAA-invariant Beer-factor hoisting in the native coupled angular LUT;
- stage-2 water-output RAA convention correction: reported `Rrs(0+)` and `rrs(0-)` now use the same public OCRT RAA and positive-sine Stokes-U reconstruction as the atmospheric output;
- explicit separation of public output reconstruction from the local single-scatter `pi-RAA` propagation-vector geometry;
- source-level and runtime regression tests for the water RAA convention, full-grid/single-geometry consistency, and the wind-zero branch;
- canonical `--surface black_fresnel_ocean` naming for a rough Fresnel interface over a black ocean;
- deprecated `--surface coxmunk` compatibility alias;
- explicit `--sigma-model ocrt-floor|nakajima-tanaka` naming, with the default OCRT low-wind floor law documented separately from the surface boundary condition;

- water Fourier external-bottom-source mode-bound and shape fix;
- shape-preserving PCHIP wavelength interpolation for Mie P11/P12/P33;
- RT-derived upward transmittance output with explicit validity flag;
- physically corrected Kd(0-) using unscattered transmitted skylight at
  the surface and first water level;
- full-grid CSV constituent IOP output, including phytoplankton-only `a_chl`;
- 16 SnF, 80 Ahmad/AccuRT, and 80 A2010ver aerosol Mie models in the full package.

### Stage-2 air-to-water interface correction (2026-07-24)

- Uses signed photon-propagation cosines and the corresponding pi-shifted relative azimuth in the rough air-to-water Stokes rotation.
- Stores the m>0 incoming-U Fourier-column sign in the TAW operator so the downstream contraction remains an ordinary 3x3 matrix-vector product.
- Adds independent polar-limit and Fourier-column regression tests.
- Closes the diffuse-top primary-source incoming-U column against the production SOS operator for every mode/node/basis direction tested. The three molecular U-input coefficients in `ocrt_add_diffuse_top_primary` carry the sign required by the existing U-output wrapper.
- Adds a reference-free source/operator contraction regression for pure-Rayleigh and mixed Rayleigh–aerosol media.
- The correction materially reduces high-SZA Q/U residuals, while exposing separate low-SZA/intensity residuals that remain under investigation; it is not an empirical fit to OSOAA.

### Exact-pole surface rotation and reproducible build (2026-07-24)

- Applies the analytical meridian-rotation limit when either incident or outgoing direction is exactly vertical in all four rough-interface Mueller kernels (`R_air`, `T_wa`, `R_ww`, and `T_aw`).
- Preserves the already validated non-pole equations, FIX1–FIX4 conventions, and atmospheric harness.
- Adds exact-pole/one-sided-limit continuity and vertical spin-2 covariance regressions.
- Documents that the simultaneous double-pole ambiguity is multiplied by zero quadrature weight in production boundary sums; a deliberately perturbed double-pole rotation leaves the audited full-grid outputs unchanged.
- Replaces host-dependent `-march=native` release builds with a fixed default `-march=cascadelake`. Set `OCRT_MARCH=<target>` only when a different documented CPU baseline is required; changing it intentionally changes the binary fingerprint and may change last-bit floating-point results.

### Water-to-air coupling interpolation clamp correction (2026-07-25)

- Removes the silent 64-direction cap in the water-to-air coupling interpolation path.
- Preserves appended zero-weight view, solar, and nadir directions used by the water solver; these directions must participate in the sorted interpolation grid.
- Uses a 320-entry stack workspace for the validated solver range and a fail-loud heap fallback for future larger grids.
- Restores the exact nadir spin-2 constraint at `n_mu_water=64` and the LUT default `n_mu_water=96`; `rrs(0-)` is unchanged because the defect is downstream of the in-water solve.
- Adds regressions for `n_mu_water=48/64/96`, forbidden nadir `m=0` polarization, full-grid/single-geometry consistency, and the no-trigger byte-identical path.
- Does not modify the atmospheric solver, atmospheric harness, FIX1–FIX4 interface kernels, exact-pole equations, or in-water SOS physics.

## Build

```bash
./scripts/build_release_v1.2.sh build/ocrt
```

The default build target is `cascadelake` for deterministic code generation on the validated toolchain. Reproducibility claims require the same source, compiler/linker versions, flags, environment, and `OCRT_MARCH`. To use another explicit target:

```bash
OCRT_MARCH=x86-64-v3 ./scripts/build_release_v1.2.sh build/ocrt
```

For bit-level CSV comparisons, use raw byte comparison or `scripts/compare_csv_bitexact.py`. Do not use values already parsed through `pandas` as the source of a bit-exact assertion.

## Verification

Run from the full package root:

```bash
./verify_package.sh
./verify_package.sh --full-rebuild
```

## Repository and distribution

Canonical repository: https://github.com/brtnt/OCRT

Session deliverables are provided directly in the session by default. GitHub upload is performed only when explicitly requested by the user.

## License

Academic and non-commercial use only. See `LICENSE`. Commercial licensing inquiries: brtnt@kiost.ac.kr.

## EAP 17-species phase generator status

The 17-species coated-sphere EAP generator and FR631 P11/P12/P33 files are included and pass the active-file physical-contract gate. Species-specific phytoplankton scattering remains intentionally disabled in the constituent production model pending the final end-to-end scientific acceptance gate for the EAP closure. It is not silently selected by chlorophyll input.

- Chl without `--ocrt-phyto-group`: the frozen validated phytoplankton absorption table is used; phytoplankton `b=bb=0`; detritus scattering remains.
- Chl with explicit `--ocrt-phyto-group`: fail-loud, exit code 2.
- The 17 EAP files remain available for generator validation and explicit fixed-bulk/direct-kernel diagnostics.
- The direct-kernel and forward-truncation implementation is validated independently; enabling a species in the constituent closure requires a separate scientific-model decision.

## Direct particle kernel, truncation, and spectral scope

The production particle-phase path is the direct theta-linear vector kernel. The legacy moment representation is retained only for controlled diagnostics:

- default/direct: `OCRT_WATER_PARTICLE_KERNEL=direct`
- diagnostic legacy comparison: `OCRT_WATER_PARTICLE_KERNEL=moment`
- azimuth quadrature used while building a cold direct kernel: `OCRT_WATER_VALUE_NPHI=N`

Once a complete Fourier kernel is built, it is cached in worker-private RAM and reused for the case. The requested output RAA grid does not trigger repeated SOS solves.

The active AHN TSM files were remapped to FR631 while preserving every historical 0.5-degree node at common phase wavelengths. The temporary 330-nm TSM phase contract is a 350-nm endpoint hold because a canonical AHN microphysical generator was not present in the supplied package. Source-informed regenerated candidates were retained under validation and were not installed.

The active legacy detritus file remains on its original variable angle grid and is evaluated by the same theta-linear direct evaluator. The rejected 330–1100-nm detritus candidate remains excluded from production.

At exact VZA=0, the polarized direct-kernel output extraction retains the documented singularity. Use VZA=0.001 degree for polarized diagnostic output at nominal nadir; no silent geometry substitution is applied.

The OCRT constituent spectral contract is 330–1100 nm subject to each selected component's explicit data contract. Pure-water optical properties are tabulated over the wider Z09 interval; the separate temperature-correction table must also cover the requested wavelength when temperature differs from 20 degrees C.
