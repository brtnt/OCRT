# OCRT v1.2 clean-room release notes

## Scope

OCRT v1.2 replaces the production atmospheric numerical kernels that had strong
external-code lineage markers with independently structured OCRT modules while
preserving the OCRT numerical contract. The Mie converter and Python utilities
are outside this clean-room scope.

## Main changes

- US Standard Atmosphere 1962 cumulative molecular profile for Rayleigh layering.
- Fixed 2 km exponential aerosol vertical profile; legacy ODA550/an23 dispatch removed.
- Independent angle-grid, generalized angular-basis, vector phase-moment,
  Fourier-kernel, primary-source, layer-transport, SOS-operator, molecular-phase,
  and surface-boundary modules.
- Removed production 6SV Rayleigh compatibility mode, OSOAA phase-moment input,
  AF1982/OSOAA experimental surface switches, and obsolete parity probes.
- Fixed a sanitizer-detected lifetime defect in the atmosphere-to-water diffuse
  boundary quadrature (`stack-use-after-scope`) without changing numerical output.

## Numerical verification

- 47 canonical conditions: byte-exact against the pre-fix clean-room candidate.
- 63 expanded surface/geometry conditions: byte-exact.
- Water branch contract: 14/14 PASS.
- CCRR Chl: 11/11 PASS.
- OCRT organic Chl: 10/10 PASS.
- Ahn TSM: 6/6 PASS.
- TSM phase-cache: 168 array comparisons PASS; build trace unchanged.
- ASan/UBSan representative OCRT/CCRR/IOP runs: 0 findings after lifetime fix.
- Representative performance: median runtime unchanged (0.23 s vs 0.23 s);
  peak RSS difference below 0.1 MB in the measured case.

## Important interpretation

US62 changes the physical altitude allocation of Rayleigh optical depth. Pure
Rayleigh optical-depth-coordinate solutions remain byte-identical. Rayleigh +
aerosol or Rayleigh + absorbing-gas cases may change slightly because their
vertical overlap is now physically different.

References to OSOAA, 6SV, and AF1982 retained in comments or validation documents
identify comparison targets and historical validation, not linked source modules.
A legal clean-room opinion is outside the scope of this technical release.
