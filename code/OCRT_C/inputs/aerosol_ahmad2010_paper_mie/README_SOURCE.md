# A2010ver Ahmad et al. (2010) 80-model Mie family

This download contains 80 OCRT-compatible paper-centric Ahmad et al. (2010)
aerosol Mie caches. Every aerosol model filename has the requested `A2010ver_`
prefix.

Examples:

- `A2010ver_r30f95v01.mie`: RH 30%, fine-mode volume fraction 0.95
- `A2010ver_r80f20v01.mie`: RH 80%, fine-mode volume fraction 0.20
- `A2010ver_r95f00v01.mie`: RH 95%, coarse-only endpoint

Numerical format:

- 20 wavelengths, 0.350–3.750 µm
- 361 scattering angles, 180° to 0° at 0.5° spacing
- P11, P12, and P33 phase-matrix blocks

Scientific boundary:

- 0.412–0.865 µm is constrained by Ahmad et al. (2010) Tables 3 and 4.
- Outside that interval, the package uses the documented constituent-level
  engineering extension from the prior OCRT paper-centric reconstruction.
- This family is intentionally distinct from the simplified AccuRT effective-mode family.
