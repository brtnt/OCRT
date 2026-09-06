# OCRT distribution, repository, and license policy

## Canonical repository

- Repository: https://github.com/brtnt/OCRT
- Repository identifier: `brtnt/OCRT`
- Default integration branch: `main`
- Do not use the former or similarly named `brtnt/OC-RT` repository.

## Distribution policy

1. Session deliverables are provided directly in the current session by default.
2. GitHub upload is performed only when the user explicitly requests Git backup or upload.
3. Every future OCRT Git operation must target `brtnt/OCRT`.
4. Preserve the latest validated source baseline and do not overwrite it with an older session package.
5. Large generated packages must be versioned and accompanied by SHA-256 manifests.
6. Every future OCRT migration or handoff document must repeat the canonical repository and this policy.

## License policy

OCRT is distributed under the custom `ACADEMIC AND NON-COMMERCIAL USE LICENSE` in the root `LICENSE` file.

- Permitted: academic, educational, and non-profit research use, copying, modification, and redistribution.
- Prohibited without a separate license: commercial use, private or for-profit company use, incorporation into commercial products, NDA-bound research, or work intended for commercial gain.
- Commercial licensing contact: `brtnt@kiost.ac.kr`.
- Redistribution must preserve the copyright notice, the complete `LICENSE` file, the canonical repository attribution, and the non-commercial restrictions.

## Current validated baseline

- Release identifier: `OCRT-v1.2-2026-07-21-KST-iop-kd-fullgrid-csv-fix`
- Base: `OCRT-v1.2-2026-07-20-KST-water-mmax-mode-bound-fix`
- Integrated corrections:
  - water Fourier external-bottom-source mode bound and shape contract;
  - PCHIP interpolation of Mie P11/P12/P33 across wavelength;
  - RT-derived upward atmospheric transmittance output and explicit invalid-state handling.
- Integrated aerosol inventory:
  - SnF canonical models: 16;
  - Ahmad/AccuRT models: 80;
  - Ahmad et al. (2010) paper-centric models: 80, with `A2010ver_` prefix.
