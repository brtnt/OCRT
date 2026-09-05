# OCRT value-kernel convention, polarization, and NIR extension

Version: `OCRT-v1.2-2026-07-27-KST-value-kernel-pol-nir`

- `OCRT_VALUE_PHASE_SPLINE=1`: scalar water value kernel uses cubic spline in `mu=cos(theta)` plus Gauss normalization (`--water-mie-moment-n-mu`, default 400).
- `OCRT_WATER_VALUE_KERNEL_POL=1`: opt-in six-component polarized direct value kernel (`pfm,gr,gt,arr,art,att`). Direct tables use a dedicated cache and never overwrite the coefficient/moment cache.
- Default gate-OFF path is unchanged.
- OCRT Chl wavelength gate extends to 1100 nm. Above the actual 850-nm absorption-table limit, Chl absorption is zero with a one-time warning. Pure-water data remain authoritative only through 900 nm.
- Exact VZA=0 remains a known polarized value-kernel output-extraction singularity. No internal numerical fix or silent input substitution is applied. Use VZA=0.001 deg (safe range 0.001--0.01 deg; avoid <1e-6 deg).
