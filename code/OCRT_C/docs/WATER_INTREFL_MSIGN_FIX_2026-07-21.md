# Water internal-reflection Fourier-sign fix

Version: `OCRT-v1.2-2026-07-21-KST-water-intrefl-msign-fix`

Water-side specular reflection preserves physical azimuth.  The internal-
reflection top boundary therefore uses `msign = +1` for every Fourier order.
The previous `(-1)^m` factor corrupted odd-m contributions and depressed
off-nadir Stokes U/DoLP.

Changed C paths in `src/rt_solver.c`:
- `add_flat_intrefl_downfield_order_complete`
- `rt_solver_sos_pol_intrefl`
- `rt_solver_sos_pol_intrefl_rough`

`rt_fourier_pi_shift_sign()` remains available for other conventions and was
not modified.
