# Python water internal-reflection Fourier-sign fix

Version: `pyOCRT-v1.2-2026-07-21-water-intrefl-msign-fix`

The water-side internal-reflection boundary now uses `msign = +1` for all
Fourier orders in the native single solver, NumPy/CuPy batch solver, and LUT
batch path.  Specular reflection preserves azimuth; the previous `(-1)^m`
incorrectly reversed odd-m Stokes-U contributions.
