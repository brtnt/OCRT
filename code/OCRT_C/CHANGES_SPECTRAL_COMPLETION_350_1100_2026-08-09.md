# Changes — spectral completion review, 2026-08-09

- Explicit Chl absorption zero extension for 850 < wavelength <= 1100 nm.
- Source cutoff kept at 850 nm to avoid interpolation across the boundary.
- User-supplied pure-water Z09 table confirmed byte-identical and retained.
- Removed three unused legacy narrow Mie files.
- Added machine-readable spectral support manifest and current data inventory.
- Kept detritus 850-nm endpoint behavior unchanged; documented as provisional.
- Corrected psi_T header units to m^-1 degC^-1; values unchanged.
