# Data for the calculated azimuth figure

These are actual OCRT C results, generated with `examples/generate_figure_data.py`.
They are typesetting inputs; ordinary PDF builds do not rerun the solver.

- `glint_w1.csv`, `glint_w5.csv`, `decoupled_w5.csv`: native 144-row, 13-column grids.
- `rayleigh_azimuth.csv`: VZA=40 subset, combined for plotting; the final 360-degree row repeats zero only to close the plotted curves.
- `rayleigh_azimuth_metadata.json`: source/executable provenance, full argument templates, and numerical checks.

Conditions: SZA=40 degrees, wavelength=555 nm, pressure=1013.25 hPa, no aerosol,
all six gas columns zero, IPSS off. The retained full angular grid has VZA=0/40
degrees and RAA=0..355 degrees at 5-degree intervals. Wind is 1 or 5 m/s;
the decoupled case uses wind=5 m/s. Numerical settings are the current C defaults.
All three runs exit successfully. The two curves including direct glint peak
at RAA=180 degrees. The decoupled mirror residuals are zero at printed precision
for I/Q and about 1.84e-17 for signed U.

To update data, run the generator from the full repository root with a fresh
output directory, inspect its logs and metadata, then replace the matching data
files here. Do not alter individual curve points by hand. This is a worked
configuration, not a production convergence certificate.
