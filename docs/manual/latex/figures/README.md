# Editable manual figures

All twelve figures are native LaTeX/TikZ or PGFPlots. They rebuild with the manual,
remain sharp when zoomed, and use `\FigLang{Korean}{English}` for shared labels.

| Files | Topic |
|---|---|
| `geometry_zenith_refraction.tex` | Pixel-based SZA/VZA and water refraction |
| `geometry_curvature.tex` | Pixel VZA versus satellite off-nadir |
| `geometry_glint_branches.tex` | RAA=0/180 branches, photon and look directions |
| `geometry_azimuth.tex` | SAA/VAA, full-circle RAA and signed-U mirrors |
| `options_surfaces.tex`, `options_glint.tex` | Surface modes and direct-glint separation |
| `options_gas.tex`, `options_aod.tex` | Pressure, gas profiles/columns and AOD reference wavelength |
| `options_water.tex`, `options_grids.tex` | Water inputs/IOPs and output grid versus numerical resolution |
| `measured_azimuth.tex` | Calculated azimuth curves from the CSV files in `data/` |
| `output_levels.tex` | TOA, above/below water, and output normalization |

The geometry and option panels are labelled schematics, not retrieved numerical
fields or full ray tracers. The measured-azimuth figure is explicitly identified
as real C output; its generator, metadata, original grids and selected plotting
data are included. See [data/README.md](data/README.md) before updating those data.
