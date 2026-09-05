# Validation reference policy after the 2026-07-28 pure-water IOP replacement

Historical output files that depend on pure-water absorption or scattering are not current numeric
anchors after this intentional data replacement. Current gates are the Python test suite, especially
`test_pure_water_z09_table.py`, `test_native_coupled_lut.py`, `test_value_kernel_native_lut.py`, and
`test_water_intrefl_msign.py`. The full suite passes with the corrected table.
