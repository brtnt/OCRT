# pyOCRT full-grid aerosol object fix

Version: `pyOCRT-v1.2-2026-07-21-fullgrid-aerosol-object-fix`

## Modified source files

- `ocrt_py/aerosol.py`
- `ocrt_py/atmos.py`
- `ocrt_py/driver.py`
- `ocrt_py/lut.py`
- `ocrt_py/batch_driver.py`
- `ocrt_py/water.py`
- `ocrt_py/__init__.py`
- `ocrt_solve.py`
- `README.md`

## New files

- `ocrt_fullgrid.py`
- `tests/test_fullgrid_aerosol_object_fix.py`
- `validation/python_fullgrid_aerosol_object_fix/*`

## Object lifecycle

1. Read one `.mie` file.
2. Prepare one frozen `AerosolRuntime` for wavelength, AOD reference, interpolation, truncation and moment settings.
3. Pass the same object explicitly to every geometry cell.
4. Pass it unchanged to atmospheric pass 1 and pass 2 of every coupled solve.
5. Reject AOD-positive calls with no aerosol object/data.

The runtime arrays are NumPy read-only and the dataclass is frozen.
