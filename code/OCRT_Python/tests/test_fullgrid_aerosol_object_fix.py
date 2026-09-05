#!/usr/bin/env python3
import os
import sys
import types
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ocrt_py.aerosol import read_mie, prepare_aerosol_runtime
from ocrt_py.lut import run_ocean_aerosol_full_grid


class FullGridAerosolObjectFixTests(unittest.TestCase):
    def test_runtime_is_frozen_and_arrays_are_read_only(self):
        mie = read_mie(ROOT / 'data' / 'C50.mie')
        rt = prepare_aerosol_runtime(mie, 490.0, 0.2,
                                     aod_ref_nm=865.0, L_max=16)
        self.assertAlmostEqual(rt.aod_target, 0.2719174901355678, places=14)
        self.assertAlmostEqual(rt.tau_a_eff, 0.24393084023659078, places=14)
        self.assertAlmostEqual(rt.ssa_a_eff, 0.95370972550096, places=14)
        self.assertFalse(rt.P11_norm_asc.flags.writeable)
        self.assertFalse(rt.betal.flags.writeable)
        with self.assertRaises(ValueError):
            rt.P11_norm_asc[0] = 0.0
        with self.assertRaises(Exception):
            rt.aod_target = 0.0

    def test_aod_positive_requires_explicit_aerosol_input(self):
        with self.assertRaisesRegex(ValueError, 'Aerosol runtime is required'):
            run_ocean_aerosol_full_grid(
                sza_deg=30.0, wavelength_nm=490.0, n_water=1.34,
                wind_speed=3.0, a_tot=0.2, b_tot=0.3, bb_tot=0.01,
                phase_lut_path=None, mie=None, aod_ref=0.2,
                aod_ref_nm=865.0, vza_values=[0.0], raa_values=[0.0])

    def test_prepare_once_and_native_solve_once_for_100_cells(self):
        import numpy as np
        calls = {'prepare': 0, 'solve': 0, 'ids': []}
        fake = types.SimpleNamespace(
            aod_ref=0.2, aod_ref_nm=865.0, aod_target=0.27,
            tau_a_eff=0.24, ssa_a_eff=0.95)

        def fake_prepare(*args, **kwargs):
            calls['prepare'] += 1
            return fake

        def fake_solve(*args, **kwargs):
            calls['solve'] += 1
            calls['ids'].append(id(kwargs['aerosol_runtime']))
            nv, nr = len(args[1]), len(args[2])
            shp = (nv, nr)
            return dict(
                Ed_0plus=1.0, Lu_0plus=np.full(shp, 0.1),
                Qu_0plus=np.full(shp, 0.01), Uu_0plus=np.full(shp, 0.001),
                Rrs0plus=np.full(shp, 0.1), Rrs0plus_Q=np.full(shp, 0.01),
                Rrs0plus_U=np.full(shp, 0.001), Ed_0minus=0.9,
                Lu_0minus=np.full(shp, 0.08), Qu_0minus=np.full(shp, 0.008),
                Uu_0minus=np.full(shp, 0.0008), rrs0minus=np.full(shp, 0.08/0.9),
                rrs0minus_Q=np.full(shp, 0.008/0.9),
                rrs0minus_U=np.full(shp, 0.0008/0.9),
                rho_I=np.full(shp, 0.2), rho_Q=np.full(shp, 0.02),
                rho_U=np.full(shp, 0.002), orders=7, conv=1,
                diagnostics={'atmos_pass1_solves': 1, 'water_solves': 1,
                             'atmos_pass2_solves': 1, 'cell_solver_calls': 0,
                             'vza_count': nv, 'raa_count': nr})

        with patch('ocrt_py.aerosol.prepare_aerosol_runtime', fake_prepare),              patch('ocrt_py.driver.solve_case_ocean_coupled_lut', fake_solve):
            out = run_ocean_aerosol_full_grid(
                sza_deg=30.0, wavelength_nm=490.0, n_water=1.34,
                wind_speed=3.0, a_tot=0.2, b_tot=0.3, bb_tot=0.01,
                phase_lut_path=None, mie=object(), aod_ref=0.2,
                aod_ref_nm=865.0,
                vza_values=list(range(10)), raa_values=list(range(10)))
        self.assertEqual(calls['prepare'], 1)
        self.assertEqual(calls['solve'], 1)
        self.assertEqual(set(calls['ids']), {id(fake)})
        self.assertEqual(out['diagnostics']['aerosol_prepare_count'], 1)
        self.assertEqual(out['diagnostics']['cell_solver_calls'], 0)
        self.assertEqual(len(out['rows']), 100)


if __name__ == '__main__':
    unittest.main()
