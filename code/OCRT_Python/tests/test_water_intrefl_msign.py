#!/usr/bin/env python3
import inspect
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ocrt_py import sos, atmos_batch, lutbatch


class WaterInternalReflectionMsignTests(unittest.TestCase):
    def test_all_python_water_paths_use_azimuth_preserving_sign(self):
        for fn in (sos.sos_pol_intrefl_rough,
                   atmos_batch.sos_water_intrefl_batch,
                   lutbatch.solve_batch):
            src = inspect.getsource(fn)
            self.assertIn('msign = 1.0', src)
            self.assertNotIn('m % 2 == 0 else -1.0', src)

    def test_batch_solver_odd_mode_boundary_is_positive(self):
        atm = SimpleNamespace(B=1, nt=1, n_mu=1,
                              mu_pos=np.array([[0.5]]),
                              gb=np.array([[2.0/np.pi, np.nan, 2.0/np.pi]]))
        shape=(1,2,3)
        pI=np.zeros(shape); pQ=np.zeros(shape); pU=np.zeros(shape)
        pI[0,0,2]=1.; pQ[0,0,2]=2.; pU[0,0,2]=3.
        K=np.zeros((1,1,1,9)); K[0,0,0,0]=K[0,0,0,4]=K[0,0,0,8]=1.
        z=np.zeros(shape)
        def integrate(_atm,_jI,_jQ,_jU,tI,tQ,tU,top=False):
            oI=np.zeros(shape); oQ=np.zeros(shape); oU=np.zeros(shape)
            oI[0,0,0]=tI[0,0]; oQ[0,0,0]=tQ[0,0]; oU[0,0,0]=tU[0,0]
            return oI,oQ,oU
        with patch.object(atmos_batch,'sos_build_source_pol_batch',return_value=(z,z,z)), \
             patch.object(atmos_batch,'integrate_bcs_iqu',side_effect=integrate):
            I,Q,U,*_=atmos_batch.sos_water_intrefl_batch(
                atm,1,None,pI,pQ,pU,None,None,None,K,2,1e9)
        self.assertGreater(float(I[0,0,0]),0.0)
        self.assertGreater(float(Q[0,0,0]),0.0)
        self.assertGreater(float(U[0,0,0]),0.0)

    def test_single_solver_odd_mode_boundary_is_positive(self):
        # One layer, one positive stream.  Choose w so C_m*mu*w=1 for m=1.
        atm = SimpleNamespace(n_layers=1, n_mu=1,
                              rm=np.array([-0.5, np.nan, 0.5]),
                              gb=np.array([1.0/np.pi/0.5, np.nan, 1.0/np.pi/0.5]))
        shape=(2,3)
        pI=np.zeros(shape); pQ=np.zeros(shape); pU=np.zeros(shape)
        pI[0,2]=1.; pQ[0,2]=2.; pU[0,2]=3.
        K=np.zeros((1,1,9)); K[0,0,0]=K[0,0,4]=K[0,0,8]=1.
        z=np.zeros(shape)
        def integrate(_atm,_src,_bot=None,top=None):
            out=np.zeros(shape)
            if top is not None: out[0,0]=top[0]
            return out
        with patch.object(sos,'sos_build_source_pol',return_value=(z,z,z)), \
             patch.object(sos,'integrate_bcs',side_effect=integrate):
            I,Q,U,*_=sos.sos_pol_intrefl_rough(atm,1,None,pI,pQ,pU,K,2,1e9)
        self.assertGreater(I[0,0],0.0)
        self.assertGreater(Q[0,0],0.0)
        self.assertGreater(U[0,0],0.0)


if __name__ == '__main__':
    unittest.main()
