#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""OCRT Python atmosphere-ocean angular full-grid CLI.

The aerosol .mie file and geometry-independent runtime are prepared exactly
once, then the same immutable object is passed to every VZA/RAA cell.
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
DATA = os.path.join(HERE, 'data')
from ocrt_py.constituent import PHYTO_GROUP_CHOICES, validate_phyto_group_option


def _float_list(text):
    return [float(x.strip()) for x in text.split(',') if x.strip()]


def main():
    ap = argparse.ArgumentParser(description='pyOCRT ocean aerosol full-grid')
    ap.add_argument('--out', required=True)
    ap.add_argument('--wl', type=float, required=True)
    ap.add_argument('--sza', type=float, required=True)
    ap.add_argument('--wind', type=float, required=True)
    ap.add_argument('--chl', type=float, required=True)
    ap.add_argument('--tsm', type=float, required=True)
    ap.add_argument('--acdom440', type=float, required=True)
    ap.add_argument('--aod865', type=float, required=True)
    ap.add_argument('--aer', default='C50')
    ap.add_argument('--vza-values', default='0,30')
    ap.add_argument('--raa-values', default='0,90,180,270')
    ap.add_argument('--phyto-group', default=None, choices=PHYTO_GROUP_CHOICES,
                    help='omit for absorption-only Chl; explicit species with Chl>0 is rejected in stage-2')
    ap.add_argument('--tsm-species', default='red_clay',
                    choices=['red_clay', 'brown_earth', 'yellow_clay', 'calcareous_sand'])
    ap.add_argument('--cdom-slope', type=float, default=0.014)
    ap.add_argument('--n-mu-water', type=int, default=48)
    ap.add_argument('--m-max', type=int, default=16)
    ap.add_argument('--n-layers-atm', type=int, default=400)
    ap.add_argument('--water-m-max', type=int, default=30)
    ap.add_argument('--water-n-layers', type=int, default=300)
    ap.add_argument('--max-it-water', type=int, default=500)
    ap.add_argument('--pressure', type=float, default=1013.25)
    ap.add_argument('--no-gas', action='store_true')
    a = ap.parse_args()
    explicit_phyto = a.phyto_group is not None
    a.phyto_group = validate_phyto_group_option(
        a.phyto_group, chl=a.chl, explicit=explicit_phyto)

    from ocrt_py.aerosol import read_mie
    from ocrt_py.constituent import OCRTConstituentModel
    from ocrt_py.absorption import Absorption
    from ocrt_py.lut import run_ocean_aerosol_full_grid
    from ocrt_solve import _constituent_greek

    mie_name = a.aer if a.aer.endswith('.mie') else a.aer + '.mie'
    mie = read_mie(os.path.join(DATA, mie_name))
    cm = OCRTConstituentModel(DATA, a.phyto_group, a.tsm_species)
    rr, constituent = _constituent_greek(
        cm, a.wl, a.chl, a.tsm, a.acdom440, cdom_slope=a.cdom_slope)
    absorption = None if a.no_gas else Absorption(DATA + '/afgl_atm', DATA + '/xsec')

    out = run_ocean_aerosol_full_grid(
        sza_deg=a.sza, wavelength_nm=a.wl, n_water=1.34,
        wind_speed=a.wind, a_tot=rr['a'], b_tot=rr['b'], bb_tot=rr['bb'],
        phase_lut_path=None, mie=mie, aod_ref=a.aod865, aod_ref_nm=865.0,
        vza_values=_float_list(a.vza_values), raa_values=_float_list(a.raa_values),
        pressure_hpa=a.pressure, absorption=absorption,
        n_mu_water=a.n_mu_water, fourier_m_max=a.m_max,
        nt_atm=a.n_layers_atm, m_max_water=a.water_m_max,
        n_layers_water=a.water_n_layers, max_it_water=a.max_it_water,
        constituent=constituent, out_csv=a.out, progress=True)
    rt = out['aerosol_runtime']
    print('written=%s cells=%d AOD_band=%.10g tau_eff=%.10g ssa_eff=%.10g' %
          (a.out, len(out['rows']), rt.aod_target, rt.tau_a_eff, rt.ssa_a_eff))
    print('lifecycle=%s' % out['diagnostics'])


if __name__ == '__main__':
    main()
