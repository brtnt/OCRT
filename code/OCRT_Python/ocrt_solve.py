#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""OCRT 파이썬 네이티브 CLI — GOCI-III 편광 논문 반사도 산출.

C --water-model ocrt 경로를 파이썬으로 재현한다(rho_I 5.6e-12 검증).
행/케이스마다 세 반사도를 차감 분리한다:
  R1 rho_TOA        = ocean + 대기(레일리+기체) + 에어로졸 + 구성모델 water-leaving
  R2 rho_R          = black Fresnel ocean surface 레일리 대기(에어로졸 없음), 선글린트 decouple
  R3 rho_(R+A,black)= black Fresnel ocean surface 레일리+에어로졸(black ocean), 선글린트 decouple
  도출: rho_RC = R1-R2, rho_A+rho_RA = R3-R2, t*rho_w = R1-R3, Rrs = R1의 Rrs(0+)

주의: 파이썬 결합 1행은 축소조건에서 ~70초, 전체 노브(n_mu_water=48/nt=400/m_max=16)는
수 분/행이다.  대량 자료생산은 C 바이너리(bin/ocrt) 권장.  이 CLI는 검증·소규모용이다.

사용:
  단일:  python3 ocrt_solve.py single --wl 555 --sza 30 --vza 30 --raa 90 --wind 3 \
              --chl 1.2 --tsm 3.5 --acdom440 0.04 --aod865 0.18 --aer C50
  배치:  python3 ocrt_solve.py grid --grid full_grid_design_v1.csv --out out.csv \
              [--bands 443 555 865] [--rows 0:10] [--which R1 R2 R3]
"""
import argparse
import csv
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
DATA = os.path.join(HERE, 'data')
from ocrt_py.constituent import PHYTO_GROUP_CHOICES, validate_phyto_group_option

BANDS_DEFAULT = [380, 412, 443, 490, 510, 555, 620, 660, 680, 709, 745, 865]

# ---- lazy singletons --------------------------------------------------------
_STATE = {}


def _absorption():
    if 'ab' not in _STATE:
        from ocrt_py.absorption import Absorption
        _STATE['ab'] = Absorption(DATA + '/afgl_atm', DATA + '/xsec')
    return _STATE['ab']


def _aer_mie(name):
    from ocrt_py.aerosol import read_mie
    key = 'aer:' + name
    if key not in _STATE:
        fn = name if name.endswith('.mie') else name + '.mie'
        _STATE[key] = read_mie(os.path.join(DATA, fn))
    return _STATE[key]


def _constituent_model(phyto_group, tsm_species):
    from ocrt_py.constituent import OCRTConstituentModel
    key = 'cm:%s:%s' % (phyto_group, tsm_species)
    if key not in _STATE:
        _STATE[key] = OCRTConstituentModel(DATA, phyto_group, tsm_species)
    return _STATE[key]


# ---- constituent Greek phase (b-weighted mix) -------------------------------
def _constituent_greek(cm, wl_nm, chl, tsm, acdom440, L=200, nmg=400,
                       cdom_slope=0.014, det_a440=0.0, det_slope=0.0109):
    from ocrt_py.aerosol import (water_component_phase_moments_gauss,
                                phase_matrix_at_wavelength_pchip)
    r = cm.evaluate(wl_nm, chl, tsm, acdom440, cdom_slope=cdom_slope,
                    det_a440=det_a440, det_slope=det_slope)
    b_phyto, b_det, b_min = r['b_phyto'], r['b_det'], r['b_min']
    b_particle = b_phyto + b_det + b_min
    be = np.zeros(L + 1); ga = np.zeros(L + 1)
    al = np.zeros(L + 1); ze = np.zeros(L + 1)
    diagnostics = []
    components = (
        ('phyto', cm.phyto_mie, b_phyto, r['bb_phyto']),
        ('detritus', cm.det_mie, b_det, r['bb_det']),
        ('tsm', cm.min_mie, b_min, r['bb_min']),
    )
    if b_particle > 0.0:
        for label, mie, bc, bb_iop in components:
            if bc <= 0.0:
                continue
            b, g, a, z, bb_phase = water_component_phase_moments_gauss(
                mie, wl_nm, L, nmg)
            w = bc / b_particle
            be += w * b; ga += w * g; al += w * a; ze += w * z
            diagnostics.append((label, mie, bc, bb_iop, bb_phase, w))

    # Match the C debug contract for a one-component constituent.  The normal
    # production path performs no additional work and emits no text.
    dump_phase = (os.environ.get('OCRT_DEBUG', '0') not in ('', '0') and
                  os.environ.get('OCRT_DUMP_IOP', '0') not in ('', '0'))
    if dump_phase and diagnostics:
        for label, mie, bc, bb_iop, bb_phase, w in diagnostics:
            p11, p12, p33 = phase_matrix_at_wavelength_pchip(mie, wl_nm)
            order = np.argsort(mie.angles, kind='stable')
            ang = np.asarray(mie.angles, float)[order]
            q11 = float(np.interp(90.0, ang, np.asarray(p11)[order]))
            q12 = float(np.interp(90.0, ang, np.asarray(p12)[order]))
            q33 = float(np.interp(90.0, ang, np.asarray(p33)[order]))
            print(
                'OCRT_PHASE_COMPONENT name=%s weight_b=%.12g b=%.12g '
                'bb_b_iop=%.12g bb_b_phase=%.12g '
                'P90=[%.8e %.8e %.8e] cache=python path=%s' %
                (label, w, bc, bb_iop / bc, bb_phase,
                 q11, q12, q33, mie.source_path),
                file=sys.stderr)
        mu_back = (-1.0, -0.93969262078591, -0.86602540378444,
                   -0.70710678118655, -0.5)
        rmin = float('inf')
        for x in mu_back:
            rmin = min(rmin, float(np.polynomial.legendre.legval(x, be)))
        print(
            'OCRT_PHASE_MIX components=%d L=%d beta0=%.12g beta2=%.12g '
            'gamma2=%.12g alpha2=%.12g zeta2=%.12g '
            'recon_back_min=%.12g cache=python' %
            (len(diagnostics), L, be[0], be[2] if L >= 2 else 0.0,
             ga[2] if L >= 2 else 0.0, al[2] if L >= 2 else 0.0,
             ze[2] if L >= 2 else 0.0, rmin),
            file=sys.stderr)

    constituent = dict(betal=be, gammal=ga, alphal=al, zetal=ze,
                       b_w=r['b_w'], b_particle=b_particle)
    return r, constituent


# ---- three reflectances -----------------------------------------------------
def solve_reflectances(wl_nm, sza, vza, raa, wind, chl, tsm, acdom440,
                       aod865, aer_name, which=('R1', 'R2', 'R3'),
                       phyto_group='micro', tsm_species='red_clay',
                       n_mu_water=48, fourier_m_max=16, nt_atm=400,
                       max_it_water=500, pressure_hpa=1013.25,
                       gas_on=True, cdom_slope=0.014,
                       det_a440=0.0, det_slope=0.0109):
    """Return dict with rho_TOA/rho_R/rho_RA_black (I/Q/U), Rrs, aod_band, and
    derived rho_RC / rho_AplusRA / t_rho_w.  Only the requested subset runs."""
    from ocrt_py.atmos import solve_aerosol_black_fresnel_ocean_value
    from ocrt_py.driver import solve_case_ocean_coupled
    from ocrt_py.aerosol import aod_at_wavelength, prepare_aerosol_runtime
    mie = _aer_mie(aer_name)
    ab = _absorption() if gas_on else None
    aod_band = aod_at_wavelength(mie, wl_nm, aod865, 865.0)
    aer_runtime = prepare_aerosol_runtime(
        mie, wl_nm, aod865, aod_ref_nm=865.0, L_max=80)
    out = dict(wavelength_nm=wl_nm, sza=sza, vza=vza, raa=raa, wind=wind,
               chl=chl, tsm=tsm, acdom440=acdom440, aod865=aod865,
               aer_model=aer_name, aod_band=aod_band,
               adom_slope=cdom_slope, detritus_a440=det_a440,
               detritus_slope=det_slope)

    if 'R2' in which:
        r2 = solve_aerosol_black_fresnel_ocean_value(
            sza, vza, raa, wl_nm, wind, mie, 0.0, L_max=80,
            pressure_hpa=pressure_hpa, n_mu_gl=24, nt=nt_atm,
            m_max=2, max_iterations=100, absorption=ab, surface='black_fresnel_ocean')
        # decouple sunglint: solve_aerosol returns decoupled rho by default
        out['rho_R_I'], out['rho_R_Q'], out['rho_R_U'] = r2[0], r2[1], r2[2]

    if 'R3' in which:
        r3 = solve_aerosol_black_fresnel_ocean_value(
            sza, vza, raa, wl_nm, wind, mie, aod_band, L_max=80,
            pressure_hpa=pressure_hpa, n_mu_gl=24, nt=nt_atm,
            m_max=fourier_m_max, max_iterations=100, absorption=ab,
            surface='black_fresnel_ocean', aerosol_runtime=aer_runtime)
        out['rho_RA_black_I'], out['rho_RA_black_Q'], out['rho_RA_black_U'] = \
            r3[0], r3[1], r3[2]

    if 'R1' in which:
        cm = _constituent_model(phyto_group, tsm_species)
        _, constituent = _constituent_greek(cm, wl_nm, chl, tsm, acdom440,
                                             cdom_slope=cdom_slope,
                                             det_a440=det_a440, det_slope=det_slope)
        rr = cm.evaluate(wl_nm, chl, tsm, acdom440, cdom_slope=cdom_slope,
                         det_a440=det_a440, det_slope=det_slope)
        res = solve_case_ocean_coupled(
            sza, vza, raa, wl_nm, n_water=1.34, wind_speed=wind,
            a_tot=rr['a'], b_tot=rr['b'], bb_tot=rr['bb'],
            phase_lut_path=None, mie=mie, user_aod=aod_band,
            pressure_hpa=pressure_hpa, absorption=ab,
            n_mu_water=n_mu_water, fourier_m_max=fourier_m_max, nt_atm=nt_atm,
            max_it_water=max_it_water, constituent=constituent,
            aerosol_runtime=aer_runtime)
        out['rho_TOA_I'], out['rho_TOA_Q'], out['rho_TOA_U'] = \
            res['rho_I'], res['rho_Q'], res['rho_U']
        out['Rrs'] = res['Rrs0plus']
        out['Rrs_Q'] = res['Rrs0plus_Q']
        out['Rrs_U'] = res['Rrs0plus_U']
        out['rrs0minus'] = res['rrs0minus']
        out['rrs0minus_Q'] = res['rrs0minus_Q']
        out['rrs0minus_U'] = res['rrs0minus_U']
        out['orders'] = res['orders']
        out['conv'] = res['conv']

    # derived (only where both operands present)
    def g(k):
        return out.get(k, float('nan'))
    if 'R1' in which and 'R2' in which:
        out['rho_RC_I'] = g('rho_TOA_I') - g('rho_R_I')
    if 'R3' in which and 'R2' in which:
        out['rho_AplusRA_I'] = g('rho_RA_black_I') - g('rho_R_I')
    if 'R1' in which and 'R3' in which:
        out['t_rho_w_I'] = g('rho_TOA_I') - g('rho_RA_black_I')
    return out


# ---- grid CSV loader (shared aliases) ---------------------------------------
ALIASES = {
    'row_id': ['case_id', 'row_id', 'id', 'case', 'index', 'row'],
    'tsm': ['tsm', 'tsm_g_m3', 'tsm_gm3'],
    'chl': ['chl', 'chl_mg_m3', 'chl_mgm3', 'chla'],
    'acdom': ['acdom440', 'adom440', 'a_cdom_440', 'acdom440_m_inv', 'acdom_440',
              'cdom440', 'a_cdom440'],
    'aer': ['aer_model', 'aerosol_model', 'aerosol', 'mie', 'model', 'aer'],
    'aod865': ['aod865', 'aod_865'],
    'wind': ['wind', 'wind_ms', 'wind_speed', 'wind_speed_ms', 'wind_m_s'],
    'sza': ['sza', 'sza_deg'], 'vza': ['vza', 'vza_deg'], 'raa': ['raa', 'raa_deg'],
}


# v1.11 (2026-08-29): optional per-row constituent-slope columns.  Absent
# columns fall back to the CLI value, which defaults to the C defaults
# (S_cdom 0.014 nm^-1, S_detritus 0.0109 nm^-1, detritus a440 0), so existing
# grids keep bit-identical results.  Names accept the C batch spellings
# (ocrt_adom_slope / ocrt_detritus_*) as well as short forms.
OPTIONAL_ALIASES = {
    'adom_slope': ['ocrt_adom_slope', 'ccrr_adom_slope', 'adom_slope',
                   'cdom_slope', 's_cdom', 'scdom'],
    'det_a440': ['ocrt_detritus_a440', 'detritus_a440', 'det_a440'],
    'det_slope': ['ocrt_detritus_slope', 'detritus_slope', 'det_slope',
                  's_detritus', 'sdetritus', 's_det'],
}


def _norm(s):
    return re.sub(r'[^a-z0-9]', '', s.lower())


def load_grid(path):
    with open(path, newline='') as fp:
        rd = csv.reader(fp)
        header = next(rd)
        rows = [r for r in rd if r and any(c.strip() for c in r)]
    nmap = {_norm(h): i for i, h in enumerate(header)}
    col = {}
    for key, names in ALIASES.items():
        for n in names:
            if _norm(n) in nmap:
                col[key] = nmap[_norm(n)]
                break
    miss = [k for k in ALIASES if k not in col and k != 'row_id']
    if miss:
        raise SystemExit('grid CSV 필수 열 없음: %s (헤더 %s)' % (miss, header))
    ocol = {}
    for key, names in OPTIONAL_ALIASES.items():
        for n in names:
            if _norm(n) in nmap:
                ocol[key] = nmap[_norm(n)]
                break
    out = []
    for i, r in enumerate(rows):
        rid = r[col['row_id']].strip() if 'row_id' in col else str(i)
        rec = dict(row_id=rid, tsm=float(r[col['tsm']]),
                   chl=float(r[col['chl']]), acdom=float(r[col['acdom']]),
                   aer=r[col['aer']].strip(), aod865=float(r[col['aod865']]),
                   wind=float(r[col['wind']]), sza=float(r[col['sza']]),
                   vza=float(r[col['vza']]), raa=float(r[col['raa']]))
        for key, ci in ocol.items():
            cell = r[ci].strip() if ci < len(r) else ''
            if cell:
                rec[key] = float(cell)
        out.append(rec)
    return out


# ---- CLI --------------------------------------------------------------------
def cmd_single(a):
    t = time.time()
    r = solve_reflectances(
        a.wl, a.sza, a.vza, a.raa, a.wind, a.chl, a.tsm, a.acdom440,
        a.aod865, a.aer, which=tuple(a.which),
        phyto_group=a.phyto_group, tsm_species=a.tsm_species,
        n_mu_water=a.n_mu_water, fourier_m_max=a.m_max, nt_atm=a.n_layers,
        max_it_water=a.max_it_water, pressure_hpa=a.pressure, gas_on=not a.no_gas,
        cdom_slope=a.adom_slope, det_a440=a.detritus_a440,
        det_slope=a.detritus_slope)
    dt = time.time() - t
    print('# %.1f s, conv=%s orders=%s' % (dt, r.get('conv'), r.get('orders')))
    for k in ('rho_TOA_I', 'rho_R_I', 'rho_RA_black_I', 'rho_RC_I',
              'rho_AplusRA_I', 't_rho_w_I', 'Rrs', 'rrs0minus', 'aod_band',
              'rho_TOA_Q', 'rho_TOA_U'):
        if k in r:
            print('%-16s = %.10e' % (k, r[k]))


ASSEMBLE_COLS = ['row_id', 'wavelength_nm', 'tsm', 'chl', 'acdom440', 'aer_model',
                 'aod865', 'wind', 'sza', 'vza', 'raa', 'aod_band',
                 'rho_TOA_I', 'rho_R_I', 'rho_RA_black_I',
                 'rho_RC_I', 'rho_AplusRA_I', 't_rho_w_I', 'Rrs', 'rrs0minus',
                 'rho_TOA_Q', 'rho_TOA_U', 'rho_R_Q', 'rho_R_U',
                 'rho_RA_black_Q', 'rho_RA_black_U', 'orders', 'conv',
                 'adom_slope', 'detritus_a440', 'detritus_slope']


def cmd_grid(a):
    rows = load_grid(a.grid)
    if a.rows:
        lo, hi = a.rows.split(':')
        rows = rows[int(lo):int(hi)]
    bands = a.bands
    total = len(rows) * len(bands)
    print('rows=%d bands=%d -> %d cases  which=%s' %
          (len(rows), len(bands), total, a.which))
    write_header = not (a.resume and os.path.exists(a.out))
    done = set()
    if a.resume and os.path.exists(a.out):
        with open(a.out, newline='') as fp:
            for rec in csv.DictReader(fp):
                done.add((rec['row_id'], int(float(rec['wavelength_nm']))))
    mode = 'a' if (a.resume and os.path.exists(a.out)) else 'w'
    t0 = time.time()
    n = 0
    with open(a.out, mode, newline='') as fo:
        w = csv.writer(fo)
        if write_header:
            w.writerow(ASSEMBLE_COLS)
        for row in rows:
            for wl in bands:
                n += 1
                if (row['row_id'], wl) in done:
                    continue
                try:
                    r = solve_reflectances(
                        wl, row['sza'], row['vza'], row['raa'], row['wind'],
                        row['chl'], row['tsm'], row['acdom'], row['aod865'],
                        row['aer'], which=tuple(a.which),
                        phyto_group=a.phyto_group, tsm_species=a.tsm_species,
                        n_mu_water=a.n_mu_water, fourier_m_max=a.m_max,
                        nt_atm=a.n_layers, max_it_water=a.max_it_water,
                        pressure_hpa=a.pressure, gas_on=not a.no_gas,
                        cdom_slope=row.get('adom_slope', a.adom_slope),
                        det_a440=row.get('det_a440', a.detritus_a440),
                        det_slope=row.get('det_slope', a.detritus_slope))
                    r['row_id'] = row['row_id']
                    w.writerow([r.get(c, '') for c in ASSEMBLE_COLS])
                    fo.flush()
                except Exception as e:
                    print('FAIL row=%s wl=%d: %s' % (row['row_id'], wl, e),
                          flush=True)
                if n % a.progress == 0 or n == total:
                    el = time.time() - t0
                    print('[%d/%d] %.1f min' % (n, total, el / 60.0), flush=True)
    print('완료 -> %s' % a.out)


def main():
    ap = argparse.ArgumentParser(description='OCRT 파이썬 네이티브 반사도 산출')
    sub = ap.add_subparsers(dest='cmd', required=True)

    def common(p):
        p.add_argument('--which', nargs='+', default=['R1', 'R2', 'R3'],
                       choices=['R1', 'R2', 'R3'])
        p.add_argument('--phyto-group', default=None,
                       choices=PHYTO_GROUP_CHOICES,
                       help='omit for absorption-only Chl; explicit species with Chl>0 is rejected in stage-2')
        p.add_argument('--tsm-species', default='red_clay',
                       choices=['red_clay', 'brown_earth', 'yellow_clay',
                                'calcareous_sand'])
        p.add_argument('--n-mu-water', type=int, default=48)
        p.add_argument('--m-max', type=int, default=16)
        p.add_argument('--n-layers', type=int, default=400)
        p.add_argument('--max-it-water', type=int, default=500)
        p.add_argument('--pressure', type=float, default=1013.25)
        p.add_argument('--no-gas', action='store_true',
                       help='기체흡수 끄기(기본 켜짐)')
        # v1.11 (2026-08-29): constituent spectral-slope controls, mirroring the
        # C options --ocrt-adom-slope / --ocrt-detritus-a440 / --ocrt-detritus-slope.
        p.add_argument('--adom-slope', '--cdom-slope', dest='adom_slope',
                       type=float, default=0.014,
                       help='aDOM/CDOM 지수 기울기 S [nm^-1] (기본 0.014, Bricaud 1981)')
        p.add_argument('--detritus-a440', dest='detritus_a440',
                       type=float, default=0.0,
                       help='detritus a(440) [m^-1] (기본 0; Chl>0 에서만 의미)')
        p.add_argument('--detritus-slope', dest='detritus_slope',
                       type=float, default=0.0109,
                       help='detritus 지수 기울기 S [nm^-1] (기본 0.0109, Bricaud-Stramski)')

    ps = sub.add_parser('single')
    common(ps)
    ps.add_argument('--wl', type=float, required=True)
    ps.add_argument('--sza', type=float, required=True)
    ps.add_argument('--vza', type=float, required=True)
    ps.add_argument('--raa', type=float, required=True)
    ps.add_argument('--wind', type=float, required=True)
    ps.add_argument('--chl', type=float, default=0.0)
    ps.add_argument('--tsm', type=float, default=0.0)
    ps.add_argument('--acdom440', type=float, default=0.0)
    ps.add_argument('--aod865', type=float, default=0.0)
    ps.add_argument('--aer', default='C50')

    pg = sub.add_parser('grid')
    common(pg)
    pg.add_argument('--grid', required=True)
    pg.add_argument('--out', required=True)
    pg.add_argument('--bands', type=int, nargs='+', default=BANDS_DEFAULT)
    pg.add_argument('--rows', default='')
    pg.add_argument('--resume', action='store_true')
    pg.add_argument('--progress', type=int, default=1)

    a = ap.parse_args()
    explicit_phyto = a.phyto_group is not None
    if a.cmd == 'single':
        a.phyto_group = validate_phyto_group_option(
            a.phyto_group, chl=a.chl, explicit=explicit_phyto)
        cmd_single(a)
    else:
        rows_for_gate = load_grid(a.grid)
        if explicit_phyto and any(float(r['chl']) > 0.0 for r in rows_for_gate):
            validate_phyto_group_option(a.phyto_group, chl=1.0, explicit=True)
        a.phyto_group = validate_phyto_group_option(
            a.phyto_group, chl=0.0, explicit=explicit_phyto)
        cmd_grid(a)


if __name__ == '__main__':
    main()
