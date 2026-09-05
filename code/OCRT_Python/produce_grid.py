#!/usr/bin/env python3
"""GOCI-III polarized atmospheric-correction production driver.

Reads a grid CSV (one row per water/atmosphere/geometry condition) and, for
each of the 12 GOCI-III bands, runs the three radiance components as batched
(optionally GPU) solves:

  R1  rho_TOA   : ocean surface + atmosphere + aerosol + water-leaving
  R2  rho_R     : Cox-Munk fresnel (black) + Rayleigh, aod=0
  R3  rho_R+A   : Cox-Munk fresnel (black) + Rayleigh + aerosol, aod>0

and writes the difference products per (case, band):

  rho_RC       = rho_TOA - rho_R
  rho_A+rho_RA = rho_R+A - rho_R
  t_rho_w      = rho_TOA - rho_R+A
  Rrs_I/Q/U    = water-leaving Stokes / Ed(0+) from R1
  rrs_I/Q/U    = in-water upwelling Stokes / Ed(0-) from R1
  Rrs          = backward-compatible alias of Rrs_I

GPU: set env OCRT_PY_GPU=1 to use CuPy (backend.py).  Resume: rows already in
the output CSV are skipped.  Progress is printed per band and per chunk.

Grid CSV columns (header, case-insensitive):
  case_id, sza, vza, raa, wind, aerosol, aod865, chl, tsm, acdom440
where aerosol is a .mie basename in data/ (e.g. C50, r80f50v01).

NOTE (accuracy items still pending, see --strict-check):
  * gas absorption is not yet applied in the batch path
  * R2/R3 use the Greek-coefficient kernel; C production uses the aerosol
    value kernel (aer_use_value_kernel=1) for the Cox-Munk-only runs.
  These affect absolute agreement with the C production binary and must be
  wired in before the CSV is treated as final.
"""
import os
import sys
import csv
import argparse
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from ocrt_py import aerosol as AER
from ocrt_py import batch_driver as BD
from ocrt_py.constituent import OCRTConstituentModel

BANDS_NM = [380, 412, 443, 490, 510, 555, 620, 660, 680, 709, 745, 865]
AOD_REF_NM = 865.0   # grid specifies AOD at 865 nm

OUT_FIELDS = ['case_id', 'band_nm',
              'sza', 'vza', 'raa', 'wind', 'aerosol', 'aod865',
              'chl', 'tsm', 'acdom440',
              'rho_TOA', 'rho_R', 'rho_RpA', 'rho_RC', 'rho_A_RA', 't_rho_w',
              'T_dir_dn', 'T_diff_dn_hemi', 'T_total_dn_hemi',
              'T_dir_up_view', 'T_diff_up_view', 'T_total_up_view',
              'TOA_water_signal_I', 'T_up_rt_valid',
              'rho_wn', 'Rrs', 'Rrs_I', 'Rrs_Q', 'Rrs_U',
              'rrs_I', 'rrs_Q', 'rrs_U',
              'a_total', 'b_total', 'bb_total', 'a_w', 'b_w', 'bb_w',
              'a_chl', 'a_phyto_detritus', 'b_phyto_detritus', 'bb_phyto_detritus',
              'a_dom', 'a_min', 'b_min', 'bb_min', 'Kd0minus']

from ocrt_py.constituent import validate_phyto_group_option


def _read_grid(path):
    rows = []
    with open(path, newline='') as fh:
        rd = csv.DictReader(fh)
        cmap = {k.lower(): k for k in rd.fieldnames}
        req = ['case_id', 'sza', 'vza', 'raa', 'wind', 'aerosol',
               'aod865', 'chl', 'tsm', 'acdom440']
        miss = [c for c in req if c not in cmap]
        if miss:
            raise SystemExit(f'grid CSV missing columns: {miss}\n'
                             f'have: {rd.fieldnames}')
        for r in rd:
            rec = dict(
                case_id=str(r[cmap['case_id']]),
                sza=float(r[cmap['sza']]), vza=float(r[cmap['vza']]),
                raa=float(r[cmap['raa']]), wind=float(r[cmap['wind']]),
                aerosol=str(r[cmap['aerosol']]).strip(),
                aod865=float(r[cmap['aod865']]), chl=float(r[cmap['chl']]),
                tsm=float(r[cmap['tsm']]), acdom440=float(r[cmap['acdom440']]))
            # v1.11 (2026-08-29): optional constituent spectral-slope columns.
            # Absent columns keep the C defaults (S_cdom 0.014, S_det 0.0109,
            # detritus a440 0), so existing grids stay bit-identical.
            for key, names in (('cdom_slope', ('ocrt_adom_slope', 'adom_slope',
                                               'cdom_slope')),
                               ('det_a440', ('ocrt_detritus_a440',
                                             'detritus_a440', 'det_a440')),
                               ('det_slope', ('ocrt_detritus_slope',
                                              'detritus_slope', 'det_slope'))):
                for nm in names:
                    if nm in cmap and str(r[cmap[nm]]).strip():
                        rec[key] = float(r[cmap[nm]])
                        break
            rows.append(rec)
    return rows


def _load_done(out_path):
    """Return completed keys and reject incompatible legacy output schemas."""
    done = set()
    if not os.path.exists(out_path):
        return done
    with open(out_path, newline='') as fh:
        rd = csv.DictReader(fh)
        if rd.fieldnames != OUT_FIELDS:
            raise SystemExit(
                'output CSV schema differs from this pyOCRT version; '
                'use a new --out file (Rrs/rrs I/Q/U columns were added).')
        for r in rd:
            done.add((r['case_id'], int(float(r['band_nm']))))
    return done


def _mie_cache(rows, data_dir):
    names = sorted({r['aerosol'] for r in rows})
    cache = {}
    for nm in names:
        cache[nm] = AER.read_mie(os.path.join(data_dir, nm + '.mie'))
    return cache


def _band_cases(rows, band_nm, mie_cache, kind):
    """Build per-band case dicts for solve_r{1,2,3}_grid.  kind in
    {'r1','r2','r3'}.  AOD is scaled from 865 nm to the band; R2 uses aod=0."""
    cases = []
    for r in rows:
        if kind == 'r2':
            aod = 0.0
        else:
            mie = mie_cache[r['aerosol']]
            aod = AER.aod_at_wavelength(mie, band_nm, r['aod865'], AOD_REF_NM)
        c = dict(sza=r['sza'], vza=r['vza'], raa=r['raa'], wind=r['wind'],
                 aer=r['aerosol'], wl=band_nm, aod=aod)
        if kind == 'r1':
            c.update(chl=r['chl'], tsm=r['tsm'], ad=r['acdom440'])
            for key in ('cdom_slope', 'det_a440', 'det_slope'):
                if key in r:
                    c[key] = r[key]
        cases.append(c)
    return cases


def _chunks(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def produce(grid_path, out_path, data_dir, cm_phyto, cm_min,
            chunk=64, n_mu_water=24, nt_atm=400, fourier_m_max=16,
            L_max=80, Lmix=200, max_it_atm=100, max_it_water=500,
            tol=1.0e-7, n_water=1.34, pressure_hpa=1013.25, bands=None,
            gas=True):
    rows = _read_grid(grid_path)
    if cm_phyto is None:
        cm_phyto = 'micro'
    bands = bands or BANDS_NM
    mie_cache = _mie_cache(rows, data_dir)
    cm = OCRTConstituentModel(data_dir, cm_phyto, cm_min)
    absorption = None
    if gas:
        from ocrt_py.absorption import Absorption
        absorption = Absorption(os.path.join(data_dir, 'afgl_atm'),
                                os.path.join(data_dir, 'xsec'))
    done = _load_done(out_path)
    new_file = not os.path.exists(out_path)
    fh = open(out_path, 'a', newline='')
    wr = csv.DictWriter(fh, fieldnames=OUT_FIELDS)
    if new_file:
        wr.writeheader(); fh.flush()

    total = len(rows) * len(bands)
    done_ct = len(done)
    t0 = time.time()
    print(f'[produce] {len(rows)} cases x {len(bands)} bands = {total} '
          f'(already done: {done_ct}); GPU={os.environ.get("OCRT_PY_GPU","0")} '
          f'gas={"on" if gas else "off"}')

    for band in bands:
        pend = [r for r in rows if (r['case_id'], band) not in done]
        if not pend:
            print(f'  band {band} nm: all done, skip')
            continue
        print(f'  band {band} nm: {len(pend)} pending', flush=True)
        for chunk_rows in _chunks(pend, chunk):
            lo = done_ct + 1
            hi = done_ct + len(chunk_rows)
            print(f'    rows {lo}-{hi}/{total} 처리 중...', end='', flush=True)
            c2 = _band_cases(chunk_rows, band, mie_cache, 'r2')
            c3 = _band_cases(chunk_rows, band, mie_cache, 'r3')
            c1 = _band_cases(chunk_rows, band, mie_cache, 'r1')
            rR = BD.solve_r2_grid(c2, pressure_hpa=pressure_hpa,
                                  n_mu_gl=n_mu_water, absorption=absorption)
            rRA = BD.solve_r3_grid(c3, mie_cache, L_max=L_max,
                                   pressure_hpa=pressure_hpa, n_mu_gl=n_mu_water,
                                   nt=nt_atm, m_max=fourier_m_max,
                                   max_iterations=max_it_atm, tolerance=tol,
                                   n_water=n_water, absorption=absorption)
            rTOA, rrsTOA, trA, iopA = BD.solve_r1_grid(c1, mie_cache, cm, L_max=L_max, Lmix=Lmix,
                                            pressure_hpa=pressure_hpa,
                                            n_mu_water=n_mu_water, nt_atm=nt_atm,
                                            fourier_m_max=fourier_m_max,
                                            max_it_atm=max_it_atm, tol_atm=tol,
                                            max_it_water=max_it_water, tol_water=tol,
                                            n_water=n_water, absorption=absorption)
            for k, r in enumerate(chunk_rows):
                rc = rTOA[k] - rR[k]
                ara = rRA[k] - rR[k]
                trw = rTOA[k] - rRA[k]
                tdn = trA['T_total_dn_hemi'][k]; tup = trA['T_total_up_view'][k]
                denom = tdn * tup
                rwn = (trw / denom) if (denom == denom and denom != 0.0) else float('nan')
                wr.writerow(dict(case_id=r['case_id'], band_nm=band,
                                 sza=f"{r['sza']:.4f}", vza=f"{r['vza']:.4f}",
                                 raa=f"{r['raa']:.4f}", wind=f"{r['wind']:.4f}",
                                 aerosol=r['aerosol'], aod865=f"{r['aod865']:.6f}",
                                 chl=f"{r['chl']:.6f}", tsm=f"{r['tsm']:.6f}",
                                 acdom440=f"{r['acdom440']:.6f}",
                                 rho_TOA=f'{rTOA[k]:.10e}', rho_R=f'{rR[k]:.10e}',
                                 rho_RpA=f'{rRA[k]:.10e}', rho_RC=f'{rc:.10e}',
                                 rho_A_RA=f'{ara:.10e}', t_rho_w=f'{trw:.10e}',
                                 T_dir_dn=f"{trA['T_dir_dn'][k]:.10e}",
                                 T_diff_dn_hemi=f"{trA['T_diff_dn_hemi'][k]:.10e}",
                                 T_total_dn_hemi=f'{tdn:.10e}',
                                 T_dir_up_view=f"{trA['T_dir_up_view'][k]:.10e}",
                                 T_diff_up_view=f"{trA['T_diff_up_view'][k]:.10e}",
                                 T_total_up_view=f'{tup:.10e}',
                                 TOA_water_signal_I=f"{trA['TOA_water_signal_I'][k]:.10e}",
                                 T_up_rt_valid=int(trA['T_up_rt_valid'][k]),
                                 rho_wn=f'{rwn:.10e}', Rrs=f'{rrsTOA[k]:.10e}',
                                 Rrs_I=f"{trA['Rrs_I'][k]:.10e}",
                                 Rrs_Q=f"{trA['Rrs_Q'][k]:.10e}",
                                 Rrs_U=f"{trA['Rrs_U'][k]:.10e}",
                                 rrs_I=f"{trA['rrs_I'][k]:.10e}",
                                 rrs_Q=f"{trA['rrs_Q'][k]:.10e}",
                                 rrs_U=f"{trA['rrs_U'][k]:.10e}",
                                 a_total=f"{iopA['a_total'][k]:.8e}",
                                 b_total=f"{iopA['b_total'][k]:.8e}",
                                 bb_total=f"{iopA['bb_total'][k]:.8e}",
                                 a_w=f"{iopA['a_w'][k]:.8e}",
                                 b_w=f"{iopA['b_w'][k]:.8e}",
                                 bb_w=f"{iopA['bb_w'][k]:.8e}",
                                 a_chl=f"{iopA['a_chl'][k]:.8e}",
                                 a_phyto_detritus=f"{iopA['a_phyto_detritus'][k]:.8e}",
                                 b_phyto_detritus=f"{iopA['b_phyto_detritus'][k]:.8e}",
                                 bb_phyto_detritus=f"{iopA['bb_phyto_detritus'][k]:.8e}",
                                 a_dom=f"{iopA['a_dom'][k]:.8e}",
                                 a_min=f"{iopA['a_min'][k]:.8e}",
                                 b_min=f"{iopA['b_min'][k]:.8e}",
                                 bb_min=f"{iopA['bb_min'][k]:.8e}",
                                 Kd0minus=f"{iopA['Kd0minus'][k]:.6e}"))
            fh.flush()
            done_ct += len(chunk_rows)
            el = time.time() - t0
            rate = done_ct / el if el > 0 else 0
            remain = total - done_ct
            eta = remain / rate if rate > 0 else 0
            # \r로 같은 줄 덮어써 청크 완료 요약 표시
            print(f'\r    {done_ct}/{total} 완료  경과 {el:.0f}s  '
                  f'{rate:.2f} rows/s  ETA {eta/60:.1f}분        ', flush=True)
    fh.close()
    print(f'[produce] done: {done_ct}/{total} in {time.time()-t0:.0f}s -> {out_path}')


def main():
    ap = argparse.ArgumentParser(description='GOCI-III AC production driver')
    ap.add_argument('--grid', required=True, help='grid CSV path')
    ap.add_argument('--out', required=True, help='output CSV path (append/resume)')
    ap.add_argument('--data', default=os.path.join(_HERE, 'data'),
                    help='data dir with .mie and constituent files')
    ap.add_argument('--phyto', default=None,
                    help='omit for absorption-only Chl; explicit species with Chl>0 is rejected in stage-2')
    ap.add_argument('--mineral', default='red_clay', help='mineral model key')
    ap.add_argument('--chunk', type=int, default=64, help='rows per batch chunk')
    ap.add_argument('--n-mu-water', type=int, default=24)
    ap.add_argument('--nt-atm', type=int, default=400)
    ap.add_argument('--m-max', type=int, default=16)
    ap.add_argument('--max-it-water', type=int, default=500)
    ap.add_argument('--bands', default='', help='comma list of band nm (default all 12)')
    ap.add_argument('--no-gas', action='store_true', help='disable gas absorption')
    args = ap.parse_args()
    explicit_phyto = args.phyto is not None
    gate_rows = _read_grid(args.grid)
    if explicit_phyto and any(float(r['chl']) > 0.0 for r in gate_rows):
        validate_phyto_group_option(args.phyto, chl=1.0, explicit=True)
    args.phyto = validate_phyto_group_option(
        args.phyto, chl=0.0, explicit=explicit_phyto)
    bands = [int(x) for x in args.bands.split(',')] if args.bands else None
    produce(args.grid, args.out, args.data, args.phyto, args.mineral,
            chunk=args.chunk, n_mu_water=args.n_mu_water, nt_atm=args.nt_atm,
            fourier_m_max=args.m_max, max_it_water=args.max_it_water, bands=bands,
            gas=not args.no_gas)


if __name__ == '__main__':
    main()
