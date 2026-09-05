# -*- coding: utf-8 -*-
"""OCRT/CCRR constituent-model scalar IOP assembly.

Transliteration of rt_water_iop.c (pure water, CDOM), rt_iop_organic.c
(EAP phytoplankton, detritus), rt_iop_ahn_mineral.c (mineral), and the
orchestration in rt_water_rt.c (rt_iop_total + Chl+TSM mixed bb override).

Scalar path only (a, b, bb).  Phase-function mixing is in a separate module.
The bb/b ratio matches build_aerosol_interpolators (mie_io.c): P11 linterp
over phase wavelengths per angle, PCHIP in angle, 200-node midpoint
integration of P11*sin(theta) over the backward vs full hemisphere.
"""
import numpy as np

from .aerosol import read_mie, Pchip
from .spectral_contract import (SPECTRAL_MIN_NM, SPECTRAL_MAX_NM,
                                require_wavelength, require_table_coverage,
                                require_query_in_table)

ORGANIC_DETRITUS_SLOPE_DEFAULT = 0.0109
CDOM_SLOPE_DEFAULT = 0.014
PHYTO_ABS_ZERO_START_NM = 800.0
PHYTO_ABS_TABLE_MAX_NM = SPECTRAL_MAX_NM

# canonical file names (rt_iop_organic.c / rt_iop_ahn_mineral.c)
PHYTO_FILE = {
    # Backward-compatible aggregate aliases, all backed by the new catalog.
    'pico': 'eap/EAP_15_Synechococcus_D1p2.mie',
    'nano': 'eap/EAP_12_Hapto_Prymnesiaceae_D4.mie',
    'micro': 'eap/EAP_02_Diatoms_centric_D6.mie',
    # Canonical 17-species catalog.
    'eap_diatoms_pennate': 'eap/EAP_00_Diatoms_pennate_D6.mie',
    'eap_chlorophytes': 'eap/EAP_01_Chlorophytes_D8.mie',
    'eap_diatoms_centric': 'eap/EAP_02_Diatoms_centric_D6.mie',
    'eap_cryptophytes': 'eap/EAP_03_Cryptophytes_D6.mie',
    'eap_cyano_blue': 'eap/EAP_04_Cyano_blue_D6.mie',
    'eap_cyano_red': 'eap/EAP_05_Cyano_red_D6.mie',
    'eap_dinoflagellates': 'eap/EAP_06_Dinoflagellates_D24.mie',
    'eap_eustigmatophytes': 'eap/EAP_07_Eustigmatophytes_D6.mie',
    'eap_hapto_pavlovaceae': 'eap/EAP_08_Hapto_Pavlovaceae_D6.mie',
    'eap_pelagophytes': 'eap/EAP_09_Pelagophytes_D3.mie',
    'eap_prasinophytes': 'eap/EAP_10_Prasinophytes_D3.mie',
    'eap_prochlorococcus': 'eap/EAP_11_Prochlorococcus_D0p5.mie',
    'eap_hapto_prymnesiaceae': 'eap/EAP_12_Hapto_Prymnesiaceae_D4.mie',
    'eap_raphidophytes': 'eap/EAP_13_Raphidophytes_D24.mie',
    'eap_rhodophytes': 'eap/EAP_14_Rhodophytes_D6.mie',
    'eap_synechococcus': 'eap/EAP_15_Synechococcus_D1p2.mie',
    'eap_microcystis': 'eap/EAP_16_Microcystis_D5.mie',
}
PHYTO_GROUP_CHOICES = tuple(PHYTO_FILE.keys())
DEFAULT_PHYTO_GROUP = 'micro'
DEFAULT_PHYTO_CANONICAL = 'eap_diatoms_centric'

def advanced_mode_enabled():
    import os
    value = os.environ.get('OCRT_ADVANCED', '')
    return bool(value and value != '0')


def phyto_group_is_default(group):
    return group in (DEFAULT_PHYTO_GROUP, DEFAULT_PHYTO_CANONICAL)


def validate_phyto_group_option(group, chl=None, explicit=True):
    """Validate/resolve a CLI phytoplankton selection.

    Stage-2 production contract: the 17-species EAP generator/catalog is
    available, but species-specific phytoplankton scattering is not yet
    enabled in the OCRT constituent water model.  A CLI selection supplied
    explicitly while Chl>0 therefore fails loudly.  With no explicit option,
    the historical ``micro`` file remains the absorption-spectrum source and
    phytoplankton b/bb are forced to zero.
    """
    if group is None:
        return DEFAULT_PHYTO_GROUP
    if group not in PHYTO_FILE:
        raise ValueError('unsupported phyto_group: %s' % group)
    if explicit and chl is not None and float(chl) > 0.0:
        raise ValueError(
            '--phyto-group is not supported yet for Chl>0 in the OCRT '
            'constituent water model; species-specific EAP scattering is '
            'disabled (requested %r). Omit the option to use absorption-only '
            'Chl, or use the fixed-bulk IOP phase path.' % group)
    if not phyto_group_is_default(group) and not advanced_mode_enabled():
        raise ValueError(
            'phyto_group %s is an advanced EAP catalog selection; set '
            'OCRT_ADVANCED=1 (the constituent Chl scattering path remains '
            'disabled)' % group)
    return group

DETRITUS_FILE = 'Detritus_Stramski2001.mie'
AHN_FILE = {'red_clay': 'Red_clay_AHN.mie',
            'brown_earth': 'Brown_earth_AHN.mie',
            'yellow_clay': 'Yellow_clay_AHN.mie',
            'calcareous_sand': 'Calcareous_sand_AHN.mie'}


# ---- pure water + CDOM ------------------------------------------------------
def _load_end_header(path, ncol):
    """Read a '/end_header'-terminated whitespace table with ncol columns."""
    rows = []
    in_data = False
    with open(path) as fp:
        for ln in fp:
            if not in_data:
                if ln.strip().startswith('/end_header'):
                    in_data = True
                continue
            parts = ln.split()
            if len(parts) >= ncol:
                try:
                    rows.append([float(parts[i]) for i in range(ncol)])
                except ValueError:
                    continue
    return np.array(rows)


class WaterCoefLUT:
    """pure-water a_w, b_w LUT (water_coef_z09_1nm.txt: lambda_nm aw bw)."""
    def __init__(self, path):
        t = _load_end_header(path, 3)
        self.wl = require_table_coverage(t[:, 0], context='pure-water IOP')
        self.aw = t[:, 1]
        self.bw = t[:, 2]


class PsiTLUT:
    """psi_T LUT (lambda_nm psi_T), T_ref = 20 C."""
    def __init__(self, path, t_ref=20.0):
        t = _load_end_header(path, 2)
        self.wl = require_table_coverage(t[:, 0], context='psi_T')
        self.psi = t[:, 1]
        self.t_ref = t_ref


class PhytoAbsorptionLUT:
    """Validated Stage-2 default phytoplankton absorption spectrum.

    This is deliberately separate from the EAP phase-function catalog.  The
    catalog/generator is retained for future vector-scattering validation,
    while the production constituent model uses this frozen a* spectrum and
    forces phytoplankton b=bb=0.
    """
    def __init__(self, path):
        t = np.loadtxt(path, delimiter=',', comments='#', dtype=float)
        if t.ndim != 2 or t.shape[1] < 2 or t.shape[0] < 2:
            raise ValueError('invalid phytoplankton absorption LUT: %s' % path)
        if not np.all(np.diff(t[:, 0]) > 0.0):
            raise ValueError('non-monotonic phytoplankton absorption LUT: %s' % path)
        self.wl = require_table_coverage(
            t[:, 0], context='phytoplankton absorption')
        self.astar = t[:, 1]
        if not np.any(np.isclose(self.wl, PHYTO_ABS_ZERO_START_NM,
                                 atol=1e-9, rtol=0.0)):
            raise ValueError('missing 800-nm phytoplankton zero anchor: %s' % path)
        if np.any(self.astar[self.wl >= PHYTO_ABS_ZERO_START_NM - 1e-9] != 0.0):
            raise ValueError('non-zero phytoplankton absorption at or above 800 nm: %s' % path)

    def eval(self, wl_nm):
        wl = require_query_in_table(
            wl_nm, self.wl, context='phytoplankton absorption')
        return max(0.0, float(np.interp(wl, self.wl, self.astar)))


def pure_water(aw_lut, psi_lut, wl_nm, t_c=20.0):
    """a = aw20 + psi*(T-T_ref); b = bw; bb = 0.5*b (exact)."""
    wl_nm = require_wavelength(wl_nm, context='pure water')
    require_query_in_table(wl_nm, aw_lut.wl, context='pure-water IOP')
    aw20 = float(np.interp(wl_nm, aw_lut.wl, aw_lut.aw))
    bw = float(np.interp(wl_nm, aw_lut.wl, aw_lut.bw))
    a = aw20
    if psi_lut is not None and len(psi_lut.wl) > 0:
        require_query_in_table(wl_nm, psi_lut.wl, context='psi_T')
        psi = float(np.interp(wl_nm, psi_lut.wl, psi_lut.psi))
        a = aw20 + psi * (t_c - psi_lut.t_ref)
    return a, bw, 0.5 * bw


def cdom(a440, wl_nm, slope=CDOM_SLOPE_DEFAULT, ref_nm=440.0):
    """a_cdom = a440*exp(-S*(lambda-ref)); b=bb=0.  a440<=0 => off."""
    if not (a440 > 0.0):
        return 0.0, 0.0, 0.0
    a = a440 * np.exp(-slope * (wl_nm - ref_nm))
    return max(0.0, a), 0.0, 0.0


# ---- phase backscatter ratio (build_aerosol_interpolators) ------------------
def bb_b_ratio(mie, wl_nm, _cache={}):
    """Backscattering ratio from the delivered FR631 phase representation.

    Phase wavelength interpolation is linear. The ratio is integrated directly
    on the native nonuniform angle grid with a trapezoid in theta, rather than
    rebuilding a cubic/PCHIP angular curve. Cache key follows the stable
    (source_path, mtime_ns, size) identity delivered by the FR631 loader.
    Component gating checks from the 2026-08-14 tree are retained.
    """
    source_key = getattr(mie, 'cache_key', ('memory', id(mie)))
    key = (tuple(source_key), round(wl_nm, 6))
    if key in _cache:
        return _cache[key]
    wl_nm = require_wavelength(wl_nm, context='particle phase')
    wl_um = wl_nm / 1000.0
    require_query_in_table(wl_um, mie.phase_wavelengths,
                           context='particle phase grid', unit='um')
    order = np.argsort(mie.angles, kind='stable')
    ang = np.asarray(mie.angles, float)[order]
    P = np.asarray(mie.P11, float)[order, :]
    P_at = np.array([np.interp(wl_um, mie.phase_wavelengths, P[k, :])
                     for k in range(len(ang))], dtype=float)
    theta = np.deg2rad(ang)
    integrand = P_at * np.sin(theta)
    b_int = float(np.trapezoid(integrand, theta))
    mask = ang >= 90.0
    bb_int = float(np.trapezoid(integrand[mask], theta[mask]))
    r = (bb_int / b_int) if b_int > 0.0 else 0.0
    _cache[key] = r
    return r


def _star_ab(mie, wl_nm):
    """astar = Extinct_Co - Scatter_Co, bstar = Scatter_Co, linterp over
    the spectral wavelength grid (um).  spectral cols: 4=Extinct, 5=Scatter."""
    wl_nm = require_wavelength(wl_nm, context='particle bulk IOP')
    wl_um = wl_nm / 1000.0
    require_query_in_table(wl_um, mie.wavelengths,
                           context='particle bulk grid', unit='um')
    ext = float(np.interp(wl_um, mie.wavelengths, mie.spectral[:, 4]))
    sca = float(np.interp(wl_um, mie.wavelengths, mie.spectral[:, 5]))
    return max(0.0, ext - sca), sca


# ---- Huot / f_ph closures ---------------------------------------------------
def fph_fraction(chl):
    if not (chl > 0.0):
        return 0.020
    return 0.035 + 0.015 * np.tanh(np.log10(chl))


def huot_bbp(wl_nm, chl):
    if not (chl > 0.0):
        return 0.0
    a1 = 2.267e-3 - 5.058e-6 * (wl_nm - 550.0)
    b1 = 0.565 + 0.000486 * (wl_nm - 550.0)
    if not (a1 > 0.0):
        return 0.0
    return a1 * chl ** b1


# ---- EAP phytoplankton / detritus / AHN mineral -----------------------------
def eap_phyto(phyto_abs_lut, wl_nm, chl, group='micro'):
    """Stage-2 phytoplankton absorption-only production closure.

    The EAP catalog/generator remains available, but it is *not* the source of
    production constituent-model absorption or scattering in Stage 2.  The
    frozen validated a* spectrum is read from
    ``phyto_absorption_default.csv`` and phytoplankton b=bb=0.
    """
    if not (chl > 0.0):
        return 0.0, 0.0, 0.0
    return phyto_abs_lut.eval(wl_nm) * chl, 0.0, 0.0


def detritus(det_mie, wl_nm, chl, a_d440=0.0, slope=ORGANIC_DETRITUS_SLOPE_DEFAULT):
    """a = a_d440*exp(-S*(lambda-440)); bb = (1-f_ph)*huot_bbp(550)*(r_l/r_550);
    b = bb/r_l.

    The active production phase is the previously validated 350--850 nm
    Stramski table.  The recipe-regenerated 330--1100 candidate is retained
    under ``data/water_iop/candidates`` but is not activated because it failed
    the L=200/SOS matched-input acceptance test.
    """
    if not (chl > 0.0):
        return 0.0, 0.0, 0.0
    require_query_in_table(wl_nm, det_mie.phase_wavelengths * 1000.0,
                           context='active Chl-linked detritus phase')
    r = bb_b_ratio(det_mie, wl_nm)
    r550 = bb_b_ratio(det_mie, 550.0)
    a = a_d440 * np.exp(-slope * (wl_nm - 440.0))
    bb = (1.0 - fph_fraction(chl)) * huot_bbp(550.0, chl) * (r / r550)
    b = bb / r
    return a, b, bb


def ahn_mineral(min_mie, wl_nm, tsm, chl_positive):
    """a = a*_min*TSM; b = b*_min*TSM.  Mixed (Chl>0 & TSM>0): bb = b * r
    (PCHIP ratio, overriding the TSM-only direct-integral path).  Our grid
    always has Chl>0 & TSM>0 so the mixed branch is used."""
    if not (tsm > 0.0):
        return 0.0, 0.0, 0.0
    astar, bstar = _star_ab(min_mie, wl_nm)
    a = astar * tsm
    b = bstar * tsm
    r = bb_b_ratio(min_mie, wl_nm)
    bb = b * r
    return a, b, bb


# ---- assembly ---------------------------------------------------------------
class OCRTConstituentModel:
    """OCRT constituent scalar-IOP assembler.  Loads pure-water/psi LUTs and
    the requested particle .mie files once; evaluate() returns total a/b/bb
    plus component diagnostics at a wavelength."""
    def __init__(self, inputs_dir, phyto_group='micro', tsm_species='red_clay'):
        wi = inputs_dir + '/water_iop'
        ta = inputs_dir + '/tsm_ahn'
        self.aw_lut = WaterCoefLUT(wi + '/water_coef_z09_1nm.txt')
        try:
            self.psi_lut = PsiTLUT(wi + '/psi_T_rottgers2014_OCRT.txt')
        except Exception:
            self.psi_lut = None
        if phyto_group not in PHYTO_FILE:
            raise ValueError('unsupported phyto_group: %s' % phyto_group)
        # Keep the selected EAP file available for the catalog/generator and
        # future fixed-phase work, but do not use it in the Stage-2 production
        # constituent closure.
        self.phyto_mie = read_mie(wi + '/' + PHYTO_FILE[phyto_group])
        # EAP species-specific phase is intentionally outside the official
        # 330--1100 nm production contract.  It remains loaded for the
        # advanced catalog/generator but may only be queried inside its native
        # wavelength range; generic phase helpers enforce that fail-loud rule.
        self.phyto_abs_lut = PhytoAbsorptionLUT(
            wi + '/phyto_absorption_default.csv')
        self.det_mie = read_mie(wi + '/' + DETRITUS_FILE)
        if len(self.det_mie.phase_wavelengths) < 2:
            raise ValueError('invalid active detritus phase table')
        # Intentionally do not require 330--1100 coverage here.  The active
        # validated detritus phase is native 350--850 nm and each Chl query is
        # checked fail-loud by detritus(); the rejected 330--1100 candidate is
        # distributed separately for further research.
        self.min_mie = read_mie(ta + '/' + AHN_FILE[tsm_species])
        require_table_coverage(self.min_mie.wavelengths * 1000.0,
                               context='AHN mineral bulk Mie')
        require_table_coverage(self.min_mie.phase_wavelengths * 1000.0,
                               context='AHN mineral phase Mie')
        self.phyto_group = phyto_group
        self.tsm_species = tsm_species

    def evaluate(self, wl_nm, chl, tsm, acdom440,
                 t_c=20.0, cdom_slope=CDOM_SLOPE_DEFAULT,
                 det_a440=0.0, det_slope=ORGANIC_DETRITUS_SLOPE_DEFAULT):
        wl_nm = require_wavelength(wl_nm, context='OCRT constituent model')
        aw, bw, bbw = pure_water(self.aw_lut, self.psi_lut, wl_nm, t_c)
        ac, bc, bbc = cdom(acdom440, wl_nm, cdom_slope)
        ap, bp, bbp = eap_phyto(self.phyto_abs_lut, wl_nm, chl, self.phyto_group)
        ad, bd, bbd = detritus(self.det_mie, wl_nm, chl, det_a440, det_slope)
        am, bm, bbm = ahn_mineral(self.min_mie, wl_nm, tsm, chl > 0.0)
        a = aw + ac + ap + ad + am
        b = bw + bc + bp + bd + bm
        bb = bbw + bbc + bbp + bbd + bbm
        return dict(a=a, b=b, bb=bb,
                    a_w=aw, b_w=bw, bb_w=bbw, a_cdom=ac,
                    a_phyto=ap, b_phyto=bp, bb_phyto=bbp,
                    a_det=ad, b_det=bd, bb_det=bbd,
                    a_min=am, b_min=bm, bb_min=bbm)

    def evaluate_batch(self, wl_nm, chl, tsm, acdom440, t_c=20.0,
                       cdom_slope=CDOM_SLOPE_DEFAULT, det_a440=0.0,
                       det_slope=ORGANIC_DETRITUS_SLOPE_DEFAULT):
        """Batched evaluate over the case axis.  wl_nm scalar (shared within a
        band); chl, tsm, acdom440 are (B,) per-case.  Returns a dict of (B,)
        backend arrays.  All wavelength-only terms (star_ab, bb/b ratios, pure
        water) are computed once on host; only the constituent-dependent
        formulas carry the case axis.  Bit-identical to evaluate() applied case
        by case for grids where chl,tsm>0 (the mixed branch)."""
        from .backend import xp
        wl_nm = require_wavelength(wl_nm, context='OCRT constituent batch')
        chl = xp.asarray(chl, dtype=float)
        tsm = xp.asarray(tsm, dtype=float)
        acd = xp.asarray(acdom440, dtype=float)
        B = chl.shape[0]
        # --- wavelength-only scalars (host, computed once) ---
        aw, bw, bbw = pure_water(self.aw_lut, self.psi_lut, wl_nm, t_c)
        cdom_fac = float(np.exp(-cdom_slope * (wl_nm - 440.0)))
        astar_p = self.phyto_abs_lut.eval(wl_nm)
        r_d = bb_b_ratio(self.det_mie, wl_nm); r_d550 = bb_b_ratio(self.det_mie, 550.0)
        det_fac = float(np.exp(-det_slope * (wl_nm - 440.0)))
        astar_m, bstar_m = _star_ab(self.min_mie, wl_nm)
        r_m = bb_b_ratio(self.min_mie, wl_nm)
        # --- per-case (vectorized over B) ---
        pos_chl = chl > 0.0; pos_tsm = tsm > 0.0; pos_acd = acd > 0.0
        chl_safe = xp.where(pos_chl, chl, 1.0)
        fph = 0.035 + 0.015 * xp.tanh(xp.log10(chl_safe))
        huot550 = 2.267e-3 * chl_safe ** 0.565
        a_cdom = xp.maximum(0.0, xp.where(pos_acd, acd * cdom_fac, 0.0))
        a_ph = xp.where(pos_chl, astar_p * chl, 0.0)
        # Stage-2: species-specific phytoplankton scattering disabled.
        bb_ph = xp.zeros_like(chl, dtype=float)
        b_ph = xp.zeros_like(chl, dtype=float)
        a_det = xp.where(pos_chl, det_a440 * det_fac, 0.0)
        bb_det = xp.where(pos_chl, (1.0 - fph) * huot550 * (r_d / r_d550), 0.0)
        b_det = xp.where(pos_chl, bb_det / r_d, 0.0)
        a_min = xp.where(pos_tsm, astar_m * tsm, 0.0)
        b_min = xp.where(pos_tsm, bstar_m * tsm, 0.0)
        bb_min = xp.where(pos_tsm, b_min * r_m, 0.0)
        aw_a = xp.full(B, aw); bw_a = xp.full(B, bw); bbw_a = xp.full(B, bbw)
        a = aw_a + a_cdom + a_ph + a_det + a_min
        b = bw_a + b_ph + b_det + b_min
        bb = bbw_a + bb_ph + bb_det + bb_min
        return dict(a=a, b=b, bb=bb, a_w=aw_a, b_w=bw_a, bb_w=bbw_a, a_cdom=a_cdom,
                    a_phyto=a_ph, b_phyto=b_ph, bb_phyto=bb_ph,
                    a_det=a_det, b_det=b_det, bb_det=bb_det,
                    a_min=a_min, b_min=b_min, bb_min=bb_min)
