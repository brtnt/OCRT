"""rt_absorption.c transliteration: AFGL profiles + xsec LUTs -> tau_abs."""
import numpy as np
import os

from .spectral_contract import (require_wavelength, require_table_coverage,
                                require_query_in_table)

from .backend import xp

GAS_NAMES = ["h2o", "o3", "o2", "co2", "no2", "ch4"]
N_GAS = 6


class AfglAtm:
    pass


def afgl_load(afgl_dir, profile_name="usstd76"):
    path = os.path.join(afgl_dir, f"afgl_{profile_name}.dat")
    z, P, T, nt = [], [], [], []
    mr = [[] for _ in range(N_GAS)]
    with open(path, encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            v = line.split()
            if len(v) != 10:
                continue
            z.append(float(v[0])); P.append(float(v[1]))
            T.append(float(v[2])); nt.append(float(v[3]))
            for g in range(N_GAS):
                mr[g].append(float(v[4 + g]))
    atm = AfglAtm()
    atm.z_km = np.array(z)
    atm.P_mbar = np.array(P)
    atm.T_K = np.array(T)
    atm.n_total = np.array(nt)
    atm.mr = [np.array(m) for m in mr]
    atm.n_levels = len(z)
    atm.column_default = np.zeros(N_GAS)
    for g in range(N_GAS):
        col = 0.0
        for i in range(atm.n_levels - 1):
            n1 = atm.mr[g][i] * 1e-6 * atm.n_total[i]
            n2 = atm.mr[g][i + 1] * 1e-6 * atm.n_total[i + 1]
            dz_cm = (atm.z_km[i + 1] - atm.z_km[i]) * 1.0e5
            col += 0.5 * (n1 + n2) * dz_cm
        atm.column_default[g] = col
    return atm


class Xsec:
    def __init__(self):
        self.available = False
        self.n_layers = 0
        self.wl = None
        self.sigma = None  # (n_layers, n_wl)


def xsec_load(xsec_dir, gas_idx):
    xs = Xsec()
    path = os.path.join(xsec_dir, f"xsec_{GAS_NAMES[gas_idx]}.dat")
    if not os.path.exists(path):
        return xs
    n_layers = 1
    rows = []
    with open(path, encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('#'):
                if 'n_layers' in line:
                    try:
                        n_layers = int(line.split('n_layers =')[1].split()[0])
                    except Exception:
                        pass
                continue
            if not line.strip():
                continue
            rows.append(line.split())
    if n_layers == 1:
        wl = [float(r[0]) for r in rows]
        sg = [float(r[1]) for r in rows]
        xs.n_layers = 1
        xs.wl = require_table_coverage(np.array(wl),
                                       context=f'xsec {GAS_NAMES[gas_idx]}')
        xs.sigma = np.array(sg)[None, :]
        xs.available = True
        return xs
    # layered: row0 = wavelengths, rows 1..K = per-layer sigma
    xs.n_layers = n_layers
    xs.wl = require_table_coverage(
        np.array([float(v) for v in rows[0]]),
        context=f'xsec {GAS_NAMES[gas_idx]}')
    sig = np.zeros((n_layers, len(xs.wl)))
    for L in range(n_layers):
        if L + 1 < len(rows):
            vals = [float(v) for v in rows[L + 1]]
            sig[L, :len(vals)] = vals
    xs.sigma = sig
    xs.available = True
    return xs


def xsec_interp_layer(xs, layer_idx, wl_nm):
    if not xs.available or len(xs.wl) < 2:
        return 0.0
    wl_nm = require_wavelength(wl_nm, context='gas cross section')
    require_query_in_table(wl_nm, xs.wl, context='gas cross section')
    layer_idx = min(max(layer_idx, 0), xs.n_layers - 1)
    row = xs.sigma[layer_idx]
    if wl_nm == xs.wl[0]:
        return row[0]
    if wl_nm == xs.wl[-1]:
        return row[-1]
    hi = int(np.searchsorted(xs.wl, wl_nm, side='right'))
    lo = hi - 1
    t = (wl_nm - xs.wl[lo]) / (xs.wl[hi] - xs.wl[lo])
    return (1.0 - t) * row[lo] + t * row[hi]


def xsec_interp_layer_z(xs, z_grid, z_km, wl_nm):
    if not xs.available:
        return 0.0
    n = xs.n_layers
    if n == 1:
        return xsec_interp_layer(xs, 0, wl_nm)
    if z_km <= z_grid[0]:
        return xsec_interp_layer(xs, 0, wl_nm)
    if z_km >= z_grid[n - 1]:
        return xsec_interp_layer(xs, n - 1, wl_nm)
    hi = int(np.searchsorted(z_grid[:n], z_km, side='right'))
    lo = hi - 1
    dz = z_grid[hi] - z_grid[lo]
    if dz <= 0.0:
        return xsec_interp_layer(xs, lo, wl_nm)
    t = (z_km - z_grid[lo]) / dz
    return (1.0 - t) * xsec_interp_layer(xs, lo, wl_nm) + \
        t * xsec_interp_layer(xs, hi, wl_nm)


def afgl_cumulative_column(atm, g, z_km):
    cum = 0.0
    for i in range(atm.n_levels - 1):
        z1 = atm.z_km[i]
        z2 = atm.z_km[i + 1]
        if z2 <= z_km:
            continue
        zlo = z_km if z1 < z_km else z1
        zhi = z2
        if zhi <= zlo:
            continue
        f1 = (zlo - z1) / (z2 - z1)
        f2 = (zhi - z1) / (z2 - z1)
        n1 = atm.mr[g][i] * 1e-6 * atm.n_total[i]
        n2 = atm.mr[g][i + 1] * 1e-6 * atm.n_total[i + 1]
        nlo = (1.0 - f1) * n1 + f1 * n2
        nhi = (1.0 - f2) * n1 + f2 * n2
        cum += 0.5 * (nlo + nhi) * (zhi - zlo) * 1.0e5
    return cum


def afgl_layer_column(atm, g, z_lo, z_hi):
    if z_hi <= z_lo:
        return 0.0
    return afgl_cumulative_column(atm, g, z_lo) - afgl_cumulative_column(atm, g, z_hi)


def afgl_cumulative_column_vec(atm, g, z_arr):
    """Vectorized afgl_cumulative_column: molecular column above each z in
    z_arr (M,).  Trapezoid over AFGL levels with partial top layer.  Returns
    (M,) array."""
    z_arr = xp.asarray(z_arr, dtype=float)
    z_afgl = xp.asarray(atm.z_km, dtype=float)                    # (nlev,)
    n_g = (xp.asarray(atm.mr[g], dtype=float) * 1e-6
           * xp.asarray(atm.n_total, dtype=float))                # (nlev,)
    z1 = z_afgl[:-1][None, :]                                     # (1,nlev-1)
    z2 = z_afgl[1:][None, :]
    n1 = n_g[:-1][None, :]
    n2 = n_g[1:][None, :]
    z = z_arr[:, None]                                            # (M,1)
    zlo = xp.maximum(z, z1)
    zhi = z2
    dzl = z2 - z1
    f1 = (zlo - z1) / dzl
    f2 = (zhi - z1) / dzl
    nlo = (1.0 - f1) * n1 + f1 * n2
    nhi = (1.0 - f2) * n1 + f2 * n2
    contrib = 0.5 * (nlo + nhi) * (zhi - zlo) * 1.0e5
    valid = (z2 > z) & (zhi > zlo)
    contrib = xp.where(valid, contrib, 0.0)
    return contrib.sum(axis=1)                                   # (M,)


def xsec_layers_at_wl(xs, wl_nm):
    """xsec at wl for all layers: (n_layers,).  Linear in wavelength."""
    if not xs.available or len(xs.wl) < 2:
        return xp.zeros(xs.n_layers)
    sig = xp.asarray(xs.sigma, dtype=float)                       # (n_layers,n_wl)
    if wl_nm <= xs.wl[0]:
        return sig[:, 0]
    if wl_nm >= xs.wl[-1]:
        return sig[:, -1]
    hi = int(np.searchsorted(xs.wl, wl_nm, side='right'))
    lo = hi - 1
    t = (wl_nm - xs.wl[lo]) / (xs.wl[hi] - xs.wl[lo])
    return (1.0 - t) * sig[:, lo] + t * sig[:, hi]


def xsec_interp_layer_z_vec(xs, z_grid, z_arr, wl_nm):
    """Vectorized xsec_interp_layer_z: xsec at each z_mid in z_arr (M,) and
    wavelength wl_nm.  Returns (M,) array."""
    z_arr = xp.asarray(z_arr, dtype=float)
    if not xs.available:
        return xp.zeros_like(z_arr)
    n = xs.n_layers
    sig = xsec_layers_at_wl(xs, wl_nm)                            # (n_layers,)
    if n == 1:
        return xp.full_like(z_arr, float(sig[0]))
    zg = xp.asarray(z_grid[:n], dtype=float)
    hi = xp.clip(xp.searchsorted(zg, z_arr, side='right'), 1, n - 1)
    lo = hi - 1
    dz = zg[hi] - zg[lo]
    t = xp.where(dz > 0.0, (z_arr - zg[lo]) / dz, 0.0)
    res = (1.0 - t) * sig[lo] + t * sig[hi]
    res = xp.where(z_arr <= zg[0], sig[0], res)
    res = xp.where(z_arr >= zg[n - 1], sig[n - 1], res)
    return res


class Absorption:
    def __init__(self, afgl_dir, xsec_dir, profile="usstd76", overrides=None):
        self.atm = afgl_load(afgl_dir, profile)
        self.xsec = [xsec_load(xsec_dir, g) for g in range(N_GAS)]
        self.column_eff = self.atm.column_default.copy()
        self.overridden = [False] * N_GAS
        if overrides:
            for g in range(N_GAS):
                if overrides[g] >= 0.0:
                    self.column_eff[g] = overrides[g]
                    self.overridden[g] = True

    def tau_per_gas(self, g, wl_nm):
        xs = self.xsec[g]
        if not xs.available:
            return 0.0
        atm = self.atm
        Ndef = atm.column_default[g]
        scale = (self.column_eff[g] / Ndef) if (self.overridden[g] and Ndef > 0.0) else 1.0
        tau = 0.0
        for k in range(atm.n_levels - 1):
            z_lo = atm.z_km[k]
            z_hi = atm.z_km[k + 1]
            z_mid = 0.5 * (z_lo + z_hi)
            sigma = xsec_interp_layer_z(xs, atm.z_km, z_mid, wl_nm)
            N = afgl_layer_column(atm, g, z_lo, z_hi)
            if self.overridden[g] and Ndef > 0.0:
                N *= scale
            tau += sigma * N
        return tau

    def tau_total(self, wl_nm):
        return sum(self.tau_per_gas(g, wl_nm) for g in range(N_GAS))
