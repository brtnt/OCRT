"""fixed-bulk phase LUT machinery (rt_water_rt.c transliteration)."""
import numpy as np
from math import pi, log, exp, cos, sin, log10


class PhaseTable:
    def __init__(self, theta_deg, p11):
        self.theta = np.asarray(theta_deg, dtype=float)
        self.p11 = np.asarray(p11, dtype=float)
        self.n = len(self.theta)


def phase_lut_load(path, case_id=None, wl_nm=0.0):
    """fixed_bulk_phase_lut_load: CSV with theta_deg,P11 (+optional case_id,
    wavelength_nm filters). Rows sorted by theta ascending."""
    with open(path, encoding='utf-8', errors='ignore') as f:
        header = [c.strip() for c in f.readline().split(',')]
        def col(*names):
            for nm in names:
                if nm in header:
                    return header.index(nm)
            return -1
        it = col('theta_deg', 'theta')
        ip = col('P11', 'p11')
        ic = col('case_id')
        iw = col('wavelength_nm')
        if it < 0 or ip < 0:
            raise ValueError('phase LUT header missing theta/P11')
        pairs = []
        for line in f:
            if not line.strip():
                continue
            fld = [c.strip() for c in line.split(',')]
            if len(fld) <= max(it, ip):
                continue
            if ic >= 0 and case_id:
                if len(fld) <= ic or fld[ic] != case_id:
                    continue
            if iw >= 0 and wl_nm > 0.0:
                if len(fld) <= iw or abs(float(fld[iw]) - wl_nm) > 1e-6:
                    continue
            t = float(fld[it])
            p = float(fld[ip])
            if not (0.0 <= t <= 180.0) or not (p >= 0.0):
                continue
            pairs.append((t, p))
    if len(pairs) < 2:
        raise ValueError('phase LUT: <2 rows after filter')
    pairs.sort(key=lambda x: x[0])
    return PhaseTable([p[0] for p in pairs], [p[1] for p in pairs])


def phase_lut_eval(tab, mu):
    """Scalar eval: log-theta/log-P11 interpolation with endpoint clamps."""
    mu = min(1.0, max(-1.0, mu))
    theta = np.arccos(mu) * 180.0 / pi
    if theta <= tab.theta[0]:
        return tab.p11[0]
    if theta >= tab.theta[-1]:
        return tab.p11[-1]
    hi = int(np.searchsorted(tab.theta, theta, side='right'))
    lo = hi - 1
    if tab.theta[lo] == theta:
        # C bisection lands lo at the equal node; interp with f computed below
        pass
    t0, t1 = tab.theta[lo], tab.theta[hi]
    p0, p1 = tab.p11[lo], tab.p11[hi]
    if p0 > 0.0 and p1 > 0.0:
        if theta > 0.0 and t0 > 0.0 and t1 > 0.0:
            f = (log(theta) - log(t0)) / (log(t1) - log(t0))
        else:
            f = (theta - t0) / (t1 - t0) if t1 > t0 else 0.0
        f = min(1.0, max(0.0, f))
        return exp(log(p0) * (1.0 - f) + log(p1) * f)
    f = (theta - t0) / (t1 - t0) if t1 > t0 else 0.0
    return p0 * (1.0 - f) + p1 * f


def phase_lut_eval_vec(tab, mu):
    """Vectorized eval (same math as phase_lut_eval)."""
    mu = np.clip(np.asarray(mu, dtype=float), -1.0, 1.0)
    theta = np.degrees(np.arccos(mu))
    out = np.empty_like(theta)
    lo_mask = theta <= tab.theta[0]
    hi_mask = theta >= tab.theta[-1]
    out[lo_mask] = tab.p11[0]
    out[hi_mask] = tab.p11[-1]
    mid = ~(lo_mask | hi_mask)
    th = theta[mid]
    hi = np.searchsorted(tab.theta, th, side='right')
    lo = hi - 1
    t0 = tab.theta[lo]
    t1 = tab.theta[hi]
    p0 = tab.p11[lo]
    p1 = tab.p11[hi]
    pos = (p0 > 0.0) & (p1 > 0.0)
    logok = pos & (th > 0.0) & (t0 > 0.0) & (t1 > 0.0)
    f = np.where(t1 > t0, (th - t0) / np.where(t1 > t0, t1 - t0, 1.0), 0.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        flog = (np.log(th) - np.log(np.where(t0 > 0, t0, 1.0))) / \
               (np.log(np.where(t1 > 0, t1, 1.0)) - np.log(np.where(t0 > 0, t0, 1.0)))
    f = np.where(logok, flog, f)
    f = np.clip(f, 0.0, 1.0)
    with np.errstate(divide='ignore'):
        vlog = np.exp(np.log(np.where(p0 > 0, p0, 1.0)) * (1.0 - f) +
                      np.log(np.where(p1 > 0, p1, 1.0)) * f)
    vlin = p0 * (1.0 - f) + p1 * f
    out[mid] = np.where(pos, vlog, vlin)
    return out


def phase_table_integral_theta(tab, lo_deg, hi_deg):
    """(1/2)∫ P11 sinθ dθ on the table grid, trapezoid per segment with
    log-log endpoint evaluation (matches C)."""
    if tab.n < 2 or hi_deg <= lo_deg:
        return 0.0
    acc = 0.0
    for i in range(tab.n - 1):
        t0 = tab.theta[i]
        t1 = tab.theta[i + 1]
        if t1 <= lo_deg or t0 >= hi_deg:
            continue
        a = max(t0, lo_deg)
        b = min(t1, hi_deg)
        if b <= a:
            continue
        f0 = (a - t0) / (t1 - t0) if t1 > t0 else 0.0
        f1 = (b - t0) / (t1 - t0) if t1 > t0 else 0.0
        p0 = tab.p11[i]
        p1 = tab.p11[i + 1]
        if p0 > 0.0 and p1 > 0.0:
            pa = exp(log(p0) * (1.0 - f0) + log(p1) * f0)
            pb = exp(log(p0) * (1.0 - f1) + log(p1) * f1)
        else:
            pa = p0 * (1.0 - f0) + p1 * f0
            pb = p0 * (1.0 - f1) + p1 * f1
        ar = a * pi / 180.0
        br = b * pi / 180.0
        acc += 0.25 * (pa * sin(ar) + pb * sin(br)) * (br - ar)
    return acc


def phase_table_osoaa_cap(tab, T1_deg, T2_deg):
    """fixed_bulk_phase_table_osoaa_cap: log-linear forward cap below T2 using
    the slope between T1/T2, then renormalize total=1. Returns removed mass A.
    Modifies tab in place."""
    if tab.n < 2 or not (T1_deg > T2_deg > 0.0):
        raise ValueError('cap args')
    total0 = phase_table_integral_theta(tab, 0.0, 180.0)
    if not total0 > 0.0:
        raise ValueError('cap total0')
    T1r = T1_deg * pi / 180.0
    T2r = T2_deg * pi / 180.0
    P1 = phase_lut_eval(tab, cos(T1r))
    P2 = phase_lut_eval(tab, cos(T2r))
    if not (P1 > 0.0 and P2 > 0.0):
        raise ValueError('cap P1/P2')
    AA = (log10(P2) - log10(P1)) / (T2r - T1r)
    for i in range(tab.n):
        if tab.theta[i] < T2_deg:
            thr = tab.theta[i] * pi / 180.0
            tab.p11[i] = 10.0 ** (log10(P2) + AA * (thr - T2r))
    total1 = phase_table_integral_theta(tab, 0.0, 180.0)
    if not total1 > 0.0:
        raise ValueError('cap total1')
    A = 1.0 - total1 / total0
    A = min(1.0 - 1.0e-6, max(0.0, A))
    tab.p11 /= total1
    return A


def direct_phase_fourier_allm(rm, n_mu, m_count, tab, nphi):
    """fixed_bulk_direct_phase_fourier_allm: value-kernel azimuth Fourier
    coefficients for all m.  rm: signed direction array indexable rm[j],
    j in [-n_mu..n_mu] (pass a dict-like or use offset array).
    Returns out[m][j][k] with j=0..n_mu, k index offset +n_mu (kw=2n_mu+1)."""
    nphi = max(16, min(20000, nphi))
    kw = 2 * n_mu + 1
    mu_j = np.array([rm[j] for j in range(0, n_mu + 1)])            # (J,)
    mu_k = np.array([rm[k] for k in range(-n_mu, n_mu + 1)])        # (K,)
    sj = np.sqrt(np.maximum(0.0, 1.0 - mu_j * mu_j))
    sk = np.sqrt(np.maximum(0.0, 1.0 - mu_k * mu_k))
    phi = 2.0 * pi * (np.arange(nphi) + 0.5) / nphi                 # (P,)
    cosphi = np.cos(phi)
    # cth[j,k,p]
    cth = mu_j[:, None, None] * mu_k[None, :, None] + \
        sj[:, None, None] * sk[None, :, None] * cosphi[None, None, :]
    np.clip(cth, -1.0, 1.0, out=cth)
    p11 = phase_lut_eval_vec(tab, cth.ravel()).reshape(cth.shape)   # (J,K,P)
    cosm = np.cos(np.arange(m_count)[:, None] * phi[None, :])       # (M,P)
    out = np.einsum('jkp,mp->mjk', p11, cosm) / nphi
    # reorder to [m][j][k_offset]
    return out  # shape (m_count, n_mu+1, kw); k axis already -n_mu..n_mu order
