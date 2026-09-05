"""Direct angle-space scalar/vector phase kernels.

This module mirrors ``src/rt_value_phase.c``.  It is intentionally separate
from the Legendre/moment caches: enabling the value kernel cannot overwrite or
poison the coefficient path.
"""
from math import pi, cos, sin, sqrt
import numpy as np

from .kernel import gauss_legendre_pos
from .surface import build_rotation_L, mat3_mul


class ValuePhaseInterp:
    def __init__(self, theta_deg, p11, p12=None, p33=None, norm_n_mu=400):
        th = np.asarray(theta_deg, float)
        p11 = np.asarray(p11, float)
        if p12 is None:
            p12 = np.zeros_like(p11)
        if p33 is None:
            p33 = p11.copy()
        rows = sorted(zip(np.cos(np.deg2rad(th)), p11,
                          np.asarray(p12, float), np.asarray(p33, float)),
                      key=lambda r: r[0])
        uniq = []
        for r in rows:
            if uniq and abs(r[0] - uniq[-1][0]) < 1.0e-14:
                uniq[-1] = r
            else:
                uniq.append(r)
        if len(uniq) < 3:
            raise ValueError('value phase table needs >=3 unique mu nodes')
        a = np.asarray(uniq, float)
        self.mu = a[:, 0]
        self.p11 = a[:, 1]
        self.p12 = a[:, 2]
        self.p33 = a[:, 3]
        self.y2_11 = _spline_build(self.mu, self.p11)
        self.y2_12 = _spline_build(self.mu, self.p12)
        self.y2_33 = _spline_build(self.mu, self.p33)
        n = max(2, int(norm_n_mu or 400))
        x, w = gauss_legendre_pos(n)
        norm = 0.5 * np.sum(w * (self.eval_p11(x) + self.eval_p11(-x)))
        if not np.isfinite(norm) or not norm > 0.0:
            raise ValueError('invalid value-phase normalization')
        inv = 1.0 / float(norm)
        self.p11 *= inv; self.p12 *= inv; self.p33 *= inv
        self.y2_11 = _spline_build(self.mu, self.p11)
        self.y2_12 = _spline_build(self.mu, self.p12)
        self.y2_33 = _spline_build(self.mu, self.p33)

    def eval_p11(self, mu):
        return _spline_eval_vec(self.mu, self.p11, self.y2_11, mu)

    def eval(self, mu):
        return (_spline_eval_vec(self.mu, self.p11, self.y2_11, mu),
                _spline_eval_vec(self.mu, self.p12, self.y2_12, mu),
                _spline_eval_vec(self.mu, self.p33, self.y2_33, mu))


def _spline_build(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    n = len(x); y2 = np.zeros(n); u = np.zeros(n - 1)
    dx0 = x[1] - x[0]; dxn = x[-1] - x[-2]
    yp1 = (y[1] - y[0]) / dx0
    ypn = (y[-1] - y[-2]) / dxn
    y2[0] = -0.5
    u[0] = (3.0 / dx0) * ((y[1] - y[0]) / dx0 - yp1)
    for i in range(1, n - 1):
        dxm = x[i] - x[i - 1]; dxp = x[i + 1] - x[i]
        sig = dxm / (x[i + 1] - x[i - 1])
        p = sig * y2[i - 1] + 2.0
        y2[i] = (sig - 1.0) / p
        dd = (y[i + 1] - y[i]) / dxp - (y[i] - y[i - 1]) / dxm
        u[i] = (6.0 * dd / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p
    qn = 0.5
    un = (3.0 / dxn) * (ypn - (y[-1] - y[-2]) / dxn)
    y2[-1] = (un - qn * u[-1]) / (qn * y2[-2] + 1.0)
    for k in range(n - 2, -1, -1):
        y2[k] = y2[k] * y2[k + 1] + u[k]
    return y2


def _spline_eval_vec(x, y, y2, q):
    scalar = np.ndim(q) == 0
    q = np.asarray(q, float)
    qc = np.clip(q, x[0], x[-1])
    hi = np.searchsorted(x, qc, side='right')
    hi = np.clip(hi, 1, len(x) - 1)
    lo = hi - 1
    h = x[hi] - x[lo]
    a = (x[hi] - qc) / h
    b = (qc - x[lo]) / h
    out = a*y[lo] + b*y[hi] + ((a**3-a)*y2[lo] + (b**3-b)*y2[hi])*h*h/6.0
    out = np.where(q <= x[0], y[0], np.where(q >= x[-1], y[-1], out))
    return float(out) if scalar else out


def scalar_fourier_allm(rm, n_mu, m_count, phase, nphi):
    nphi = max(16, min(20000, int(nphi)))
    kw = 2*n_mu + 1
    out = np.zeros((m_count, n_mu+1, kw))
    phis = 2*pi*(np.arange(nphi)+0.5)/nphi
    cm = np.cos(np.arange(m_count)[:, None]*phis[None, :])
    cp = np.cos(phis)
    for j in range(n_mu+1):
        mui = rm[j+n_mu]; si = sqrt(max(0.0, 1.0-mui*mui))
        for ko, k in enumerate(range(-n_mu, n_mu+1)):
            muo = rm[k+n_mu]; so = sqrt(max(0.0, 1.0-muo*muo))
            cth = np.clip(mui*muo + si*so*cp, -1.0, 1.0)
            out[:, j, ko] = cm @ phase.eval_p11(cth) / nphi
    return out


def _hovenier(muo, mui, so, si, cT, sT, sph):
    d1 = max(1e-12, si*sT); d2 = max(1e-12, so*sT); sts=max(1e-12,sT)
    c1=(muo-mui*cT)/d1; s1r=so*sph/sts
    n1=sqrt(c1*c1+s1r*s1r)
    ci1,si1=(1.0,0.0) if n1<1e-12 else (c1/n1,s1r/n1)
    c2=(mui-muo*cT)/d2; s2r=si*sph/sts
    n2=sqrt(c2*c2+s2r*s2r)
    ci2,si2=(1.0,0.0) if n2<1e-12 else (c2/n2,s2r/n2)
    return ci1,si1,ci2,si2


def vector_fourier_allm(rm, n_mu, m_count, phase, nphi):
    """Return pfm,gr,gt,arr,art,att, each [m,j,k+n_mu]."""
    nphi=max(16,min(20000,int(nphi))); kw=2*n_mu+1
    out=[np.zeros((m_count,n_mu+1,kw)) for _ in range(6)]
    phis=2*pi*(np.arange(nphi)+0.5)/nphi
    for j in range(n_mu+1):
        mui=rm[j+n_mu]; si=sqrt(max(0.0,1.0-mui*mui))
        for ko,k in enumerate(range(-n_mu,n_mu+1)):
            muo=rm[k+n_mu]; so=sqrt(max(0.0,1.0-muo*muo))
            ca=np.zeros((4,m_count)); sa=np.zeros((2,m_count))
            for phi in phis:
                cp=cos(phi); sp=sin(phi)
                cT=max(-1.0,min(1.0,mui*muo+si*so*cp)); sT=sqrt(max(0.0,1.0-cT*cT))
                p11,p12,p33=(float(v) for v in phase.eval(cT))
                P=np.array([p11,p12,0.0,p12,p11,0.0,0.0,0.0,p33])
                ci1,si1,ci2,si2=_hovenier(muo,mui,so,si,cT,sT,sp)
                Z=mat3_mul(build_rotation_L(ci2,si2),mat3_mul(P,build_rotation_L(ci1,si1)))
                ms=np.arange(m_count); c=np.cos(ms*phi); ss=np.sin(ms*phi)
                ca[0]+=Z[0]*c; ca[1]+=Z[3]*c; ca[2]+=Z[4]*c; ca[3]+=Z[8]*c
                sa[0]+=(-Z[6])*ss; sa[1]+=Z[5]*ss
            out[0][:,j,ko]=ca[0]/nphi; out[1][:,j,ko]=ca[1]/nphi
            out[2][:,j,ko]=sa[0]/nphi; out[3][:,j,ko]=ca[2]/nphi
            out[4][:,j,ko]=sa[1]/nphi; out[5][:,j,ko]=ca[3]/nphi
    return tuple(out)
