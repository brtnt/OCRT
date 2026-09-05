/* rt_eap_coated_mie.c
 *
 * Aden-Kerker stratified-sphere scattering, well-conditioned formulation of
 * Toon & Ackerman (Appl. Opt. 20, 3657, 1981).
 *
 * Differences from the legacy DMiLay Fortran routine that OCRT's EAP tables
 * were originally produced with:
 *
 *   1. The downward-recurrence start order is derived from the size
 *      parameters instead of being hard-wired.  DMiLay fixes NMX1 = 1500 and
 *      therefore diverges once |m|x exceeds roughly 1350, which corresponds to
 *      a 180 um sphere at 400 nm.  The EAP catalog needs Raphidophytes
 *      (D_eff = 60 um, distribution tail beyond 300 um), so the fixed start
 *      order is not usable.
 *   2. Workspace is heap-allocated per call, so there is no LL = 2000 ceiling
 *      and the routine is re-entrant.
 *   3. Series termination is by relative contribution, with a hard bound at
 *      the recurrence start order.
 *
 * The scattering-amplitude bookkeeping, the recurrences and the Mueller
 * assembly are a faithful transcription of the reference routine so that a
 * bit-level comparison against it remains possible wherever it is valid.
 */
#include "rt_eap_coated_mie.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Extra orders carried above the physical stop order before the downward
 * recurrence is started.  The reference routine used 1500 against a working
 * ceiling of about 1350, i.e. roughly 11%.  A fixed absolute margin is used
 * here in addition, so that small particles keep a generous head-room. */
#define RT_EAP_MIE_ORDER_MARGIN_REL  0.15
#define RT_EAP_MIE_ORDER_MARGIN_ABS  150
#define RT_EAP_MIE_ORDER_MIN         64
#define RT_EAP_MIE_TERM_TOL          1.0e-11

static int start_order(double zmax)
{
    /* Physical stop order plus margin. */
    double nstop = zmax + 4.0 * cbrt(zmax) + 2.0;
    double n = nstop * (1.0 + RT_EAP_MIE_ORDER_MARGIN_REL)
             + RT_EAP_MIE_ORDER_MARGIN_ABS;
    if (n < RT_EAP_MIE_ORDER_MIN) n = RT_EAP_MIE_ORDER_MIN;
    if (!(n < 2.0e6)) return -1;
    return (int)n;
}

int rt_eap_coated_mie(double rcore, double rshell, double wvno,
                      double _Complex rindsh, double _Complex rindco,
                      const double *mu, size_t numang,
                      double *m1, double *m2, double *s21, double *d21,
                      rt_eap_mie_eff_t *eff)
{
    if (!eff) return RT_EAP_MIE_EINVAL;
    const int want_phase = numang > 0;
    if (want_phase && (!mu || !m1 || !m2 || !s21 || !d21))
        return RT_EAP_MIE_EINVAL;
    if (!(wvno > 0.0) || !(rshell > 0.0) || !(rcore > 0.0) || rcore > rshell)
        return RT_EAP_MIE_EINVAL;
    if (!(creal(rindsh) > 0.0) || cimag(rindsh) > 0.0) return RT_EAP_MIE_EINVAL;
    if (!(creal(rindco) > 0.0) || cimag(rindco) > 0.0) return RT_EAP_MIE_EINVAL;
    if (want_phase) {
        for (size_t j = 0; j < numang; ++j) {
            if (!isfinite(mu[j]) || mu[j] < -1.0e-6 || mu[j] > 1.0 + 1.0e-6)
                return RT_EAP_MIE_EINVAL;
        }
    }

    const double xshell = rshell * wvno;
    const double xcore  = rcore * wvno;

    const double _Complex ci = I;
    const double _Complex k1 = rindco * wvno;
    const double _Complex k2 = rindsh * wvno;
    const double _Complex k3 = wvno + 0.0 * I;

    double _Complex z[5];
    z[1] = rindsh * xshell;
    z[2] = xshell + 0.0 * I;
    z[3] = rindco * xcore;
    z[4] = rindsh * xcore;

    double zmax = cabs(z[1]);
    for (int i = 2; i <= 4; ++i) if (cabs(z[i]) > zmax) zmax = cabs(z[i]);

    const int nmx1 = start_order(zmax);
    if (nmx1 < 0) return RT_EAP_MIE_ERANGE;

    double _Complex *acap = malloc((size_t)(nmx1 + 2) * sizeof *acap);
    double _Complex *wa   = malloc(3u * (size_t)(nmx1 + 2) * sizeof *wa);
    double *pi0 = want_phase ? malloc(numang * sizeof *pi0) : NULL;
    double *pi1 = want_phase ? malloc(numang * sizeof *pi1) : NULL;
    double *pi2 = want_phase ? malloc(numang * sizeof *pi2) : NULL;
    double *ta0 = want_phase ? malloc(numang * sizeof *ta0) : NULL;
    double *ta1 = want_phase ? malloc(numang * sizeof *ta1) : NULL;
    double *ta2 = want_phase ? malloc(numang * sizeof *ta2) : NULL;
    double *si2 = want_phase ? malloc(numang * sizeof *si2) : NULL;
    double _Complex *s1 = want_phase ? malloc(2u * numang * sizeof *s1) : NULL;
    double _Complex *s2 = want_phase ? malloc(2u * numang * sizeof *s2) : NULL;
    if (!acap || !wa || (want_phase && (!pi0 || !pi1 || !pi2 || !ta0 ||
        !ta1 || !ta2 || !si2 || !s1 || !s2))) {
        free(acap); free(wa); free(pi0); free(pi1); free(pi2);
        free(ta0); free(ta1); free(ta2); free(si2); free(s1); free(s2);
        return RT_EAP_MIE_ENOMEM;
    }
#define W(m, n) wa[(size_t)((m) - 1) * (size_t)(nmx1 + 2) + (size_t)(n)]

    /* Downward recurrences for the logarithmic derivatives. */
    acap[nmx1 + 1] = 0.0;
    for (int m = 1; m <= 3; ++m) W(m, nmx1 + 1) = 0.0;
    const double _Complex rrfx = 1.0 / (rindsh * xshell);
    for (int nn = nmx1; nn >= 1; --nn) {
        acap[nn] = ((double)(nn + 1) * rrfx)
                 - 1.0 / (((double)(nn + 1) * rrfx) + acap[nn + 1]);
        for (int m = 1; m <= 3; ++m) {
            double _Complex t = (double)(nn + 1) / z[m + 1];
            W(m, nn) = t - 1.0 / (t + W(m, nn + 1));
        }
    }

    if (want_phase) {
        for (size_t j = 0; j < numang; ++j) {
            si2[j] = 1.0 - mu[j] * mu[j];
            pi0[j] = 0.0; pi1[j] = 1.0;
            ta0[j] = 0.0; ta1[j] = mu[j];
        }
    }

    const double rx = 1.0 / xshell;
    double t1 = cos(xshell), t2 = sin(xshell);
    double _Complex wm1 = t1 - t2 * I;
    double _Complex wfn0 = t2 + t1 * I;
    double _Complex wfn1 = rx * wfn0 - wm1;
    double ta3 = creal(wfn1);

    const double x1 = creal(z[1]), y1 = cimag(z[1]);
    const double x4 = creal(z[4]), y4 = cimag(z[4]);
    const double sinx1 = sin(x1), cosx1 = cos(x1);
    const double sinx4 = sin(x4), cosx4 = cos(x4);
    const double ey1 = exp(y1), e2y1 = ey1 * ey1, ey4 = exp(y4);
    const double ey1my4 = exp(y1 - y4), ey1py4 = ey1 * ey4;
    const double aa = sinx4 * (ey1py4 + ey1my4);
    const double bb = cosx4 * (ey1py4 - ey1my4);
    const double cc = sinx1 * (e2y1 + 1.0);
    const double dd = cosx1 * (e2y1 - 1.0);
    const double denom = 1.0 + e2y1 * (4.0 * sinx1 * sinx1 - 2.0 + e2y1);
    if (!(fabs(denom) > 0.0) || !isfinite(denom)) {
        free(acap); free(wa); free(pi0); free(pi1); free(pi2);
        free(ta0); free(ta1); free(ta2); free(si2); free(s1); free(s2);
        return RT_EAP_MIE_ENUMERIC;
    }
    double _Complex dummy = ((aa * cc + bb * dd) / denom)
                          + ((bb * cc - aa * dd) / denom) * I;

    int n = 1;
    dummy = dummy * (acap[n] + (double)n / z[1]) / (W(3, n) + (double)n / z[4]);
    double _Complex dumsq = dummy * dummy;

    double _Complex p24h24 = 0.5
        + ((sinx4 * sinx4 - 0.5) + (cosx4 * sinx4) * I) * (ey4 * ey4);
    double _Complex p24h21 =
        0.5 * ((sinx1 * sinx4 - cosx1 * cosx4)
               + (sinx1 * cosx4 + cosx1 * sinx4) * I) * ey1py4
      + 0.5 * ((sinx1 * sinx4 + cosx1 * cosx4)
               + (-sinx1 * cosx4 + cosx1 * sinx4) * I) * ey1my4;

    double _Complex dh1 = z[1] / (1.0 + ci * z[1]) - 1.0 / z[1];
    double _Complex dh2 = z[2] / (1.0 + ci * z[2]) - 1.0 / z[2];
    double _Complex dh4 = z[4] / (1.0 + ci * z[4]) - 1.0 / z[4];

    p24h24 = p24h24 / ((dh4 + (double)n / z[4]) * (W(3, n) + (double)n / z[4]));
    p24h21 = p24h21 / ((dh1 + (double)n / z[1]) * (W(3, n) + (double)n / z[4]));

    double _Complex u[9];
    u[1] = k3 * acap[n] - k2 * W(1, n);
    u[2] = k3 * acap[n] - k2 * dh2;
    u[3] = k2 * acap[n] - k3 * W(1, n);
    u[4] = k2 * acap[n] - k3 * dh2;
    u[5] = k1 * W(3, n) - k2 * W(2, n);
    u[6] = k2 * W(3, n) - k1 * W(2, n);
    u[7] = -ci * (dummy * p24h21 - p24h24);
    u[8] = ta3 / wfn1;

    double _Complex acoe = u[8] * (u[1] * u[5] * u[7] + k1 * u[1] - dumsq * k3 * u[5])
                                / (u[2] * u[5] * u[7] + k1 * u[2] - dumsq * k3 * u[5]);
    double _Complex bcoe = u[8] * (u[3] * u[6] * u[7] + k2 * u[3] - dumsq * k2 * u[6])
                                / (u[4] * u[6] * u[7] + k2 * u[4] - dumsq * k2 * u[6]);
    double _Complex acoem1 = acoe, bcoem1 = bcoe;

    double dqext = 3.0 * (creal(acoe) + creal(bcoe));
    double dqsca = 3.0 * (creal(acoe) * creal(acoe) + cimag(acoe) * cimag(acoe)
                        + creal(bcoe) * creal(bcoe) + cimag(bcoe) * cimag(bcoe));
    double dgqsc = 0.0;
    double _Complex sback = 3.0 * (acoe - bcoe);
    double rmm = 1.0;

    double _Complex ac = 1.5 * acoe, bc = 1.5 * bcoe;
    if (want_phase) {
        for (size_t j = 0; j < numang; ++j) {
            s1[j * 2 + 0] = ac * pi1[j] + bc * ta1[j];
            s1[j * 2 + 1] = ac * pi1[j] - bc * ta1[j];
            s2[j * 2 + 0] = bc * pi1[j] + ac * ta1[j];
            s2[j * 2 + 1] = bc * pi1[j] - ac * ta1[j];
        }
    }

    double t4 = 0.0;
    for (n = 2; n <= nmx1; ++n) {
        const double c1 = (double)(2 * n - 1);
        const double c2 = (double)(n - 1);
        if (want_phase) {
            for (size_t j = 0; j < numang; ++j) {
                pi2[j] = (c1 * pi1[j] * mu[j] - (double)n * pi0[j]) / c2;
                ta2[j] = mu[j] * (pi2[j] - pi0[j]) - c1 * si2[j] * pi1[j] + ta0[j];
            }
        }
        wm1 = wfn0; wfn0 = wfn1;
        wfn1 = c1 * rx * wfn0 - wm1;
        ta3 = creal(wfn1);

        dh1 = -(double)n / z[1] + 1.0 / ((double)n / z[1] - dh1);
        dh2 = -(double)n / z[2] + 1.0 / ((double)n / z[2] - dh2);
        dh4 = -(double)n / z[4] + 1.0 / ((double)n / z[4] - dh4);
        p24h24 = p24h24 / ((dh4 + (double)n / z[4]) * (W(3, n) + (double)n / z[4]));
        p24h21 = p24h21 / ((dh1 + (double)n / z[1]) * (W(3, n) + (double)n / z[4]));
        dummy = dummy * (acap[n] + (double)n / z[1]) / (W(3, n) + (double)n / z[4]);
        dumsq = dummy * dummy;

        u[1] = k3 * acap[n] - k2 * W(1, n);
        u[2] = k3 * acap[n] - k2 * dh2;
        u[3] = k2 * acap[n] - k3 * W(1, n);
        u[4] = k2 * acap[n] - k3 * dh2;
        u[5] = k1 * W(3, n) - k2 * W(2, n);
        u[6] = k2 * W(3, n) - k1 * W(2, n);
        u[7] = -ci * (dummy * p24h21 - p24h24);
        u[8] = ta3 / wfn1;
        acoe = u[8] * (u[1] * u[5] * u[7] + k1 * u[1] - dumsq * k3 * u[5])
                    / (u[2] * u[5] * u[7] + k1 * u[2] - dumsq * k3 * u[5]);
        bcoe = u[8] * (u[3] * u[6] * u[7] + k2 * u[3] - dumsq * k2 * u[6])
                    / (u[4] * u[6] * u[7] + k2 * u[4] - dumsq * k2 * u[6]);

        const double are = creal(acoe), aim = cimag(acoe);
        const double bre = creal(bcoe), bim = cimag(bcoe);
        const double am1re = creal(acoem1), am1im = cimag(acoem1);
        const double bm1re = creal(bcoem1), bm1im = cimag(bcoem1);

        const double g4 = (2.0 * n - 1.0) / ((double)n * (n - 1.0));
        const double g2 = (n - 1.0) * (n + 1.0) / (double)n;
        dgqsc += g2 * (am1re * are + am1im * aim + bm1re * bre + bm1im * bim)
               + g4 * (am1re * bm1re + am1im * bm1im);

        const double c3 = (double)(2 * n + 1);
        dqext += c3 * (are + bre);
        t4 = are * are + aim * aim + bre * bre + bim * bim;
        dqsca += c3 * t4;
        rmm = -rmm;
        sback += c3 * rmm * (acoe - bcoe);

        const double f = c3 / ((double)n * (n + 1.0));
        ac = f * acoe; bc = f * bcoe;
        if (want_phase) {
            for (size_t j = 0; j < numang; ++j) {
                s1[j * 2 + 0] += ac * pi2[j] + bc * ta2[j];
                s2[j * 2 + 0] += bc * pi2[j] + ac * ta2[j];
            }
            if (n % 2 == 0) {
                for (size_t j = 0; j < numang; ++j) {
                    s1[j * 2 + 1] += -ac * pi2[j] + bc * ta2[j];
                    s2[j * 2 + 1] += -bc * pi2[j] + ac * ta2[j];
                }
            } else {
                for (size_t j = 0; j < numang; ++j) {
                    s1[j * 2 + 1] += ac * pi2[j] - bc * ta2[j];
                    s2[j * 2 + 1] += bc * pi2[j] - ac * ta2[j];
                }
            }
        }

        if (!(t4 >= RT_EAP_MIE_TERM_TOL) || !isfinite(t4)) break;

        if (want_phase) {
            double *sw;
            sw = pi0; pi0 = pi1; pi1 = pi2; pi2 = sw;
            sw = ta0; ta0 = ta1; ta1 = ta2; ta2 = sw;
        }
        acoem1 = acoe; bcoem1 = bcoe;
    }

    int rc = RT_EAP_MIE_OK;
    if (!isfinite(t4)) rc = RT_EAP_MIE_ENUMERIC;
    if (n >= nmx1 && t4 >= RT_EAP_MIE_TERM_TOL) rc = RT_EAP_MIE_ENOCONV;

    if (want_phase) {
        for (size_t j = 0; j < numang; ++j) {
            for (int k = 0; k < 2; ++k) {
                const double _Complex a = s1[j * 2 + k];
                const double _Complex b = s2[j * 2 + k];
                m1[j * 2 + k]  = creal(a) * creal(a) + cimag(a) * cimag(a);
                m2[j * 2 + k]  = creal(b) * creal(b) + cimag(b) * cimag(b);
                s21[j * 2 + k] = creal(a) * creal(b) + cimag(a) * cimag(b);
                d21[j * 2 + k] = cimag(a) * creal(b) - cimag(b) * creal(a);
            }
        }
    }
    const double sc = 2.0 * rx * rx;
    eff->qext = sc * dqext;
    eff->qsca = sc * dqsca;
    eff->gqsc = 2.0 * sc * dgqsc;
    eff->qbs  = rx * rx * (creal(sback) * creal(sback) + cimag(sback) * cimag(sback));
    eff->n_terms = n;

#undef W
    free(acap); free(wa); free(pi0); free(pi1); free(pi2);
    free(ta0); free(ta1); free(ta2); free(si2); free(s1); free(s2);
    return rc;
}
