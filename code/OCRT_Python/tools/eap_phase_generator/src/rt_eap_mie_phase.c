/* rt_eap_mie_phase.c
 *
 * Public OCRT adapter: normalized 3-Stokes scattering-plane phase matrix for
 * one EAP phytoplankton species.  See rt_eap_mie_phase.h for the contract.
 *
 * The routine produces P11/P12/P33 only.  Absorption, scattering and
 * backscattering coefficients, bb/b, generalized moments and file output all
 * stay in their existing OCRT modules.
 */
#include <complex.h>
#include "rt_eap_mie_phase.h"
#include "internal/rt_eap_species_catalog.h"
#include "internal/rt_eap_coated_mie.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Frozen numerical settings.  Both are convergence results, not preferences:
 * the size step is the value at which the phase matrix stops changing (see
 * the validation package), and the upper size bound covers the
 * cross-section-weighted tail of the distribution. */
#define RT_EAP_PSD_STEP_UM   0.05
#define RT_EAP_PSD_MAX_MULT  5.5
#define RT_EAP_V_EFF         0.6
#define RT_EAP_N_MEDIA       1.334
#define RT_EAP_ENDPOINT_TOL  1.0e-12

static int validate_request(rt_eap_species_id_t species_id,
                            const double *wavelength_nm, size_t n_wavelength,
                            const double *theta_deg, size_t n_theta,
                            double *p11, double *p12, double *p33)
{
    if (!wavelength_nm || !theta_deg || !p11 || !p12 || !p33)
        return RT_EAP_PHASE_EINVAL;
    if (n_wavelength == 0 || n_theta < 3)
        return RT_EAP_PHASE_EINVAL;
    if (n_theta > (size_t)-1 / n_wavelength)
        return RT_EAP_PHASE_EINVAL;
    if (species_id >= RT_EAP_SPECIES_COUNT)
        return RT_EAP_PHASE_ESPECIES;

    if (fabs(theta_deg[0]) > RT_EAP_ENDPOINT_TOL) return RT_EAP_PHASE_ERANGE;
    if (fabs(theta_deg[n_theta - 1] - 180.0) > RT_EAP_ENDPOINT_TOL)
        return RT_EAP_PHASE_ERANGE;
    for (size_t i = 0; i < n_theta; ++i) {
        if (!isfinite(theta_deg[i])) return RT_EAP_PHASE_EINVAL;
        if (i > 0 && !(theta_deg[i] > theta_deg[i - 1]))
            return RT_EAP_PHASE_EINVAL;
    }
    const double lo = rt_eap_catalog_wl_lo(), hi = rt_eap_catalog_wl_hi();
    for (size_t i = 0; i < n_wavelength; ++i) {
        if (!isfinite(wavelength_nm[i]) || wavelength_nm[i] <= 0.0)
            return RT_EAP_PHASE_ERANGE;
        if (wavelength_nm[i] < lo || wavelength_nm[i] > hi)
            return RT_EAP_PHASE_ERANGE;
    }
    return RT_EAP_PHASE_OK;
}

/* Deirmendjian-form assemblage size distribution used by the EAP model.
 * Constant prefactors cancel in the phase-matrix normalization. */
static double psd_weight(double d_um, double deff_um)
{
    const double v = RT_EAP_V_EFF;
    const double r = 0.5 * d_um;
    return pow(r, (1.0 - 3.0 * v) / v) * exp(-r / (0.5 * deff_um * v));
}

int rt_eap_mie_phase_compute(rt_eap_species_id_t species_id,
                             const double *wavelength_nm, size_t n_wavelength,
                             const double *theta_deg, size_t n_theta,
                             double *p11, double *p12, double *p33)
{
    int rc = validate_request(species_id, wavelength_nm, n_wavelength,
                              theta_deg, n_theta, p11, p12, p33);
    if (rc != RT_EAP_PHASE_OK) return rc;

    const rt_eap_species_entry_t *sp = rt_eap_catalog_entry(species_id);
    if (!sp) return RT_EAP_PHASE_ESPECIES;
    const double deff = *sp->deff_um;
    const double vs   = *sp->vs;
    const double fr   = pow(1.0 - vs, 1.0 / 3.0);

    /* Fold the caller's angles onto the forward hemisphere; the Mie routine
     * returns each mu together with its supplement. */
    double *mu  = malloc(n_theta * sizeof *mu);
    int    *col = malloc(n_theta * sizeof *col);
    double *snt = malloc(n_theta * sizeof *snt);
    if (!mu || !col || !snt) { free(mu); free(col); free(snt);
                               return RT_EAP_PHASE_ENUMERIC; }
    for (size_t j = 0; j < n_theta; ++j) {
        double th = theta_deg[j];
        if (th > 90.0) { mu[j] = cos((180.0 - th) * M_PI / 180.0); col[j] = 1; }
        else           { mu[j] = cos(th * M_PI / 180.0);           col[j] = 0; }
        if (mu[j] < 0.0) mu[j] = 0.0;
        if (mu[j] > 1.0) mu[j] = 1.0;
        snt[j] = sin(th * M_PI / 180.0);
    }

    const size_t nd = (size_t)floor(RT_EAP_PSD_MAX_MULT * deff /
                                    RT_EAP_PSD_STEP_UM + 0.5);
    const size_t nout = n_wavelength * n_theta;
    double *a11 = calloc(nout, sizeof *a11);
    double *a12 = calloc(nout, sizeof *a12);
    double *a33 = calloc(nout, sizeof *a33);
    double *m1 = malloc(2 * n_theta * sizeof *m1);
    double *m2 = malloc(2 * n_theta * sizeof *m2);
    double *s21 = malloc(2 * n_theta * sizeof *s21);
    double *d21 = malloc(2 * n_theta * sizeof *d21);
    if (!a11 || !a12 || !a33 || !m1 || !m2 || !s21 || !d21 || nd < 4) {
        free(mu); free(col); free(snt); free(a11); free(a12); free(a33);
        free(m1); free(m2); free(s21); free(d21);
        return RT_EAP_PHASE_ENUMERIC;
    }

    rc = RT_EAP_PHASE_OK;
    for (size_t iw = 0; iw < n_wavelength && rc == RT_EAP_PHASE_OK; ++iw) {
        double nsh, ksh, nco, kco;
        if (rt_eap_catalog_ri(species_id, wavelength_nm[iw],
                              &nsh, &ksh, &nco, &kco) != 0) {
            rc = RT_EAP_PHASE_ERANGE; break;
        }
        /* Wave number in the medium; the catalog indices are relative to water. */
        const double wvno = 2.0 * M_PI * RT_EAP_N_MEDIA / (wavelength_nm[iw] * 1.0e-3);
        double *o11 = a11 + iw * n_theta;
        double *o12 = a12 + iw * n_theta;
        double *o33 = a33 + iw * n_theta;

        for (size_t k = 1; k <= nd; ++k) {
            const double d = (double)k * RT_EAP_PSD_STEP_UM;
            const double w = psd_weight(d, deff);
            if (!(w > 0.0) || !isfinite(w)) continue;
            rt_eap_mie_eff_t eff;
            int mrc = rt_eap_coated_mie(0.5 * d * fr, 0.5 * d, wvno,
                                        nsh - ksh * _Complex_I,
                                        nco - kco * _Complex_I,
                                        mu, n_theta, m1, m2, s21, d21, &eff);
            if (mrc != RT_EAP_MIE_OK) { rc = RT_EAP_PHASE_ENUMERIC; break; }
            for (size_t j = 0; j < n_theta; ++j) {
                const size_t q = j * 2 + (size_t)col[j];
                o11[j] += w * 0.5 * (m1[q] + m2[q]);
                o12[j] += w * 0.5 * (m2[q] - m1[q]);
                o33[j] += w * s21[q];
            }
        }
    }

    /* One common normalization per wavelength, evaluated with the same
     * quadrature the caller's grid defines, so the residual is exact. */
    if (rc == RT_EAP_PHASE_OK) {
        for (size_t iw = 0; iw < n_wavelength; ++iw) {
            double *o11 = a11 + iw * n_theta;
            double s = 0.0;
            for (size_t j = 0; j + 1 < n_theta; ++j) {
                const double dth = (theta_deg[j + 1] - theta_deg[j]) * M_PI / 180.0;
                s += 0.5 * dth * (o11[j] * snt[j] + o11[j + 1] * snt[j + 1]);
            }
            s *= 0.5;
            if (!(s > 0.0) || !isfinite(s)) { rc = RT_EAP_PHASE_ENUMERIC; break; }
            const double f = 1.0 / s;
            for (size_t j = 0; j < n_theta; ++j) {
                a11[iw * n_theta + j] *= f;
                a12[iw * n_theta + j] *= f;
                a33[iw * n_theta + j] *= f;
            }
        }
    }

    /* Physical admissibility. */
    if (rc == RT_EAP_PHASE_OK) {
        for (size_t iw = 0; iw < n_wavelength && rc == RT_EAP_PHASE_OK; ++iw) {
            double pmax = 0.0;
            for (size_t j = 0; j < n_theta; ++j)
                if (a11[iw * n_theta + j] > pmax) pmax = a11[iw * n_theta + j];
            const double eps = 1.0e-12 * pmax;
            for (size_t j = 0; j < n_theta; ++j) {
                const size_t q = iw * n_theta + j;
                if (!isfinite(a11[q]) || !isfinite(a12[q]) || !isfinite(a33[q]) ||
                    a11[q] < -eps ||
                    fabs(a12[q]) > a11[q] + eps ||
                    fabs(a33[q]) > a11[q] + eps) { rc = RT_EAP_PHASE_ENUMERIC; break; }
            }
        }
    }

    if (rc == RT_EAP_PHASE_OK) {
        memcpy(p11, a11, nout * sizeof *p11);
        memcpy(p12, a12, nout * sizeof *p12);
        memcpy(p33, a33, nout * sizeof *p33);
    }
    free(mu); free(col); free(snt); free(a11); free(a12); free(a33);
    free(m1); free(m2); free(s21); free(d21);
    return rc;
}
