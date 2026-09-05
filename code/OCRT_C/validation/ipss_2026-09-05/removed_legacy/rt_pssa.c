/* =============================================================================
 * rt_pssa.c — Atmospheric Pseudo-Spherical Approximation
 * =============================================================================
 * See rt_pssa.h for the full physics/geometry documentation. Everything here
 * is per-(case, wavelength) atmospheric precomputation; nothing runs inside
 * SOS loops.
 *
 * Explicit scope boundary: this module does not correct the refracted
 * in-water beam.  OCRT intentionally keeps the underwater direct beam
 * plane-parallel; --pssa ends at 0+ (the air side of the interface).
 * ========================================================================== */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdatomic.h>
#include "rt_pssa.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* --------------------------------------------------------------------------
 * F(z) for the DOWNWARD beam toward target level at altitude z_k:
 *   F(z) = sqrt((R+z)^2 - D^2),  D = (R+z_k) sin(theta0)
 * evaluated cancellation-free (exact algebraic identity, see rt_pssa.h):
 *   F(z)^2 = ((R+z) mu0)^2 + (z - z_k)(2R + z + z_k) s2,  s2 = sin^2(theta0)
 * Valid on the downward path where z >= z_k (both terms >= 0); a max(.,0)
 * guard covers any non-monotone z profile degeneracy.
 * ------------------------------------------------------------------------ */
static double pssa_F_level(double R, double z, double z_k,
                           double mu0, double s2)
{
    const double a = (R + z) * mu0;
    double v = a * a + (z - z_k) * (2.0 * R + z + z_k) * s2;
    return sqrt(v > 0.0 ? v : 0.0);
}

/* Same F, for a ray whose impact parameter is p = R sin(alpha) (surface
 * grazing point at z = 0):  F0(z)^2 = ((R+z) cos(alpha))^2 + z(2R+z) sin^2a.
 * Used for BOTH legs of the reflected beam (they share p). */
static double pssa_F_surface_ray(double R, double z,
                                 double cos_a, double s2a)
{
    const double a = (R + z) * cos_a;
    double v = a * a + z * (2.0 * R + z) * s2a;
    return sqrt(v > 0.0 ? v : 0.0);
}

/* Solve paper Eq. 8 for the local Fresnel incidence angle alpha at target
 * altitude z_k:  f(a) = (R+z_k) sin(2a - theta0) - R sin(a) = 0 on
 * (theta0/2, theta0].  f(theta0) = z_k sin(theta0) >= 0, f(theta0/2+) < 0,
 * f monotone increasing on the bracket -> plain bisection is safe. */
static double pssa_solve_alpha(double R, double z_k, double theta0)
{
    if (z_k <= 0.0) return theta0;              /* surface level: alpha = theta0 */
    if (sin(theta0) < 1.0e-12) return theta0;   /* sun ~overhead: alpha ~ 0     */

    double lo = 0.5 * theta0 + 1.0e-14;
    double hi = theta0;
    for (int it = 0; it < 100; ++it) {
        double mid = 0.5 * (lo + hi);
        double f = (R + z_k) * sin(2.0 * mid - theta0) - R * sin(mid);
        if (f < 0.0) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

/* High-SZA PSSA uses the physical altitude assigned to each optical-depth
 * layer.  With too few layers, a single upper-atmosphere layer can span a
 * large altitude interval while the shell secant changes rapidly.  Emit an
 * advisory once per severity level; this does not alter the solve or add work
 * inside any SOS loop.  Thresholds are conservative and are based on the
 * current 412-nm US62 convergence audit (40..1600 layers). */
static void pssa_warn_layer_resolution(int n_layers, double sza_deg)
{
    int level = 0;
    if (sza_deg >= 84.0 && n_layers < 400) level = 3;
    else if (sza_deg >= 80.0 && n_layers < 200) level = 2;
    else if (sza_deg >= 75.0 && n_layers < 100) level = 1;
    if (level == 0) return;

    static _Atomic unsigned warned_mask = 0u;
    const unsigned bit = 1u << (unsigned)level;
    const unsigned old = atomic_fetch_or_explicit(&warned_mask, bit,
                                                   memory_order_relaxed);
    if (old & bit) return;

    fprintf(stderr,
            "warning: PSSA vertical resolution may be insufficient at "
            "SZA=%.3f deg with n_layers=%d.\n",
            sza_deg, n_layers);
    if (level == 1) {
        fprintf(stderr,
                "         recommendation: use --n-layers >=100 for "
                "SZA 75-80 deg.\n");
    } else if (level == 2) {
        fprintf(stderr,
                "         recommendation: use --n-layers >=200 for "
                "SZA 80-84 deg.\n");
    } else {
        fprintf(stderr,
                "         near 85 deg, use >=200 layers for conservative "
                "TOA-radiance convergence and 400-800 layers when the "
                "PSSA correction magnitude itself must be sub-percent.\n");
    }
    fprintf(stderr,
            "         non-default Rayleigh-only layer counts require "
            "OCRT_ADVANCED=1.\n");
}

int rt_pssa_apply(rt_atm_t *atm)
{
    if (!atm || !atm->h || !atm->z_km_level) return -1;
    const int nt = atm->n_layers;
    if (nt < 1) return -1;
    const double mu0 = atm->mu_sun;
    if (!(mu0 > 0.0 && mu0 <= 1.0)) return -1;

    const double R  = RT_PSSA_EARTH_RADIUS_KM;
    const double s2 = 1.0 - mu0 * mu0;                 /* sin^2(theta0)   */
    const double theta0 = acos(mu0);
    pssa_warn_layer_resolution(nt, theta0 * (180.0 / M_PI));
    const double *z = atm->z_km_level;                 /* z[0]=TOA .. z[nt]=0 */
    const double *h = atm->h;                          /* FINAL tau grid  */

    /* Lazy allocation (sizes fixed by n_layers set at rt_atm_alloc). */
    if (!atm->pssa_xi_dn)
        atm->pssa_xi_dn   = (double *)malloc((size_t)(nt + 1) * sizeof(double));
    if (!atm->pssa_xi_refl)
        atm->pssa_xi_refl = (double *)malloc((size_t)(nt + 1) * sizeof(double));
    if (!atm->pssa_alpha)
        atm->pssa_alpha   = (double *)malloc((size_t)(nt + 1) * sizeof(double));
    if (!atm->pssa_beta)
        atm->pssa_beta    = (double *)malloc((size_t)(nt + 1) * sizeof(double));
    double *Flev = (double *)malloc((size_t)(nt + 1) * sizeof(double));
    if (!atm->pssa_xi_dn || !atm->pssa_xi_refl || !atm->pssa_alpha ||
        !atm->pssa_beta || !Flev) {
        free(Flev);
        return -1;   /* pssa arrays (if any) are freed by rt_atm_free() */
    }

    /* ---- (a) DOWNWARD beam: xi_dn[k], paper Eqs. 6-7 -------------------- */
    for (int k = 0; k <= nt; ++k) {
        const double zk = z[k];
        for (int i = 0; i <= k; ++i)
            Flev[i] = pssa_F_level(R, z[i], zk, mu0, s2);
        double xi = 0.0;
        for (int i = 0; i < k; ++i) {
            const double dtau = h[i + 1] - h[i];
            /* sec_i(k) = S_i / dz_i = (2R + z_i + z_{i+1}) / (F_i + F_{i+1})
             * — exact (F_i^2 - F_{i+1}^2 = (z_i - z_{i+1})(2R+z_i+z_{i+1})),
             * finite even for dz_i = 0 (levels merged at the z-cap): the
             * local shell secant is used, which is the physical limit. */
            const double sec_i = (2.0 * R + z[i] + z[i + 1])
                                 / (Flev[i] + Flev[i + 1]);
            xi += dtau * sec_i;
        }
        atm->pssa_xi_dn[k] = xi;   /* xi_dn[0] = 0 -> ch[0] = 0.5 as before */
    }

    /* ---- (b) REFLECTED beam: alpha_k and xi_refl[k], paper Eq. 8 -------- */
    for (int k = 0; k <= nt; ++k) {
        const double alpha = pssa_solve_alpha(R, z[k], theta0);
        double beta = 2.0 * alpha - theta0;
        if (beta < 0.0) beta = 0.0;
        if (beta > 0.5 * M_PI) beta = 0.5 * M_PI;
        atm->pssa_alpha[k] = alpha;
        atm->pssa_beta[k] = beta;

        const double cos_a = cos(alpha);
        const double s2a   = 1.0 - cos_a * cos_a;
        for (int i = 0; i <= nt; ++i)
            Flev[i] = pssa_F_surface_ray(R, z[i], cos_a, s2a);

        /* Down leg TOA->surface at incidence alpha (all layers) plus up leg
         * surface->level k along the reflected ray (layers i >= k); both
         * legs share the impact parameter p = R sin(alpha), hence one F set:
         *   xi_refl[k] = sum_{i<k} dtau_i sec0_i + 2 * sum_{i>=k} dtau_i sec0_i */
        double xi = 0.0;
        for (int i = 0; i < nt; ++i) {
            const double dtau  = h[i + 1] - h[i];
            const double sec_i = (2.0 * R + z[i] + z[i + 1])
                                 / (Flev[i] + Flev[i + 1]);
            xi += (i >= k ? 2.0 : 1.0) * dtau * sec_i;
        }
        atm->pssa_xi_refl[k] = xi;   /* k = nt: up leg empty -> = xi_dn[nt] */
    }
    free(Flev);

    /* ---- Site A: overwrite the primary-source attenuation --------------- */
    for (int k = 0; k <= nt; ++k)
        atm->ch[k] = 0.5 * exp(-atm->pssa_xi_dn[k]);

    atm->pssa_active = 1;
    return 0;
}
