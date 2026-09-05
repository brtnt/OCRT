/* =============================================================================
 * rt_atm.c — Plane-parallel atmosphere discretization for the OCRT SOS solver
 * =============================================================================
 *
 * VERTICAL GRID
 *   h[k] is cumulative vertical optical depth from TOA to level k.  The
 *   transport grid is uniform in optical depth; z_km_level[k] is the physical
 *   altitude associated with that cumulative column.  Molecular scattering
 *   follows the U.S. Standard Atmosphere 1962 pressure-column profile, while
 *   aerosol follows the configured finite exponential profile (production
 *   scale height: 2 km).  Absorbing-gas columns are subsequently integrated
 *   from the AFGL profile over these physical layer boundaries.
 *
 * DIRECT-BEAM NORMALIZATION
 *   ch[k] = 0.5 exp(-h[k]/mu_sun).  OCRT's Fourier SOS formulation places the
 *   complementary factor 1/2 in the analytical layer-source integral, giving
 *   the standard 1/4 solid-angle source normalization without carrying a
 *   delta-function solar direction in the quadrature weights.
 *
 * DIRECTION GRID
 *   rm[j], j in [-n_mu,+n_mu], stores mirrored positive Gauss-Legendre nodes.
 *   Slot j=0 stores the downward solar direction and has zero quadrature
 *   weight; it is therefore available to phase/source builders but cannot
 *   contribute as a diffuse incoming stream.
 *
 * MOLECULAR PHASE
 *   The depolarized Rayleigh Mueller matrix has only l=0 and l=2 moments in
 *   the OCRT 3-Stokes generalized spherical-function basis.  Their independent
 *   derivation and the m=0 l=2 basis are implemented in rt_molecular_phase.c.
 * ========================================================================== */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "rt_atm.h"
#include "rt_quadrature.h"
#include "rt_molecular_profile.h"
#include "rt_molecular_phase.h"
#include <string.h>

/* Normalized column fraction above altitude for a finite exponential aerosol
 * profile n(z) proportional to exp(-z/H) on [0,z_top]. */
static double aerosol_exp_fraction_above(double z_km, double z_top_km,
                                         double scale_height_km)
{
    if (z_km <= 0.0) return 1.0;
    if (z_km >= z_top_km) return 0.0;
    const double e_top = exp(-z_top_km / scale_height_km);
    return (exp(-z_km / scale_height_km) - e_top) / (1.0 - e_top);
}

int rt_atm_alloc(rt_atm_t *atm, int n_layers, int n_mu) {
    if (!atm || n_layers < 1 || n_mu < 1) return -1;
    atm->n_layers = n_layers;
    atm->n_mu     = n_mu;
    int dirs = 2 * n_mu + 1;

    atm->h           = (double *)calloc((size_t)(n_layers + 1), sizeof(double));
    atm->ch          = (double *)calloc((size_t)(n_layers + 1), sizeof(double));
    atm->xdel        = (double *)calloc((size_t)(n_layers + 1), sizeof(double));
    atm->ydel        = (double *)calloc((size_t)(n_layers + 1), sizeof(double));
    atm->z_km_level  = (double *)calloc((size_t)(n_layers + 1), sizeof(double));
    atm->tau_abs_layer = (double *)calloc((size_t)n_layers, sizeof(double));
    atm->tau_abs_total = 0.0;
    atm->rm_storage  = (double *)calloc((size_t)dirs, sizeof(double));
    atm->gb_storage  = (double *)calloc((size_t)dirs, sizeof(double));
    atm->xpl_storage = (double *)calloc((size_t)dirs, sizeof(double));
    atm->xrl_storage = (double *)calloc((size_t)dirs, sizeof(double));
    atm->xtl_storage = (double *)calloc((size_t)dirs, sizeof(double));

    if (!atm->h || !atm->ch || !atm->xdel || !atm->ydel ||
        !atm->z_km_level || !atm->tau_abs_layer ||
        !atm->rm_storage || !atm->gb_storage || !atm->xpl_storage ||
        !atm->xrl_storage || !atm->xtl_storage) {
        rt_atm_free(atm);
        return -1;
    }

    /* Offset pointers: rm[-n_mu..+n_mu] addressable. */
    atm->rm  = atm->rm_storage  + n_mu;
    atm->gb  = atm->gb_storage  + n_mu;
    atm->xpl = atm->xpl_storage + n_mu;
    atm->xrl = atm->xrl_storage + n_mu;
    atm->xtl = atm->xtl_storage + n_mu;

    /* Phase-1-irrelevant defaults that nevertheless make a blank atm sane. */
    atm->tau_total      = 0.0;
    atm->depol          = 0.0;
    atm->ssa            = 1.0;
    atm->mu_sun         = 0.0;
    atm->beta0          = 1.0;
    atm->beta2          = 0.0;
    atm->gamma2         = 0.0;
    atm->alpha2         = 0.0;
    atm->rayleigh_model = RT_RAYLEIGH_MODEL_BODHAINE_1999;

    /* Phase 3: aerosol coefficients (inactive at alloc) */
    atm->aerosol_active = 0;
    atm->L_max          = 0;
    atm->betal_aer      = NULL;
    atm->gammal_aer     = NULL;
    atm->alphal_aer     = NULL;
    atm->zetal_aer      = NULL;

    atm->transport_T = NULL;
    atm->transport_h_key = NULL;
    atm->transport_rm_key = NULL;
    atm->transport_cache_ready = 0;
    atm->transport_cache_builds = 0;
    atm->transport_cache_hits = 0;

    return 0;
}

void rt_atm_free(rt_atm_t *atm) {
    if (!atm) return;
    free(atm->h);            atm->h    = NULL;
    free(atm->ch);           atm->ch   = NULL;
    free(atm->xdel);         atm->xdel = NULL;
    free(atm->ydel);         atm->ydel = NULL;
    free(atm->z_km_level);   atm->z_km_level   = NULL;
    free(atm->tau_abs_layer); atm->tau_abs_layer = NULL;
    atm->tau_abs_total = 0.0;
    free(atm->rm_storage);   atm->rm_storage  = NULL; atm->rm  = NULL;
    free(atm->gb_storage);   atm->gb_storage  = NULL; atm->gb  = NULL;
    free(atm->xpl_storage);  atm->xpl_storage = NULL; atm->xpl = NULL;
    free(atm->xrl_storage);  atm->xrl_storage = NULL; atm->xrl = NULL;
    free(atm->xtl_storage);  atm->xtl_storage = NULL; atm->xtl = NULL;
    free(atm->transport_T);      atm->transport_T = NULL;
    free(atm->transport_h_key);  atm->transport_h_key = NULL;
    free(atm->transport_rm_key); atm->transport_rm_key = NULL;
    atm->transport_cache_ready = 0;
    atm->transport_cache_builds = 0;
    atm->transport_cache_hits = 0;

    /* Phase 3 aerosol coefficients (NULL-safe free) */
    free(atm->betal_aer);   atm->betal_aer  = NULL;
    free(atm->gammal_aer);  atm->gammal_aer = NULL;
    free(atm->alphal_aer);  atm->alphal_aer = NULL;
    free(atm->zetal_aer);   atm->zetal_aer  = NULL;
    atm->aerosol_active = 0;
    atm->L_max          = 0;
}

int rt_atm_build_extinction_grid_us62(rt_atm_t *atm, double mu_sun)
{
    if (!atm || !atm->h || !atm->ch || !atm->xdel || !atm->ydel ||
        !atm->z_km_level || !atm->tau_abs_layer) return -1;
    if (!(mu_sun > 0.0 && mu_sun <= 1.0)) return -1;

    const int nt = atm->n_layers;
    atm->tau_total = 0.0;
    atm->tau_abs_total = 0.0;
    atm->depol = 0.0;
    atm->ssa = 0.0;
    atm->mu_sun = mu_sun;
    atm->rayleigh_model = RT_RAYLEIGH_MODEL_BODHAINE_1999;
    atm->beta0 = 0.0;
    atm->beta2 = 0.0;
    atm->gamma2 = 0.0;
    atm->alpha2 = 0.0;
    atm->aerosol_active = 0;

    for (int k = 0; k <= nt; ++k) {
        atm->h[k] = 0.0;
        atm->ch[k] = 0.5;
        atm->xdel[k] = 0.0;
        atm->ydel[k] = 0.0;
        atm->z_km_level[k] = (k == 0)
            ? 100.0
            : rt_molecular_us62_altitude_from_grid_fraction(
                  (double)k / (double)nt);
        if (k < nt) atm->tau_abs_layer[k] = 0.0;
    }
    return 0;
}

int rt_atm_build_rayleigh(rt_atm_t *atm, double tau_total,
                          double depol, double mu_sun,
                          rt_rayleigh_model_t rayleigh_model) {
    if (!atm) return -1;
    if (!atm->h || !atm->ch || !atm->rm_storage) return -1;   /* not allocated */
    if (!(tau_total > 0.0))  return -1;
    if (!(depol >= 0.0 && depol < 1.0)) return -1;
    if (!(mu_sun > 0.0 && mu_sun <= 1.0)) return -1;

    const int nt   = atm->n_layers;
    const int n_mu = atm->n_mu;

    atm->tau_total      = tau_total;
    atm->depol          = depol;
    atm->ssa            = 1.0;
    atm->mu_sun         = mu_sun;
    atm->rayleigh_model = rayleigh_model;

    /* Molecular phase moments are derived from the normalized depolarized
     * Rayleigh Mueller matrix.  Only l=0 and l=2 are nonzero. */
    rt_molecular_phase_coeffs_t molecular_phase;
    if (rt_molecular_phase_coefficients(depol, &molecular_phase) != 0)
        return -1;
    atm->beta0  = molecular_phase.beta0;
    atm->beta2  = molecular_phase.beta2;
    atm->gamma2 = molecular_phase.gamma2;
    atm->alpha2 = molecular_phase.alpha2;

    /* Layer arrays. h[k] = k * tau / nt; ch[k] = exp(-h/mu_sun)/2.
     * The optical-depth grid itself is unchanged.  Only the mapping from
     * cumulative molecular optical depth to geometric altitude uses the
     * documented US Standard Atmosphere 1962 pressure-column profile:
     *
     *   h[k] / tau_R = F_mol_above(z_k)
     *
     * where F_mol_above is normalized to one at sea level and zero at the
     * US62 model top.  k=0 retains the 100 km AFGL integration cap so that
     * absorbing-gas columns above the US62 anchor top remain represented in
     * the first OCRT layer.
     */
    const double z_cap = 100.0;
    for (int k = 0; k <= nt; k++) {
        double h_k = (double)k * tau_total / (double)nt;
        atm->h[k]    = h_k;
        atm->ch[k]   = 0.5 * exp(-h_k / mu_sun);
        atm->xdel[k] = 0.0;
        atm->ydel[k] = 1.0;
        if (k == 0) {
            atm->z_km_level[k] = z_cap;
        } else {
            const double f_above = h_k / tau_total;
            atm->z_km_level[k] =
                rt_molecular_us62_altitude_from_grid_fraction(f_above);
        }
        if (k < nt) atm->tau_abs_layer[k] = 0.0;  /* no gas absorption yet */
    }
    atm->tau_abs_total = 0.0;

    /* Gauss-Legendre quadrature → fill rm[±j], gb[±j].
     * mu_quad[i] is ascending; place at j = i+1, mirror at j = -(i+1).
     */
    double mu_quad[256], w_quad[256];
    if (n_mu > 256) return -1;     /* defensive; bump if ever needed */
    int rc = rt_quadrature_gauss_legendre_pos(n_mu, mu_quad, w_quad);
    if (rc != 0) return -2;

    /* solar slot */
    atm->rm[0]  = -mu_sun;
    atm->gb[0]  = 0.0;
    /* The solar slot carries the same m=0, l=2 basis as the discrete
     * directions, but has zero quadrature weight. */
    rt_molecular_phase_l2_m0(atm->rm[0],
                             &atm->xpl[0], &atm->xrl[0], &atm->xtl[0]);

    for (int j = 1; j <= n_mu; j++) {
        double m = mu_quad[j - 1];
        double w = w_quad[j - 1];
        atm->rm[+j] = +m;
        atm->rm[-j] = -m;
        atm->gb[+j] = w;
        atm->gb[-j] = w;

        double p2, r2, t2;
        rt_molecular_phase_l2_m0(m, &p2, &r2, &t2);
        atm->xpl[+j] = p2;
        atm->xpl[-j] = p2;
        atm->xrl[+j] = r2;
        atm->xrl[-j] = r2;
        atm->xtl[+j] = t2;
        atm->xtl[-j] = t2;
    }

    return 0;
}

int rt_atm_build_aerosol_rayleigh(
    rt_atm_t *atm,
    double tau_R, double depol,
    double tau_a, double ssa_a,
    int L_max,
    const double *betal_aer, const double *gammal_aer,
    const double *alphal_aer, const double *zetal_aer,
    double mu_sun,
    rt_rayleigh_model_t rayleigh_model,
    double aer_h_km)
{
    if (!atm) return -1;
    if (!atm->h || !atm->ch || !atm->rm_storage) return -1;   /* not allocated */
    if (!(tau_R >= 0.0 && tau_a >= 0.0 && tau_R + tau_a > 0.0)) return -1;
    if (!(ssa_a >= 0.0 && ssa_a <= 1.0)) return -1;
    if (L_max < 2) return -1;
    if (!betal_aer || !gammal_aer || !alphal_aer || !zetal_aer) return -1;
    if (!(depol >= 0.0 && depol < 1.0)) return -1;
    if (!(mu_sun > 0.0 && mu_sun <= 1.0)) return -1;

    /* First step: build the Rayleigh atmosphere with τ_total = τ_R + τ_a.
     * This populates h[k], ch[k], rm/gb/xpl/xrl/xtl correctly (geometry
     * and quadrature are unaffected by aerosol presence).  It also sets
     * xdel[k]=0, ydel[k]=1 and ssa=1 — we will overwrite those below.   */
    double tau_total = tau_R + tau_a;
    int rc = rt_atm_build_rayleigh(atm, tau_total, depol, mu_sun, rayleigh_model);
    if (rc != 0) return rc;

    /* Second step: distribute the fixed column optical depths over physical
     * altitude.  Rayleigh follows the U.S. Standard Atmosphere 1962
     * molecular column.  Aerosol follows a finite exponential profile; the
     * production model is H=2 km.  A positive advanced override remains
     * accepted, while non-positive legacy values now select 2 km rather than
     * the removed ODA550/an23 reference-code table. */
    const int nt = atm->n_layers;
    const double z_top = 100.0;
    const double aerosol_h_km = (aer_h_km > 0.0) ? aer_h_km : 2.0;

    if (tau_R == 0.0 && tau_a > 0.0) {
        atm->xdel[0] = ssa_a;
        atm->ydel[0] = 0.0;
    } else {
        atm->xdel[0] = 0.0;
        atm->ydel[0] = 1.0;
    }

    double z_prev = z_top;
    double ray_above_prev = 0.0;
    double aer_above_prev = 0.0;
    atm->z_km_level[0] = z_top;

    for (int j = 1; j <= nt; ++j) {
        const double target = atm->h[j];
        double z_j;
        if (j == nt) {
            z_j = 0.0;
        } else {
            double z_lo = 0.0;
            double z_hi = z_prev;
            for (int it = 0; it < 64; ++it) {
                const double z_mid = 0.5 * (z_lo + z_hi);
                const double tau_ray_above =
                    tau_R * rt_molecular_us62_grid_fraction_above(z_mid);
                const double tau_aer_above =
                    tau_a * aerosol_exp_fraction_above(z_mid, z_top,
                                                       aerosol_h_km);
                const double tau_above = tau_ray_above + tau_aer_above;
                if (tau_above > target) z_lo = z_mid;
                else                    z_hi = z_mid;
            }
            z_j = 0.5 * (z_lo + z_hi);
        }

        atm->z_km_level[j] = z_j;
        const double ray_above_j =
            tau_R * rt_molecular_us62_grid_fraction_above(z_j);
        const double aer_above_j =
            tau_a * aerosol_exp_fraction_above(z_j, z_top, aerosol_h_km);
        const double dt_ray = ray_above_j - ray_above_prev;
        const double dt_aer = aer_above_j - aer_above_prev;
        const double dt = dt_ray + dt_aer;

        if (dt > 0.0) {
            atm->xdel[j] = dt_aer * ssa_a / dt;
            atm->ydel[j] = dt_ray / dt;
        } else {
            atm->xdel[j] = 0.0;
            atm->ydel[j] = (tau_R > 0.0) ? 1.0 : 0.0;
        }

        z_prev = z_j;
        ray_above_prev = ray_above_j;
        aer_above_prev = aer_above_j;
    }

    for (int j = 0; j < nt; ++j) atm->tau_abs_layer[j] = 0.0;
    atm->tau_abs_total = 0.0;

    /* Third step: store mixed single-scattering albedo as a metadata. */
    atm->ssa = (tau_R * 1.0 + tau_a * ssa_a) / tau_total;

    /* Fourth step: allocate + copy aerosol Legendre coefficients. */
    if (atm->betal_aer)  free(atm->betal_aer);
    if (atm->gammal_aer) free(atm->gammal_aer);
    if (atm->alphal_aer) free(atm->alphal_aer);
    if (atm->zetal_aer)  free(atm->zetal_aer);
    size_t bytes = (size_t)(L_max + 1) * sizeof(double);
    atm->betal_aer  = (double *)malloc(bytes);
    atm->gammal_aer = (double *)malloc(bytes);
    atm->alphal_aer = (double *)malloc(bytes);
    atm->zetal_aer  = (double *)malloc(bytes);
    if (!atm->betal_aer || !atm->gammal_aer || !atm->alphal_aer || !atm->zetal_aer) {
        free(atm->betal_aer);  atm->betal_aer  = NULL;
        free(atm->gammal_aer); atm->gammal_aer = NULL;
        free(atm->alphal_aer); atm->alphal_aer = NULL;
        free(atm->zetal_aer);  atm->zetal_aer  = NULL;
        atm->aerosol_active = 0;
        atm->L_max          = 0;
        return -2;
    }
    memcpy(atm->betal_aer,  betal_aer,  bytes);
    memcpy(atm->gammal_aer, gammal_aer, bytes);
    memcpy(atm->alphal_aer, alphal_aer, bytes);
    memcpy(atm->zetal_aer,  zetal_aer,  bytes);
    atm->aerosol_active = 1;
    atm->L_max          = L_max;

    return 0;
}

void rt_atm_dump(const rt_atm_t *atm, FILE *fp) {
    if (!atm || !fp) return;
    const char *model_name =
        (atm->rayleigh_model == RT_RAYLEIGH_MODEL_BODHAINE_1999)
            ? "BODHAINE_1999" : "UNKNOWN";
    fprintf(fp, "rt_atm_t dump\n");
    fprintf(fp, "  n_layers=%d  n_mu=%d  rayleigh_model=%s\n",
            atm->n_layers, atm->n_mu, model_name);
    fprintf(fp, "  tau_total=%.10g  depol=%.6g  ssa=%.6g  mu_sun=%.10g\n",
            atm->tau_total, atm->depol, atm->ssa, atm->mu_sun);
    fprintf(fp, "  beta0=%.10g  beta2=%.10g  gamma2=%.10g  alpha2=%.10g\n",
            atm->beta0, atm->beta2, atm->gamma2, atm->alpha2);

    fprintf(fp, "  k         h[k]              ch[k]              xdel[k]   ydel[k]\n");
    for (int k = 0; k <= atm->n_layers; k++) {
        fprintf(fp, "  %-3d  %18.10e  %18.10e  %.4f    %.4f\n",
                k, atm->h[k], atm->ch[k], atm->xdel[k], atm->ydel[k]);
    }

    fprintf(fp, "  j      rm[j]               gb[j]               xpl[j]              xrl[j]              xtl[j]\n");
    for (int j = -atm->n_mu; j <= atm->n_mu; j++) {
        fprintf(fp, "  %+4d  %18.10e  %18.10e  %18.10e  %18.10e  %18.10e\n",
                j, atm->rm[j], atm->gb[j], atm->xpl[j], atm->xrl[j], atm->xtl[j]);
    }
}

/* === v1.02 Gas absorption integration (SOS-coupled) =========================
 *
 * Forward-declared in rt_atm.h as taking `struct rt_absorption_state *`;
 * the actual type is `rt_absorption_t` from rt_absorption.h.
 * Since rt_atm.h and rt_absorption.h have no shared types, we just include
 * rt_absorption.h here and provide the implementation.
 * ===========================================================================*/
#include "rt_absorption.h"

int rt_atm_apply_gas_absorption(rt_atm_t *atm,
                                 const struct rt_absorption_state *abs_opaque,
                                 double wavelength_nm)
{
    /* Cast opaque to concrete type. */
    const rt_absorption_t *st = (const rt_absorption_t *)abs_opaque;
    if (!atm || !st || !st->initialized) return -1;
    if (!atm->z_km_level || !atm->tau_abs_layer || !atm->h || !atm->ch) return -1;

    const int nt = atm->n_layers;
    const double mu_sun = atm->mu_sun;
    if (!(mu_sun > 0.0)) return -1;

    /* Step 1: compute per-layer dt_abs from AFGL + xsec. */
    double dt_abs_total = 0.0;
    for (int k = 0; k < nt; k++) {
        /* OCRT layer k spans z_km_level[k] (upper, larger z) to z_km_level[k+1] (lower, smaller z).
         * For Rayleigh-only path, z_km_level[0] = 100 km (TOA cap), z_km_level[nt] = 0.
         * For aerosol path, same convention.
         */
        double z_hi = atm->z_km_level[k];
        double z_lo = atm->z_km_level[k+1];
        if (z_lo > z_hi) { double t = z_lo; z_lo = z_hi; z_hi = t; }

        /* Mid-layer altitude for xsec interpolation.
         * v1.02 (updated): Use BI-LINEAR interpolation in (z, λ) space.
         * Nearest-neighbor is FORBIDDEN — always interpolate between bracketing
         * AFGL z-levels so the xsec n_layers can differ from rt_atm nt freely.
         *
         * The xsec file's "layer index" maps 1:1 to AFGL z_km grid (st->atm.z_km).
         * That mapping is enforced by the production script
         * (scripts/generate_xsec_per_gas.py) which uses the AFGL profile to
         * select per-layer T_k, P_k for HAPI calls.
         */
        double z_mid = 0.5 * (z_lo + z_hi);

        double dt_abs_k = 0.0;
        for (int g = 0; g < RT_N_GAS; g++) {
            /* Per-gas column for this OCRT layer band */
            double N_g_layer;
            if (st->column_overridden[g]) {
                /* User-overridden column — distribute proportionally to AFGL default per layer.
                 * Use ratio: N_user / N_default_total × N_default_layer. */
                double N_def_total = st->atm.column_default[g];
                if (N_def_total > 0.0) {
                    double N_def_layer = rt_afgl_layer_column(&st->atm, g, z_lo, z_hi);
                    N_g_layer = st->column_eff[g] * (N_def_layer / N_def_total);
                } else {
                    N_g_layer = 0.0;
                }
            } else {
                N_g_layer = rt_afgl_layer_column(&st->atm, g, z_lo, z_hi);
            }
            /* BI-LINEAR interpolation: σ(z_mid, λ) from xsec table + AFGL z_km grid.
             * If xsec has only 1 layer (legacy), falls back to layer-0 σ(λ). */
            double sigma_g = rt_xsec_interp_layer_z(&st->xsec[g],
                                                     st->atm.z_km, z_mid,
                                                     wavelength_nm);
            dt_abs_k += sigma_g * N_g_layer;
        }

        atm->tau_abs_layer[k] = dt_abs_k;
        dt_abs_total += dt_abs_k;
    }
    atm->tau_abs_total = dt_abs_total;

    /* Step 2: rebuild h[k], ch[k] cumulating gas absorption into total OD.
     *
     * Convention: h[0] = 0 (TOA), h[k+1] = h[k] + dt_total_new[k].
     * dt_total_old[k] = h_old[k+1] - h_old[k]  (Rayleigh + aerosol scat+abs)
     * dt_total_new[k] = dt_total_old[k] + dt_abs_layer[k]
     */
    double cum_abs = 0.0;
    /* Capture original h[k] before rewriting */
    double *h_old = (double*)malloc((size_t)(nt+1) * sizeof(double));
    if (!h_old) return -1;
    for (int k = 0; k <= nt; k++) h_old[k] = atm->h[k];

    for (int k = 0; k <= nt; k++) {
        atm->h[k] = h_old[k] + cum_abs;
        atm->ch[k] = 0.5 * exp(-atm->h[k] / mu_sun);
        if (k < nt) cum_abs += atm->tau_abs_layer[k];
    }

    /* Step 3: re-normalize xdel[k], ydel[k] for k=1..nt.
     * xdel[k] × dt_old[k] = dt_aer_scat[k]  (absolute aerosol scattering OD)
     * ydel[k] × dt_old[k] = dt_ray[k]       (absolute Rayleigh scattering OD)
     * In Rayleigh-only build path, xdel[k]=0, ydel[k]=1, dt_old[k] = h_old[k]-h_old[k-1].
     * In aerosol+Rayleigh path, xdel/ydel already set.
     *
     * After absorption: dt_new[k] = dt_old[k] + dt_abs_layer[k-1] (note 0-indexed layer below k).
     * Scattering ODs unchanged, only denominator changes.
     */
    for (int k = 1; k <= nt; k++) {
        double dt_old = h_old[k] - h_old[k-1];
        double dt_abs_k = atm->tau_abs_layer[k-1];
        double dt_new   = dt_old + dt_abs_k;
        if (dt_new > 0.0 && dt_old > 0.0) {
            atm->xdel[k] = atm->xdel[k] * dt_old / dt_new;
            atm->ydel[k] = atm->ydel[k] * dt_old / dt_new;
        }
        /* else: degenerate, leave as-is */
    }
    /* k=0: xdel[0], ydel[0] are pure-Rayleigh or pure-aerosol top-of-atm
     * "Phase 0" values, NOT layer-integrated. Treatment matches build code. */

    /* Update tau_total to include absorption */
    atm->tau_total = atm->h[nt];

    free(h_old);
    return 0;
}
