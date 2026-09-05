/* ============================================================================
 * rt_ipss.c — IPSS geometry, spherical single scattering, Eq. (7) factor
 * Implementation notes and scope: see rt_ipss.h.
 * Oracle: IPSS_REFERENCE_PACKAGE_2026-08-25/ipss_geometry_reference.py
 * ========================================================================== */
#include "rt_ipss.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Guard threshold for Eq. (7): kappa is suppressed when the plane-parallel
 * single scattering underflows (workorder Sec. 2 step 3). */
#define IPSS_I1PP_FLOOR 1.0e-300

/* ---------------------------------------------------------------------------
 * Setup
 * ------------------------------------------------------------------------- */
int rt_ipss_init_from_atm(rt_ipss_t *G, const rt_atm_t *atm, double Re_km)
{
    if (!G || !atm || !atm->h || !atm->z_km_level) return -1;
    const int nt = atm->n_layers;
    if (nt < 1) return -1;
    if (!(atm->mu_sun > 0.0 && atm->mu_sun <= 1.0)) return -1;

    memset(G, 0, sizeof(*G));
    G->nlay   = nt;
    G->Re     = (Re_km > 0.0) ? Re_km : RT_IPSS_EARTH_RADIUS_KM;
    G->mu_sun = atm->mu_sun;
    G->n_quad = RT_IPSS_N_QUAD_DEFAULT;
    {   /* convergence knob for the gates; production uses the default */
        const char *nq = getenv("OCRT_IPSS_NQUAD");
        if (nq && nq[0]) {
            const int v = atoi(nq);
            if (v >= 8 && v <= 4096) G->n_quad = v;
        }
    }

    G->h_edge = (double *)malloc((size_t)(nt + 1) * sizeof(double));
    G->beta   = (double *)malloc((size_t)nt * sizeof(double));
    G->omega  = (double *)malloc((size_t)nt * sizeof(double));
    if (!G->h_edge || !G->beta || !G->omega) { rt_ipss_free(G); return -1; }

    /* OCRT stores k = 0 at TOA and k = nt at the ground; reverse to ascending
     * shells so that h_edge[0] = 0 km. */
    for (int k = 0; k <= nt; ++k)
        G->h_edge[k] = atm->z_km_level[nt - k];
    G->H = G->h_edge[nt];

    for (int j = 0; j < nt; ++j) {
        /* Ascending layer j spans [h_edge[j], h_edge[j+1]] and corresponds to
         * the OCRT layer between levels (nt-1-j) and (nt-j). */
        const int k    = nt - 1 - j;
        const double dz   = G->h_edge[j + 1] - G->h_edge[j];
        const double dtau = atm->h[k + 1] - atm->h[k];
        if (!(dz > 0.0)) { rt_ipss_free(G); return -2; }
        G->beta[j] = (dtau > 0.0) ? (dtau / dz) : 0.0;

        /* Scattering fraction of this layer.  atm->tau_abs_layer is indexed in
         * OCRT layer order and may be absent (pure scattering). */
        double tau_abs = 0.0;
        if (atm->tau_abs_layer) tau_abs = atm->tau_abs_layer[k];
        if (tau_abs < 0.0) tau_abs = 0.0;
        if (tau_abs > dtau) tau_abs = dtau;
        G->omega[j] = (dtau > 0.0) ? ((dtau - tau_abs) / dtau) : 0.0;
    }

    /* v1.11.1: phase-weighted source.  atm->xdel/ydel[k+1] are the aerosol /
     * Rayleigh SCATTERING fractions of the extinction of OCRT layer k (set by
     * rt_atm_build_aerosol_rayleigh and renormalized by the gas step; verified
     * omega_k == ydel + xdel/ssa_a to 1e-4 on the 400-layer production grid). */
    G->L_aer = -1; G->n_ang_aer = 0; G->depol = atm->depol;
    if (atm->xdel && atm->ydel) {
        int have_aer = 0;
        for (int k = 0; k < nt; ++k) if (atm->xdel[k + 1] > 0.0) { have_aer = 1; break; }
        const int have_tab = (atm->aer_use_value_kernel && atm->aer_n_ang_phase > 1 &&
                              atm->aer_theta_phase && atm->aer_P11_phase);
        const int have_mom = (atm->aerosol_active && atm->betal_aer && atm->L_max >= 0);
        /* Rayleigh-only profiles keep the phase-free path (bit-identical to
         * v1.11: P_R cancels exactly); the weights matter only with aerosol. */
        if (have_aer && (have_tab || have_mom)) {
            G->w_ray = (double *)malloc((size_t)nt * sizeof(double));
            G->w_aer = (double *)malloc((size_t)nt * sizeof(double));
            if (!G->w_ray || !G->w_aer) { rt_ipss_free(G); return -1; }
            for (int j = 0; j < nt; ++j) {
                const int k = nt - 1 - j;
                G->w_ray[j] = (atm->ydel[k + 1] > 0.0) ? atm->ydel[k + 1] : 0.0;
                G->w_aer[j] = (atm->xdel[k + 1] > 0.0) ? atm->xdel[k + 1] : 0.0;
            }
            if (have_tab) {
                G->n_ang_aer = atm->aer_n_ang_phase;
                G->theta_aer = atm->aer_theta_phase;
                G->p11_aer   = atm->aer_P11_phase;
            }
            if (have_mom) {
                G->L_aer = atm->L_max;
                G->betal_aer = (double *)malloc((size_t)(atm->L_max + 1) * sizeof(double));
                if (!G->betal_aer) { rt_ipss_free(G); return -1; }
                memcpy(G->betal_aer, atm->betal_aer, (size_t)(atm->L_max + 1) * sizeof(double));
            }
        } else if (have_aer) {
            static int warned = 0;
            if (!warned) {
                warned = 1;
                fprintf(stderr, "# ipss: aerosol present but no P11 data on rt_atm_t; "
                                "kappa falls back to phase-free weights\n");
            }
        }
    }
    return rt_ipss_finalize(G);
}

/* ---------------------------------------------------------------------------
 * Phase-weighted source (see rt_ipss.h, PHASE-WEIGHTED SOURCE)
 * ------------------------------------------------------------------------- */
double rt_ipss_phase_rayleigh(double depol, double cos_theta)
{
    /* P_R = 1 + b2 P2(mu), b2 = (1 - d)/(2 + d): the depolarized Rayleigh
     * phase function normalized to (1/2) int P dmu = 1 (same moment OCRT's
     * rt_molecular_phase_coefficients uses for the SOS). */
    const double d  = (depol >= 0.0 && depol < 1.0) ? depol : 0.0;
    const double b2 = (1.0 - d) / (2.0 + d);
    return 1.0 + b2 * 0.5 * (3.0 * cos_theta * cos_theta - 1.0);
}

double rt_ipss_phase_aerosol(const rt_ipss_t *G, double cos_theta)
{
    if (!G) return 1.0;
    if (G->n_ang_aer > 1 && G->theta_aer && G->p11_aer) {
        /* linear interpolation in angle on the ascending table */
        double th = acos(fmax(-1.0, fmin(1.0, cos_theta))) * 180.0 / M_PI;
        const int n = G->n_ang_aer;
        if (th <= G->theta_aer[0]) return G->p11_aer[0];
        if (th >= G->theta_aer[n - 1]) return G->p11_aer[n - 1];
        int lo = 0, hi = n - 1;
        while (hi - lo > 1) {
            const int mid = (lo + hi) >> 1;
            if (G->theta_aer[mid] <= th) lo = mid; else hi = mid;
        }
        const double t = (th - G->theta_aer[lo]) / (G->theta_aer[hi] - G->theta_aer[lo]);
        return G->p11_aer[lo] + t * (G->p11_aer[hi] - G->p11_aer[lo]);
    }
    if (G->L_aer >= 0 && G->betal_aer) {
        double p0 = 1.0, p1 = cos_theta, s = G->betal_aer[0];
        if (G->L_aer >= 1) s += G->betal_aer[1] * p1;
        for (int l = 2; l <= G->L_aer; ++l) {
            const double p2 = ((2.0 * l - 1.0) * cos_theta * p1 - (l - 1.0) * p0) / l;
            s += G->betal_aer[l] * p2; p0 = p1; p1 = p2;
        }
        return (s > 0.0) ? s : 0.0;
    }
    return 1.0;    /* isotropic fallback (never reached when init succeeded) */
}

double rt_ipss_cos_scattering_angle(double mu_sun, double vza_deg, double raa_deg)
{
    const double sn = sqrt(fmax(0.0, 1.0 - mu_sun * mu_sun));
    const double tv = vza_deg * M_PI / 180.0;
    const double fp = (180.0 - raa_deg) * M_PI / 180.0;
    /* V0 . u with V0 = (sn, 0, -mu_sun), u = (sin tv cos fp, sin tv sin fp, cos tv) */
    return sn * sin(tv) * cos(fp) - mu_sun * cos(tv);
}

void rt_ipss_set_source_weights(rt_ipss_t *G, double vza_deg, double raa_deg)
{
    if (!G || !G->w_ray || !G->w_aer || !G->w_src) return;
    const double ct = rt_ipss_cos_scattering_angle(G->mu_sun, vza_deg, raa_deg);
    const double PR = rt_ipss_phase_rayleigh(G->depol, ct);
    const double PA = rt_ipss_phase_aerosol(G, ct);
    G->cos_theta = ct; G->P_R = PR; G->P_A = PA;
    for (int j = 0; j < G->nlay; ++j)
        G->w_src[j] = G->w_ray[j] * PR + G->w_aer[j] * PA;
}

/* layer source weight used by both single-scatter integrals */
static inline double ipss_src(const rt_ipss_t *G, int j)
{
    return (G->w_src) ? G->w_src[j] : G->omega[j];
}

int rt_ipss_finalize(rt_ipss_t *G)
{
    if (!G || G->nlay < 1 || !G->h_edge || !G->beta || !G->omega) return -1;
    const int nt = G->nlay;
    free(G->scr_seg); free(G->scr_ts); free(G->tau_above);
    G->scr_seg   = (double *)malloc((size_t)nt * sizeof(double));
    G->scr_ts    = (double *)malloc((size_t)(2 * nt + 8) * sizeof(double));
    G->tau_above = (double *)malloc((size_t)(nt + 1) * sizeof(double));
    free(G->los_t); free(G->los_tau); free(G->los_lay);
    G->los_t   = (double *)malloc((size_t)(2 * nt + 8) * sizeof(double));
    G->los_tau = (double *)malloc((size_t)(2 * nt + 8) * sizeof(double));
    G->los_lay = (int *)malloc((size_t)(2 * nt + 8) * sizeof(int));
    G->los_n = 0;
    if (!G->scr_seg || !G->scr_ts || !G->tau_above ||
        !G->los_t || !G->los_tau || !G->los_lay) { rt_ipss_free(G); return -1; }
    free(G->w_src); G->w_src = NULL;
    if (G->w_ray && G->w_aer) {          /* phase-weighted profile */
        G->w_src = (double *)malloc((size_t)nt * sizeof(double));
        if (!G->w_src) { rt_ipss_free(G); return -1; }
        for (int j = 0; j < nt; ++j) G->w_src[j] = G->w_ray[j] + G->w_aer[j];
        G->cos_theta = 2.0;              /* "not set" marker */
    }
    /* vertical optical depth from each ascending edge up to the TOA */
    G->tau_above[nt] = 0.0;
    for (int j = nt - 1; j >= 0; --j)
        G->tau_above[j] = G->tau_above[j + 1]
                          + G->beta[j] * (G->h_edge[j + 1] - G->h_edge[j]);
    return 0;
}

void rt_ipss_free(rt_ipss_t *G)
{
    if (!G) return;
    free(G->h_edge); free(G->beta); free(G->omega);
    free(G->w_ray); free(G->w_aer); free(G->w_src); free(G->betal_aer);
    free(G->scr_seg); free(G->scr_ts); free(G->tau_above);
    free(G->gl_x); free(G->gl_w);
    free(G->los_t); free(G->los_tau); free(G->los_lay);
    memset(G, 0, sizeof(*G));
}

/* ---------------------------------------------------------------------------
 * Parametric ray-shell intersection (production path; no sine-rule singularity,
 * automatic non-uniform layering, automatic shadow test).
 * seg[j] receives the path length [km] inside ascending layer j.
 * hit_ground is set when the ray strikes the surface before leaving the TOA.
 * ------------------------------------------------------------------------- */
int rt_ipss_ray_segments(const rt_ipss_t *G, const double r0[3],
                         const double u[3], double *seg, int *hit_ground)
{
    if (!G || !r0 || !u || !seg) return -1;
    const int nlay = G->nlay;
    memset(seg, 0, (size_t)nlay * sizeof(double));
    if (hit_ground) *hit_ground = 0;

    const double un = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    if (!(un > 0.0)) return -1;
    const double ux = u[0] / un, uy = u[1] / un, uz = u[2] / un;
    const double b  = r0[0] * ux + r0[1] * uy + r0[2] * uz;
    const double c0 = r0[0] * r0[0] + r0[1] * r0[1] + r0[2] * r0[2];

    /* First ground intersection truncates the ray (shadow test). */
    double t_stop = INFINITY;
    const double disc_g = b * b - (c0 - G->Re * G->Re);
    if (disc_g >= 0.0) {
        const double sq = sqrt(disc_g);
        const double cand[2] = { -b - sq, -b + sq };
        for (int i = 0; i < 2; ++i)
            if (cand[i] > 1e-9 && cand[i] < t_stop) {
                t_stop = cand[i];
                if (hit_ground) *hit_ground = 1;
            }
    }

    /* Collect shell crossings; r(t)^2 = t^2 + 2 b t + c0.
     * The list is built ALREADY SORTED, which removes the O(n^2) insertion
     * sort from the hot path: for the descending root t- = -b - sqrt(.) the
     * value decreases with shell radius, and for the ascending root
     * t+ = -b + sqrt(.) it increases, with t- <= -b <= t+ always.  Walking the
     * shells downward for t- and upward for t+ therefore emits monotonically
     * increasing parameters. */
    const int cap = 2 * nlay + 8;
    double *ts = G->scr_ts;      /* preallocated; hot path allocates nothing */
    if (!ts) return -1;
    int n = 0;
    ts[n++] = 0.0;
    for (int k = nlay; k >= 0; --k) {
        const double R = G->Re + G->h_edge[k];
        const double disc = b * b - (c0 - R * R);
        if (disc < 0.0) continue;
        const double t = -b - sqrt(disc);
        if (t > 1e-9 && t <= t_stop && n < cap) ts[n++] = t;
    }
    for (int k = 0; k <= nlay; ++k) {
        const double R = G->Re + G->h_edge[k];
        const double disc = b * b - (c0 - R * R);
        if (disc < 0.0) continue;
        const double t = -b + sqrt(disc);
        if (t > 1e-9 && t <= t_stop && n < cap) ts[n++] = t;
    }
    if (hit_ground && *hit_ground && n < cap) ts[n++] = t_stop;

    /* dedupe (already ordered; ts[0] = 0 is the smallest by construction) */
    int m = 0;
    for (int i = 0; i < n; ++i)
        if (i == 0 || ts[i] > ts[m - 1] + 1e-12) ts[m++] = ts[i];

    /* Assign each interval to a layer by its midpoint radius. */
    for (int i = 0; i + 1 < m; ++i) {
        const double ta = ts[i], tb = ts[i + 1];
        const double tm = 0.5 * (ta + tb);
        const double rm = sqrt(tm * tm + 2.0 * b * tm + c0);
        const double hm = rm - G->Re;
        if (hm < 0.0 || hm > G->h_edge[nlay]) continue;
        int lo = 0, hi = nlay;              /* find j with h_edge[j] <= hm */
        while (hi - lo > 1) {
            const int mid = (lo + hi) / 2;
            if (G->h_edge[mid] <= hm) lo = mid; else hi = mid;
        }
        if (lo > nlay - 1) lo = nlay - 1;
        seg[lo] += (tb - ta);
    }
    return 0;
}

static double ipss_tau_of_segments(const rt_ipss_t *G, const double *seg)
{
    double tau = 0.0;
    for (int j = 0; j < G->nlay; ++j) tau += seg[j] * G->beta[j];
    return tau;
}

/* Slant optical depth of the solar beam from point rB up to the TOA.
 * Returns INFINITY when the point lies in the Earth's shadow. */
double rt_ipss_tau_sun_at(const rt_ipss_t *G, const double rB[3])
{
    if (!G || !rB) return INFINITY;
    const double mu = G->mu_sun;
    const double sn = sqrt(fmax(0.0, 1.0 - mu * mu));
    /* Solar beam propagates downward as V0 = (sin ts, 0, -cos ts); tracing the
     * beam back toward the sun means marching along -V0. */
    const double u[3] = { -sn, 0.0, mu };
    double *seg = G->scr_seg;
    if (!seg) return INFINITY;
    int hit = 0;
    if (rt_ipss_ray_segments(G, rB, u, seg, &hit) != 0) return INFINITY;
    return hit ? INFINITY : ipss_tau_of_segments(G, seg);
}

/* Slant view optical depth from TOA down to the ground along the line of sight
 * (Eq. 5 with the exact (Re + h_i) radius, obtained here by ray tracing). */
double rt_ipss_tau_view_toa(const rt_ipss_t *G, double vza_deg)
{
    if (!G) return 0.0;
    const double tv = vza_deg * M_PI / 180.0;
    const double r0[3] = { 0.0, 0.0, G->Re + G->H };
    /* Marching from the sensor DOWN into the atmosphere is -V. */
    const double u[3] = { -sin(tv), 0.0, -cos(tv) };
    double *seg = G->scr_seg;
    if (!seg) return 0.0;
    int hit = 0;
    rt_ipss_ray_segments(G, r0, u, seg, &hit);
    return ipss_tau_of_segments(G, seg);
}

double rt_ipss_tau_dir_dn_boa(const rt_ipss_t *G)
{
    if (!G) return -1.0;
    const double rB[3] = { 0.0, 0.0, G->Re };
    const double tau = rt_ipss_tau_sun_at(G, rB);
    return isfinite(tau) ? tau : -1.0;
}

/* ---------------------------------------------------------------------------
 * Gauss-Legendre nodes (Newton on Legendre polynomials).
 * ------------------------------------------------------------------------- */
/* Cached Gauss-Legendre nodes for the object's quadrature size. */
static int ipss_gl_cache(rt_ipss_t *G, int n);

static void ipss_gauss_legendre(int n, double *x, double *w)
{
    for (int i = 0; i < (n + 1) / 2; ++i) {
        double z = cos(M_PI * (i + 0.75) / (n + 0.5));
        double pp = 0.0;
        for (int it = 0; it < 100; ++it) {
            double p0 = 1.0, p1 = 0.0;
            for (int j = 0; j < n; ++j) {
                const double p2 = p1;
                p1 = p0;
                p0 = ((2.0 * j + 1.0) * z * p1 - j * p2) / (j + 1.0);
            }
            pp = n * (z * p0 - p1) / (z * z - 1.0);
            const double dz = -p0 / pp;
            z += dz;
            if (fabs(dz) < 1e-15) break;
        }
        x[i]         = -z;
        x[n - 1 - i] =  z;
        w[i]         = 2.0 / ((1.0 - z * z) * pp * pp);
        w[n - 1 - i] = w[i];
    }
}

static int ipss_gl_cache(rt_ipss_t *G, int n)
{
    if (G->gl_x && G->gl_n == n) return 0;
    free(G->gl_x); free(G->gl_w);
    G->gl_x = (double *)malloc((size_t)n * sizeof(double));
    G->gl_w = (double *)malloc((size_t)n * sizeof(double));
    if (!G->gl_x || !G->gl_w) { free(G->gl_x); free(G->gl_w);
                                G->gl_x = G->gl_w = NULL; G->gl_n = 0; return -1; }
    ipss_gauss_legendre(n, G->gl_x, G->gl_w);
    G->gl_n = n;
    return 0;
}

/* Altitude -> layer index (ascending grid). */
static int ipss_layer_of(const rt_ipss_t *G, double h)
{
    int lo = 0, hi = G->nlay;
    if (h <= G->h_edge[0]) return 0;
    if (h >= G->h_edge[G->nlay]) return G->nlay - 1;
    while (hi - lo > 1) {
        const int mid = (lo + hi) / 2;
        if (G->h_edge[mid] <= h) lo = mid; else hi = mid;
    }
    return (lo > G->nlay - 1) ? G->nlay - 1 : lo;
}

/* Geometric length of the line of sight from TOA to the ground. */
static double ipss_view_path_length(const rt_ipss_t *G, double vza_deg)
{
    const double tv = vza_deg * M_PI / 180.0;
    const double rA = G->Re + G->H;
    const double s  = sin(tv);
    /* Solve |r_A - d V| = Re for the near root; V is the upward view unit
     * vector so the downward march is -V. */
    const double disc = G->Re * G->Re - rA * rA * s * s;
    if (disc <= 0.0) return 0.0;                 /* limb path: no ground hit */
    return rA * cos(tv) - sqrt(disc);
}

/* Build the ascending crossing table of the line of sight r(s) = rA - s V and
 * the cumulative optical depth at each crossing.  One pass per geometry turns
 * every later tau_view(s) query into a binary search. */
static void ipss_build_los(rt_ipss_t *G, const double rA[3], const double V[3],
                           double d_max)
{
    const int nlay = G->nlay;
    const double b  = -(rA[0] * V[0] + rA[1] * V[1] + rA[2] * V[2]);
    const double c0 = rA[0] * rA[0] + rA[1] * rA[1] + rA[2] * rA[2];
    const int cap = 2 * nlay + 8;
    int n = 0;
    G->los_t[n++] = 0.0;
    for (int k = nlay; k >= 0; --k) {
        const double R = G->Re + G->h_edge[k];
        const double disc = b * b - (c0 - R * R);
        if (disc < 0.0) continue;
        const double t = -b - sqrt(disc);
        if (t > 1e-12 && t <= d_max && n < cap) G->los_t[n++] = t;
    }
    for (int k = 0; k <= nlay; ++k) {
        const double R = G->Re + G->h_edge[k];
        const double disc = b * b - (c0 - R * R);
        if (disc < 0.0) continue;
        const double t = -b + sqrt(disc);
        if (t > 1e-12 && t <= d_max && n < cap) G->los_t[n++] = t;
    }
    if (n < cap) G->los_t[n++] = d_max;
    int m = 0;
    for (int i = 0; i < n; ++i)
        if (i == 0 || G->los_t[i] > G->los_t[m - 1] + 1e-12) G->los_t[m++] = G->los_t[i];
    double tau = 0.0;
    for (int i = 0; i + 1 < m; ++i) {
        const double tm = 0.5 * (G->los_t[i] + G->los_t[i + 1]);
        const double rm = sqrt(tm * tm + 2.0 * b * tm + c0);
        const int j = ipss_layer_of(G, rm - G->Re);
        G->los_lay[i] = j;
        G->los_tau[i] = tau;
        tau += (G->los_t[i + 1] - G->los_t[i]) * G->beta[j];
    }
    G->los_tau[m - 1] = tau;
    G->los_lay[m - 1] = G->los_lay[(m > 1) ? m - 2 : 0];
    G->los_n = m;
}

static double ipss_los_tau_at(const rt_ipss_t *G, double s)
{
    int lo = 0, hi = G->los_n - 1;
    if (s <= G->los_t[0]) return 0.0;
    if (s >= G->los_t[hi]) return G->los_tau[hi];
    while (hi - lo > 1) {
        const int mid = (lo + hi) / 2;
        if (G->los_t[mid] <= s) lo = mid; else hi = mid;
    }
    return G->los_tau[lo] + (s - G->los_t[lo]) * G->beta[G->los_lay[lo]];
}

/* Set up the parametric line of sight r(s) = rA - s V, s in [0, d_max], for
 * one of the two anchor conventions (see rt_ipss.h "VZA ANCHOR").
 *   TOA anchor     : vza_deg is the zenith angle AT TOA; the sensor sits on the
 *                    nadir column and the ray reaches the ground off-axis.  The
 *                    ray does not reach the ground beyond the TOA horizon
 *                    (79.86 deg for H = 100 km), where d_max <= 0.
 *   SURFACE anchor : vza_deg is the zenith angle AT THE TARGET PIXEL (0,0,Re) —
 *                    the satellite L1B senz convention — and the TOA crossing
 *                    is displaced laterally.  Always geometrically valid for
 *                    vza < 90 deg.
 * Returns 0 on success, -1 if no ground-connected path exists. */
static int ipss_los_setup(const rt_ipss_t *G, double vza_deg, double raa_deg,
                          rt_ipss_vza_anchor_t anchor,
                          double rA[3], double V[3], double *d_max)
{
    const double tv = vza_deg * M_PI / 180.0;
    /* OCRT public RAA -> paper azimuth (glint half-plane at phi_paper = 0).
     * Verified against the OCRT direct-glint peak and the first-order Rayleigh
     * phase dependence: phi_paper = 180 deg - raa. */
    const double fp = (180.0 - raa_deg) * M_PI / 180.0;
    const double u[3] = { sin(tv) * cos(fp), sin(tv) * sin(fp), cos(tv) };
    V[0] = u[0]; V[1] = u[1]; V[2] = u[2];
    if (anchor == RT_IPSS_VZA_ANCHOR_TOA) {
        rA[0] = 0.0; rA[1] = 0.0; rA[2] = G->Re + G->H;
        *d_max = ipss_view_path_length(G, vza_deg);
    } else {
        /* Pixel P = (0,0,Re); solve |P + s u| = Re + H for the positive root. */
        const double R = G->Re + G->H;
        const double bq = G->Re * u[2];               /* P.u */
        const double disc = bq * bq + R * R - G->Re * G->Re;
        if (disc <= 0.0) { *d_max = 0.0; return -1; }
        const double s = -bq + sqrt(disc);
        if (!(s > 0.0)) { *d_max = 0.0; return -1; }
        rA[0] = s * u[0]; rA[1] = s * u[1]; rA[2] = G->Re + s * u[2];
        *d_max = s;
    }
    return (*d_max > 0.0) ? 0 : -1;
}

static double ipss_single_scatter_core(const rt_ipss_t *G, double vza_deg,
                                       double raa_deg, rt_ipss_sun_mode_t mode,
                                       rt_ipss_vza_anchor_t anchor)
{
    if (!G) return 0.0;
    /* SEGMENT-WISE quadrature.  beta is piecewise constant with a jump at every
     * shell boundary, so a single global Gauss rule over the whole line of
     * sight converges erratically on a realistic (US62) grid — measured ~1%
     * wobble between 32 and 512 nodes.  Integrating each layer crossing
     * separately makes the integrand smooth inside every panel; with OCRT's
     * equal-d(tau) layering the per-panel error is O(dtau^2) and a 2-point
     * Gauss rule is already converged.  n_quad is reinterpreted as nodes per
     * panel (OCRT_IPSS_NQUAD overrides it for the convergence gate). */
    const int nsub = (G->n_quad >= 1 && G->n_quad <= 64) ? G->n_quad : 2;
    double rA[3], V[3], d_max = 0.0;
    if (ipss_los_setup(G, vza_deg, raa_deg, anchor, rA, V, &d_max) != 0) return 0.0;
    if (!(d_max > 0.0)) return 0.0;

    rt_ipss_t *Gm = (rt_ipss_t *)G;          /* memo caches only; no state change */
    ipss_build_los(Gm, rA, V, d_max);
    if (ipss_gl_cache(Gm, nsub) != 0) return 0.0;
    const double *x = G->gl_x, *w = G->gl_w;

    const double mu_s = G->mu_sun;
    double I1 = 0.0;
    for (int i = 0; i + 1 < G->los_n; ++i) {
        const int j = G->los_lay[i];
        if (!(G->beta[j] > 0.0)) continue;
        const double sa = G->los_t[i], sb = G->los_t[i + 1];
        const double half = 0.5 * (sb - sa), mid = 0.5 * (sa + sb);
        for (int q = 0; q < nsub; ++q) {
            const double s   = mid + half * x[q];
            const double wgt = half * w[q];
            const double rB[3] = { rA[0] - s * V[0], rA[1] - s * V[1],
                                   rA[2] - s * V[2] };
            const double rn = sqrt(rB[0] * rB[0] + rB[1] * rB[1] + rB[2] * rB[2]);
            const double hB = rn - G->Re;
            if (hB < 0.0 || hB > G->H) continue;
            const double tau_view = G->los_tau[i] + (s - sa) * G->beta[j];

            double tau_sun;
            if (mode == RT_IPSS_SUN_SPHERICAL) {
                tau_sun = rt_ipss_tau_sun_at(G, rB);
            } else if (mode == RT_IPSS_SUN_NADIR) {
                const double rN[3] = { 0.0, 0.0, G->Re + hB };
                tau_sun = rt_ipss_tau_sun_at(G, rN);
            } else {
                const int jj = ipss_layer_of(G, hB);
                tau_sun = (G->tau_above[jj + 1]
                           + (G->h_edge[jj + 1] - hB) * G->beta[jj]) / mu_s;
            }
            if (!isfinite(tau_sun)) continue;     /* shadowed: no source */
            I1 += wgt * G->beta[j] * ipss_src(G, j) * exp(-tau_view - tau_sun);
        }
    }
    return I1;
}

double rt_ipss_single_scatter(const rt_ipss_t *G, double vza_deg, double raa_deg,
                              rt_ipss_sun_mode_t mode)
{
    return ipss_single_scatter_core(G, vza_deg, raa_deg, mode,
                                    RT_IPSS_VZA_ANCHOR_TOA);
}

double rt_ipss_single_scatter_surf(const rt_ipss_t *G, double vza_deg,
                                   double raa_deg, rt_ipss_sun_mode_t mode)
{
    return ipss_single_scatter_core(G, vza_deg, raa_deg, mode,
                                    RT_IPSS_VZA_ANCHOR_SURFACE);
}

double rt_ipss_single_scatter_pp(const rt_ipss_t *G, double vza_deg)
{
    if (!G) return 0.0;
    const double mu_v = cos(vza_deg * M_PI / 180.0);
    if (!(mu_v > 0.0)) return 0.0;
    const double mu_s = G->mu_sun;
    const double C = 1.0 / mu_v + 1.0 / mu_s;
    /* Plane-parallel single scattering is analytic per layer: with beta and
     * omega constant inside layer j and tau measured from the TOA,
     *   dI = (omega/mu_v) exp(-C tau) dtau  =>  layer contribution
     *      = (omega/(mu_v C)) [exp(-C tau_top) - exp(-C tau_bot)].
     * Exactness here keeps kappa free of any denominator quadrature error. */
    double I1 = 0.0;
    for (int j = 0; j < G->nlay; ++j) {
        if (!(G->beta[j] > 0.0)) continue;
        const double tau_top = G->tau_above[j + 1];
        const double tau_bot = G->tau_above[j];
        I1 += ipss_src(G, j) / (mu_v * C) * (exp(-C * tau_top) - exp(-C * tau_bot));
    }
    return I1;
}

int rt_ipss_kappa(const rt_ipss_t *G, double vza_deg, double raa_deg,
                  double *kappa, double *I1_ss_out, double *I1_pp_out, int *guard)
{
    if (!G || !kappa) return -1;
    /* SURFACE anchor: vza_deg is the zenith angle at the target pixel, exactly
     * what the plane-parallel solve used for mu_v and what satellite L1B
     * reports as senz.  Numerator and denominator then share one geometry, so
     * kappa -> 1 whenever the atmosphere is thin or collapses onto the surface
     * (gate G-I3).  Anchoring the spherical ray at TOA instead leaves a
     * spurious VZA-dependent offset (-2.9 % at VZA 55, -11 % at VZA 70) that
     * survives even with the sun at the zenith; see rt_ipss.h. */
    /* v1.11.1: real single-scattering source weights for this geometry
     * (no-op on phase-free profiles). */
    rt_ipss_set_source_weights((rt_ipss_t *)G, vza_deg, raa_deg);
    const double I1_ss = rt_ipss_single_scatter_surf(G, vza_deg, raa_deg,
                                                     RT_IPSS_SUN_SPHERICAL);
    const double I1_pp = rt_ipss_single_scatter_pp(G, vza_deg);
    if (I1_ss_out) *I1_ss_out = I1_ss;
    if (I1_pp_out) *I1_pp_out = I1_pp;
    if (guard) *guard = 0;
    if (!(I1_pp > IPSS_I1PP_FLOOR) || !isfinite(I1_ss)) {
        *kappa = 1.0;                 /* no correction rather than 0/0 */
        if (guard) *guard = 1;
        return 0;
    }
    *kappa = I1_ss / I1_pp;
    return 0;
}
