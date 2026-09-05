#include <stdio.h>
#include <stdlib.h>
/* =============================================================================
 * rt_air_water_coupling.c — atm BOA <-> in-water boundary field mapping
 * =============================================================================
 * Phase B.4 Stage 2b step2 — see header for the API description.
 *
 * ROLE
 *   Bidirectional radiometric coupling across the air-water interface for the
 *   coupled atmosphere + ocean SOS solve:
 *
 *   FORWARD  (rt_air_water_couple_atm_to_water):
 *     downwelling atm field at the bottom of atmosphere (BOA, z = 0+)
 *       -> refracted downwelling field just below the surface (z = 0-)
 *          on the in-water Gauss mu grid, per azimuth Fourier mode m.
 *   REVERSE  (rt_air_water_couple_water_to_atm):
 *     upwelling in-water field at z = 0-  ->  transmitted upwelling field at
 *     z = 0+ on the atm mu grid (used as the atm bottom boundary source).
 *
 * FIELD REPRESENTATION / INDEXING
 *   All Stokes fields are azimuth-Fourier decomposed:
 *     X(mu, phi) = X^0(mu) + 2 * sum_{m=1..m_max} X^m(mu) * cos(m phi)
 *     (U uses sin(m phi));  storage is flat: X_per_m[m * n_mu + j].
 *   mu grids are POSITIVE, sorted ascending; direction sign is implied by
 *   context (BOA export = downwelling, water 0- export = upwelling).
 *
 * PHYSICS SUMMARY (equations live at each code block below)
 *   1. Flat interface, forward: per water node mu_w, Snell-map to mu_air,
 *      linear-interp the BOA Stokes at mu_air, multiply by the flat Fresnel
 *      transmission Mueller T_aw(mu_air).  T_aw INCLUDES the n^2 radiance
 *      (basic-radiance / etendue) factor: L_water = n_w^2 * T_fresnel * L_air
 *      along the refracted ray, so no separate n^2 is applied here.
 *      Water nodes inside the total-internal-reflection (TIR) cone
 *      (mu_w < mu_crit = sqrt(1 - 1/n_w^2), ~0.6655 for n_w = 1.34) receive
 *      ZERO direct atmospheric transmission by construction.
 *   2. Rough surface, forward (FIX-B, wind > 0): the m = 0 intensity column
 *      is REPLACED by the Cox-Munk microfacet BTDF integral (block comment at
 *      the FIX-B section); m >= 1 and Q/U keep flat-Fresnel values, which is
 *      harmless because only m = 0 I is consumed downstream (Ed / Lu
 *      equivalent-beam superposition).  KNOWN OPEN DEFECT (-1.04%, "skylight
 *      Ed/Lu"): the Ed normalization integral and the Lu equivalent-beam field
 *      treat the sub-critical-cone leakage inconsistently; see
 *      CHANGELOG_FIX-SKY-EDLU doc.  Do not re-derive here.
 *   3. Reverse: flat path Snell-maps each air node to its unique specular
 *      water direction and applies T_wa (which includes the UPWARD 1/n^2
 *      radiance factor); rough path (Milestone 2b, wind > 0) replaces the
 *      specular map by the full Fourier-mode Cox-Munk transmission kernel
 *      T^m[a,w] (block comment at that section).
 *
 * INTERPOLATION RULE
 *   All grid resampling is LINEAR (never nearest-neighbor, project rule);
 *   out-of-range targets use 2-point linear extrapolation from the nearest
 *   edge pair.
 * ========================================================================== */
#include "rt_air_water_coupling.h"
#include "rt_air_water.h"
#include "rt_quadrature.h"
#include "shared/surface.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int rt_aw_coupled_field_alloc(rt_aw_coupled_field_t *f, int m_max, int n_mu_water) {
    if (!f || m_max < 0 || n_mu_water < 1) return -1;
    memset(f, 0, sizeof(*f));
    f->m_max      = m_max;
    f->n_mu_water = n_mu_water;
    size_t sz = (size_t)(m_max + 1) * (size_t)n_mu_water;
    f->mu_water_pos     = (double*)calloc((size_t)n_mu_water, sizeof(double));
    f->I_inwater_per_m  = (double*)calloc(sz, sizeof(double));
    f->Q_inwater_per_m  = (double*)calloc(sz, sizeof(double));
    f->U_inwater_per_m  = (double*)calloc(sz, sizeof(double));
    if (!f->mu_water_pos || !f->I_inwater_per_m ||
        !f->Q_inwater_per_m || !f->U_inwater_per_m) {
        rt_aw_coupled_field_free(f);
        return -1;
    }
    return 0;
}

void rt_aw_coupled_field_free(rt_aw_coupled_field_t *f) {
    if (!f) return;
    if (f->mu_water_pos)    free(f->mu_water_pos);
    if (f->I_inwater_per_m) free(f->I_inwater_per_m);
    if (f->Q_inwater_per_m) free(f->Q_inwater_per_m);
    if (f->U_inwater_per_m) free(f->U_inwater_per_m);
    memset(f, 0, sizeof(*f));
}

/* Linear interpolation on a sorted-ascending μ grid.
 *   Returns linearly-interpolated value at μ_target.
 *   Out-of-bounds uses linear *extrapolation* from the nearest two-point pair
 *   (per user rule #11: no nearest-neighbor).
 *
 * Assumes mu_grid[0..n-1] is sorted ascending. n >= 2.
 */
static double linear_interp_ascending(const double *mu_grid, const double *vals,
                                       int n, double mu_target) {
    if (n < 2) return (n == 1) ? vals[0] : 0.0;
    /* find smallest i with mu_grid[i] >= mu_target */
    int hi = 1;
    while (hi < n && mu_grid[hi] < mu_target) hi++;
    int lo = hi - 1;
    if (hi >= n) { hi = n - 1; lo = n - 2; }  /* extrapolate high */
    if (lo < 0)  { lo = 0; hi = 1; }           /* extrapolate low */
    double t = (mu_target - mu_grid[lo]) / (mu_grid[hi] - mu_grid[lo]);
    return (1.0 - t) * vals[lo] + t * vals[hi];
}

/* Build a sorted, unique local view before interpolation.  Solver grids may
 * contain zero-weight view/solar/auxiliary nodes appended after the positive
 * quadrature nodes, so raw storage order is not guaranteed to be monotone.
 * Sorting at the consumer boundary makes interpolation independent of that
 * storage detail; already sorted inputs retain the identical arithmetic path. */
/* Coupling interpolation must retain every solver direction.  The water
 * solver appends zero-weight view/solar/nadir directions after the positive
 * Gauss nodes; silently truncating the local sort buffer therefore discards
 * the very directions required for exact output interpolation and changes the
 * operation into an extrapolation.  LUT production requests 96 water nodes,
 * so the legacy 64-entry cap affected a normal production path.
 *
 * The validated solver limit is 256 requested water nodes plus a small number
 * of appended directions.  A 320-entry stack workspace covers that path.  A
 * heap fallback is retained as a fail-loud guard for future larger grids; an
 * allocation failure must never be converted into a silent numerical clamp.
 * For n <= 64 the arithmetic order is unchanged from the legacy path. */
enum { AW_INTERP_STACK_N = 320 };

/* v1.11-speed S2 (2026-08-25): the coupling loops call this once per
 * (component, mode, target) with the SAME direction grid, and the per-call
 * insertion sort + dedup of that grid dominated the calls (measured 16.5% of
 * a fair-profile coupled run, 24.6M calls).  Memoize the sorted/deduped index
 * map per exact grid content (thread-local, byte-compared — no false hits);
 * values are gathered through the memoized map, so every arithmetic input to
 * linear_interp_ascending is bit-identical to the direct path. */
typedef struct {
    int    valid, n, w;
    double mu_copy[AW_INTERP_STACK_N];
    int    keep_src[AW_INTERP_STACK_N];   /* source index per kept entry */
    double ms[AW_INTERP_STACK_N];         /* sorted, deduped mu */
} aw_sortmemo_t;
static _Thread_local aw_sortmemo_t g_aw_sortmemo;

static double aw_interp_on_unsorted(const double *mu, const double *vals,
                                    int n, double mu_target) {
    int ord_s[AW_INTERP_STACK_N];
    double ms_s[AW_INTERP_STACK_N], tmp_s[AW_INTERP_STACK_N];
    int *ord = ord_s; double *ms = ms_s, *tmp = tmp_s;
    int *heap_i = NULL; double *heap_d = NULL;
    int w = 0;
    if (n < 2) return linear_interp_ascending(mu, vals, n, mu_target);
    if (n <= AW_INTERP_STACK_N) {
        aw_sortmemo_t *M = &g_aw_sortmemo;
        if (!(M->valid && M->n == n &&
              memcmp(M->mu_copy, mu, (size_t)n * sizeof(double)) == 0)) {
            /* (re)build the sorted/deduped map for this grid content —
             * identical insertion sort + tolerance dedup as the direct path */
            for (int a = 0; a < n; ++a) ord[a] = a;
            for (int a = 1; a < n; ++a) {
                int k = ord[a]; int b = a - 1;
                while (b >= 0 && mu[ord[b]] > mu[k]) { ord[b+1] = ord[b]; --b; }
                ord[b+1] = k;
            }
            int wm = 0;
            for (int a = 0; a < n; ++a) {
                double m = mu[ord[a]];
                if (wm > 0 && fabs(m - M->ms[wm-1]) < 1e-12) continue;
                M->ms[wm] = m; M->keep_src[wm] = ord[a]; ++wm;
            }
            memcpy(M->mu_copy, mu, (size_t)n * sizeof(double));
            M->n = n; M->w = wm; M->valid = 1;
        }
        for (int k = 0; k < M->w; ++k) tmp[k] = vals[M->keep_src[k]];
        return linear_interp_ascending(M->ms, tmp, M->w, mu_target);
    }
    if (n > AW_INTERP_STACK_N) {
        heap_i = (int *)malloc((size_t)n * sizeof(int));
        heap_d = (double *)malloc((size_t)n * 2u * sizeof(double));
        if (!heap_i || !heap_d) {
            fprintf(stderr,
                    "aw_interp_on_unsorted: failed to allocate workspace for "
                    "%d directions; refusing to truncate the interpolation "
                    "grid silently.\n", n);
            free(heap_i); free(heap_d);
            abort();
        }
        ord = heap_i; ms = heap_d; tmp = heap_d + n;
    }
    for (int a = 0; a < n; ++a) ord[a] = a;
    for (int a = 1; a < n; ++a) {        /* insertion sort, ascending mu */
        int k = ord[a]; int b = a - 1;
        while (b >= 0 && mu[ord[b]] > mu[k]) { ord[b+1] = ord[b]; --b; }
        ord[b+1] = k;
    }
    for (int a = 0; a < n; ++a) {        /* tolerance-based deduplication */
        double m = mu[ord[a]];
        if (w > 0 && fabs(m - ms[w-1]) < 1e-12) continue;
        ms[w] = m; tmp[w] = vals[ord[a]]; ++w;
    }
    double out = linear_interp_ascending(ms, tmp, w, mu_target);
    free(heap_i); free(heap_d);
    return out;
}

#ifdef OCRT_TEST_AW_INTERP
double rt_test_aw_interp_on_unsorted(const double *mu, const double *vals,
                                      int n, double mu_target) {
    return aw_interp_on_unsorted(mu, vals, n, mu_target);
}
#endif

/* Forward coupling: atm BOA downwelling field -> in-water z=0- field.
 *
 * Per water quadrature node mu_w (all Fourier modes m = 0..m_max):
 *
 *   Snell (flat interface):  sin(theta_air) = n_w * sin(theta_w)
 *     => mu_air = sqrt(1 - n_w^2 * (1 - mu_w^2))         [see block below]
 *   L_air^m(mu_air)  : linear interpolation of the BOA Fourier field
 *   L_water^m(mu_w)  = M_T_aw(mu_air) . [I,Q,U]_air^m    (3x3 Mueller)
 *
 * where M_T_aw = flat Fresnel air->water transmission Mueller INCLUDING the
 * n_w^2 downward basic-radiance factor (rt_air_water_T_aw).  For wind > 0
 * the m=0 intensity column is subsequently overwritten by the Cox-Munk
 * BTDF integral (FIX-B block below).
 *
 * Inputs:
 *   boa_export   : atm BOA Fourier field {I,Q,U}_per_m on mu_quad_pos
 *                  (ascending positive air mu nodes), plus m_max, n_mu
 *   n_water      : water refractive index (> 1; e.g. 1.34)
 *   q_convention : Q sign convention selector passed through to the Mueller
 *                  builders (must match the solver's Stokes convention)
 *   mu_water_pos : ascending positive in-water mu nodes [n_mu_water]
 *   wind_speed   : m/s; > 0 activates the FIX-B rough m=0 override
 *   sigma_type   : Cox-Munk slope-variance model selector (pass-through)
 * Output:
 *   coupled->{I,Q,U}_inwater_per_m[m * n_mu_water + j_w], plus a copy of the
 *   water mu grid.  TIR-cone nodes and invalid mu get exact zeros.
 * Returns 0 on success, -1 on argument/shape errors. */
int rt_air_water_couple_atm_to_water(const rt_atm_boa_export_t *boa_export,
                                       double n_water,
                                       int q_convention,
                                       const double *mu_water_pos,
                                       int n_mu_water,
                                       const double *w_atm_pos,
                                       double wind_speed,
                                       int sigma_type,
                                       rt_aw_coupled_field_t *coupled) {
    if (!boa_export || !mu_water_pos || !coupled) return -1;
    if (n_water <= 1.0) return -1;
    if (boa_export->n_mu < 2) return -1;
    if (coupled->n_mu_water != n_mu_water ||
        coupled->m_max != boa_export->m_max) return -1;
    if (!coupled->I_inwater_per_m || !coupled->Q_inwater_per_m ||
        !coupled->U_inwater_per_m || !coupled->mu_water_pos) return -1;

    /* When atmospheric quadrature weights are available, apply the full
     * polarized Cox-Munk air-to-water Fourier transmission operator.  Nodes
     * with zero quadrature weight (view/solar/auxiliary directions) contribute
     * exactly zero.  The fallback below handles flat water and weightless
     * single-direction callers.  The Fourier contraction is
     *   S_w^m(i) = C_m sum_j mu_a_j w_j T_aw^m[i][j] . S_air^m(j),
     *   C_0 = 2*pi, C_{m>0} = pi, no (-1)^m (transmission). */
    if (w_atm_pos && wind_speed > 0.0) {
        const int n_a = boa_export->n_mu;
        const int mM  = boa_export->m_max;
        const double *mu_a = boa_export->mu_quad_pos;
        double *Tm = (double*)calloc((size_t)n_mu_water * (size_t)n_a * 9,
                                     sizeof(double));
        int ok = (Tm != NULL);
        for (int m = 0; ok && m <= mM; ++m) {
            {
                const int krc = surface_T_aw_coxmunk_fourier_kernel(
                                                    mu_water_pos, n_mu_water,
                                                    mu_a, n_a, m, 128,
                                                    wind_speed, sigma_type,
                                                    n_water, q_convention,
                                                    Tm);
                if (krc != 0) {
                    if (krc == -3) { free(Tm); return -3; }
                    ok = 0; break;
                }
            }
            const double az = (m == 0) ? (2.0 * M_PI) : M_PI;
            const double *Ia = &boa_export->I_per_m[(size_t)m * (size_t)n_a];
            const double *Qa = &boa_export->Q_per_m[(size_t)m * (size_t)n_a];
            const double *Ua = &boa_export->U_per_m[(size_t)m * (size_t)n_a];
            for (int i_w = 0; i_w < n_mu_water; ++i_w) {
                double Iw = 0.0, Qw = 0.0, Uw = 0.0;
                for (int j_a = 0; j_a < n_a; ++j_a) {
                    double wj = az * mu_a[j_a] * w_atm_pos[j_a];
                    if (wj == 0.0) continue;
                    const double *T = Tm + ((size_t)i_w * (size_t)n_a
                                            + (size_t)j_a) * 9;
                    Iw += wj * (T[0]*Ia[j_a] + T[1]*Qa[j_a] + T[2]*Ua[j_a]);
                    Qw += wj * (T[3]*Ia[j_a] + T[4]*Qa[j_a] + T[5]*Ua[j_a]);
                    Uw += wj * (T[6]*Ia[j_a] + T[7]*Qa[j_a] + T[8]*Ua[j_a]);
                }
                coupled->I_inwater_per_m[(size_t)m * (size_t)n_mu_water + i_w] = Iw;
                coupled->Q_inwater_per_m[(size_t)m * (size_t)n_mu_water + i_w] = Qw;
                coupled->U_inwater_per_m[(size_t)m * (size_t)n_mu_water + i_w] = Uw;
            }
        }
        free(Tm);
        if (ok) {
            for (int i_w = 0; i_w < n_mu_water; ++i_w)
                coupled->mu_water_pos[i_w] = mu_water_pos[i_w];
            return 0;
        }
        /* kernel failure: fall through to the legacy hybrid */
    }

    /* Copy water μ grid into output */
    for (int j = 0; j < n_mu_water; ++j) {
        coupled->mu_water_pos[j] = mu_water_pos[j];
    }

    const int m_max     = boa_export->m_max;
    const int n_mu_atm  = boa_export->n_mu;
    const double *mu_atm = boa_export->mu_quad_pos;  /* atm positive μ grid */

    /* Precompute linear interp source arrays per m. For efficiency we
     * pre-extract pointers and do interp inline below. */

    for (int j_w = 0; j_w < n_mu_water; ++j_w) {
        double mu_w = mu_water_pos[j_w];
        if (mu_w <= 0.0 || mu_w > 1.0) {
            /* invalid water μ — fill zero for this column */
            for (int m = 0; m <= m_max; ++m) {
                coupled->I_inwater_per_m[m * n_mu_water + j_w] = 0.0;
                coupled->Q_inwater_per_m[m * n_mu_water + j_w] = 0.0;
                coupled->U_inwater_per_m[m * n_mu_water + j_w] = 0.0;
            }
            continue;
        }

        /* Snell: sin(θ_water) = sin(θ_air) / n_water
         *   sin²(θ_w) = 1 - μ_w²
         *   sin²(θ_air) = n_w² × (1 - μ_w²)
         *   μ_air = sqrt(1 - n_w² × (1 - μ_w²))
         * Validity: μ_air > 0 requires n_w² × (1-μ_w²) < 1, i.e. μ_w² > 1 - 1/n_w².
         *   For n_w=1.34, μ_w_min = sqrt(1 - 1/1.7956) = sqrt(0.4429) ≈ 0.6655.
         *   In-water nodes with μ_w < 0.6655 correspond to *trapped* directions
         *   (TIR cone — these directions cannot receive direct atm transmission;
         *    they only receive in-water scattered light).
         */
        double sin2_water = 1.0 - mu_w * mu_w;
        double sin2_air   = n_water * n_water * sin2_water;
        if (sin2_air >= 1.0) {
            /* Direction inside TIR cone (sin(θ_air) > 1 would be unphysical).
             * No atmospheric radiance reaches this water direction directly.
             * Set inbound = 0 for this column. */
            for (int m = 0; m <= m_max; ++m) {
                coupled->I_inwater_per_m[m * n_mu_water + j_w] = 0.0;
                coupled->Q_inwater_per_m[m * n_mu_water + j_w] = 0.0;
                coupled->U_inwater_per_m[m * n_mu_water + j_w] = 0.0;
            }
            continue;
        }
        double mu_air = sqrt(1.0 - sin2_air);

        /* T_aw radiance Mueller for this (μ_water, n_water).
         * Note: rt_air_water_T_aw takes μ_air (air-side incidence cos).
         * The Mueller matrix already includes the n² radiance enhancement
         * (air radiance n²=1 → water radiance n_w² = 1.7956). */
        double M_T_aw[9];
        rt_air_water_T_aw(mu_air, n_water, q_convention, M_T_aw);

        /* For each m: interp atm BOA Stokes at μ_air, apply T_aw → in-water */
        for (int m = 0; m <= m_max; ++m) {
            const double *I_atm_m = &boa_export->I_per_m[m * n_mu_atm];
            const double *Q_atm_m = &boa_export->Q_per_m[m * n_mu_atm];
            const double *U_atm_m = &boa_export->U_per_m[m * n_mu_atm];

            double I_air = aw_interp_on_unsorted(mu_atm, I_atm_m, n_mu_atm, mu_air);
            double Q_air = aw_interp_on_unsorted(mu_atm, Q_atm_m, n_mu_atm, mu_air);
            double U_air = aw_interp_on_unsorted(mu_atm, U_atm_m, n_mu_atm, mu_air);

            /* Apply T_aw Mueller: L_water = M_T_aw · L_air */
            double I_w = M_T_aw[0*3+0]*I_air + M_T_aw[0*3+1]*Q_air + M_T_aw[0*3+2]*U_air;
            double Q_w = M_T_aw[1*3+0]*I_air + M_T_aw[1*3+1]*Q_air + M_T_aw[1*3+2]*U_air;
            double U_w = M_T_aw[2*3+0]*I_air + M_T_aw[2*3+1]*Q_air + M_T_aw[2*3+2]*U_air;

            coupled->I_inwater_per_m[m * n_mu_water + j_w] = I_w;
            coupled->Q_inwater_per_m[m * n_mu_water + j_w] = Q_w;
            coupled->U_inwater_per_m[m * n_mu_water + j_w] = U_w;
        }
    }

    /* [2026-05-31] FIX-B: for wind>0 the flat-Fresnel Snell coupling above is
     * physically wrong (rough sea surface).  Override the m=0 (azimuthally-
     * averaged) in-water intensity column with the air->water Cox-Munk
     * microfacet BTDF transmission of the m=0 skylight:
     *
     *   L_water^0(μ_w) = ∫_0^1 K^0(μ_w,μ_air) L_air^0(μ_air) μ_air dμ_air
     *   K^0(μ_w,μ_air) = ∫_0^{2π} f_t(μ_air,μ_w,Δφ) dΔφ
     *
     * with f_t = surface_T_aw_coxmunk_btdf_scalar (Walter microfacet, no shadow,
     * includes n² radiance enhancement; ∫f_t cosθ_o dΩ_o = T_aw_direct, validated).
     * A fine internal air-μ GL grid (Nfa) resolves the sharp low-wind lobe; the
     * coarse boa air field L_air^0 is linearly interpolated onto it (rule: no
     * nearest).  Only m=0 I is needed downstream (Ed and the Lu equivalent-beam
     * superposition use I[0] only; coupled Q/U and m>=1 are unused), so m>=1 and
     * Q/U keep their flat values (harmless) and the diffuse-transmitted irradiance
     * Ed is taken in rt_solver from the exact air-side surface_T_aw_coxmunk_direct
     * integral (which equals the hemispheric integral of this field but is
     * computed on the smooth air grid, capturing also the sub-cone leak that the
     * equivalent-beam Lu cannot represent).  Flat-Fresnel is used only at wind<=0. */
    if (wind_speed > 0.0) {
        const int Nfa = 192;   /* fine air-μ nodes (resolves sharp low-wind lobe) */
        const int Nph = 96;    /* azimuth nodes for the m=0 kernel */
        double *famu = (double *)malloc(sizeof(double) * Nfa);
        double *fawt = (double *)malloc(sizeof(double) * Nfa);
        const double *I_atm0 = &boa_export->I_per_m[0 * n_mu_atm];  /* m=0 air radiance */
        if (famu && fawt &&
            rt_quadrature_gauss_legendre_pos(Nfa, famu, fawt) == 0) {
            double dph = 2.0 * M_PI / (double)Nph;
            for (int j_w = 0; j_w < n_mu_water; ++j_w) {
                double mu_w = mu_water_pos[j_w];
                if (mu_w <= 1e-9 || mu_w > 1.0) continue;  /* keep flat-set (0) value */
                double Lw = 0.0;
                for (int a = 0; a < Nfa; ++a) {
                    double mua = famu[a];
                    double Lair0 = aw_interp_on_unsorted(mu_atm, I_atm0,
                                                          n_mu_atm, mua);
                    double K0 = 0.0;
                    for (int p = 0; p < Nph; ++p) {
                        double dphi = dph * (double)p;
                        K0 += surface_T_aw_coxmunk_btdf_scalar(mua, mu_w, dphi,
                                                               n_water, wind_speed,
                                                               sigma_type);
                    }
                    K0 *= dph;                       /* ∫_0^{2π} f_t dΔφ */
                    Lw += fawt[a] * mua * K0 * Lair0;/* ∫_0^1 K0 L_air μ dμ */
                }
                coupled->I_inwater_per_m[0 * n_mu_water + j_w] = Lw;
            }
        }
        free(famu);
        free(fawt);
    }

    return 0;
}

/* ===========================================================================
 * Reverse coupling: in-water z=0⁻ upwelling field → above-water z=0⁺ field on
 * the atmosphere μ-node grid (for injection as the atm SOS bottom source).
 *
 * Mirror of the forward flat path (rt_air_water_couple_atm_to_water):
 *   for each AIR node μ_a:  reverse-Snell μ_a → μ_w,
 *                           interp the 0⁻ field at μ_w (LINEAR; no nearest),
 *                           apply the water→air T_wa Mueller (1/n² upward).
 * The per-mode specular map preserves the azimuth Fourier mode (m_air = m_w),
 * which is exact for a FLAT interface (the flat-Fresnel Mueller has no I↔U /
 * Q↔U coupling, so the cos/cos/sin azimuth reconstruction stays consistent).
 * Rough (wind>0) needs the Cox-Munk water→air BTDF azimuth-mode coupling —
 * that is Milestone 2b and is NOT done here (wind args accepted but unused).
 *
 * Inputs:
 *   {I,Q,U}_water_per_m : z=0⁻ upwelling Fourier field, indexed [m*n_mu_water+j_w]
 *   mu_water_pos        : ascending positive water μ nodes (size n_mu_water)
 *   m_max               : highest mode index (modes 0..m_max)
 *   n_water, q_convention
 * Output (CALLER-allocated, each size (m_max+1)*n_mu_atm, [m*n_mu_atm+j_a]):
 *   {I,Q,U}_air_per_m   : z=0⁺ Fourier field at the air μ-nodes
 *   mu_atm_pos          : ascending positive air μ nodes (size n_mu_atm)
 * Returns 0 on success, negative on error.
 * =========================================================================== */

/* v1.10 OPTION-2: shared m=0 multi-encounter RE-ESCAPE addition, applied on
 * BOTH interface paths (rough kernel path and the flat fallback) just before
 * their success returns.  Energy-conserving on node grids (cell form) or the
 * continuous-density form for weightless single-view callers. */
void mb_re_escape_add_m0(const double *I_water_per_m,
                                const double *mu_water_pos,
                                const double *w_water,
                                int n_mu_water,
                                double n_water, double wind_speed,
                                const double *mu_atm_pos, int n_mu_atm,
                                const double *w_atm_pos,
                                double *I_air_per_m) {
    if (!(wind_speed > 0.0) || !w_water) return;
    /* SPRINT GUARD (2026-07-10): OPT-IN ONLY during the speed sprint
     * (OCRT_MB_CLOSURE=1); default off = pre-surgery bit-identical outputs.
     * Third active site alongside surface.c R_eff and the view-caller term -
     * keep all three gates in sync. */
    {
        static _Thread_local int mb_on = -1;
        if (mb_on < 0) mb_on = (getenv("OCRT_MB_CLOSURE") != NULL);
        if (!mb_on) return;
    }
    double _mb_dbg_sum = 0.0;
enum { MB_NB = 64 };
            double hist[MB_NB];
            for (int j_w = 0; j_w < n_mu_water; ++j_w) {
                double Th_c, Tup_c, Wret_c;
                if (surface_wa_multibounce_cascade(n_water, wind_speed,
                        mu_water_pos[j_w], 300000L,
                        &Th_c, &Tup_c, &Wret_c, hist, NULL, NULL, MB_NB) != 0)
                    continue;
                double E_inc = 2.0 * M_PI * mu_water_pos[j_w] * w_water[j_w]
                             * I_water_per_m[0 * n_mu_water + j_w];
                if (!(E_inc > 0.0)) continue;
                if (w_atm_pos) {
                    /* node-grid form: histogram mass -> nearest node's
                       irradiance cell (exactly energy conserving on the grid) */
                    for (int b = 0; b < MB_NB; ++b) {
                        if (hist[b] <= 0.0) continue;
                        double mu_b = (b + 0.5) / (double)MB_NB;
                        int best = 0; double dmin = 1e9;
                        for (int a2 = 0; a2 < n_mu_atm; ++a2) {
                            double d = fabs(mu_atm_pos[a2] - mu_b);
                            if (d < dmin) { dmin = d; best = a2; }
                        }
                        double cell = 2.0 * M_PI * mu_atm_pos[best] * w_atm_pos[best];
                        if (cell > 0.0)
                            { double _d = E_inc * hist[b] / cell;
                              I_air_per_m[0 * n_mu_atm + best] += _d;
                              _mb_dbg_sum += _d * cell; }
                    }
                } else {
                    /* weightless callers (e.g. the single-VIEW 0+ evaluation):
                       the same physics in its continuous-density form,
                       dL(mu_a) = E_inc * hist_density(mu_a) / (2*pi*mu_a).
                       This is the mathematical continuum limit of the cell
                       form above, not a physics switch. */
                    const double dmu = 1.0 / (double)MB_NB;
                    for (int a2 = 0; a2 < n_mu_atm; ++a2) {
                        double mu_a = mu_atm_pos[a2];
                        int b = (int)(mu_a * MB_NB); if (b >= MB_NB) b = MB_NB - 1;
                        if (hist[b] <= 0.0 || mu_a <= 0.0) continue;
                        { double _d = E_inc * hist[b] / (2.0 * M_PI * mu_a * dmu);
                          I_air_per_m[0 * n_mu_atm + a2] += _d;
                          _mb_dbg_sum += _d; }
                    }
                }
            }
            static _Thread_local int mb_tr = -1;   /* #21: snapshot once */
            if (mb_tr < 0) mb_tr = (getenv("OCRT_MB_TRACE") != NULL);
            if (mb_tr)
                fprintf(stderr, "[CPL] n_atm=%d w_atm=%s added=%.4e\n",
                        n_mu_atm, w_atm_pos ? "Y" : "N", _mb_dbg_sum);
}

/* [P3 / WORKORDER-6, 2026-08-18] per-(mu_a, m) kernel-row cache for the P2b
 * composite grid.  A row depends only on (mu_a, m, wind, sigma_type, n_water,
 * q_convention) plus the compile-time grid constants (keyed too, so future
 * constant changes invalidate).  View cells share vza across raa; caching
 * removes the measured 24x azimuth duplication.  Memory-for-speed per project
 * policy; rows are ~NTOT*9 doubles each. */
typedef struct {
    double mu_a, wind, n_water;
    int m, sigma_type, q_convention, ntot, nphi;
    double *Tm;
} p2b_row_t;
/* DOC-REF: docs/PERSISTENT_COXMUNK_INTERFACE_CACHE_KO_2026-08-25.md
 * Worker-local only: a batch row must never wait for, or consume mutable
 * build state from, another row.  Persistent operator files are immutable and
 * read-only; this hot in-memory row cache remains private to each OpenMP worker. */
static _Thread_local p2b_row_t *g_p2b_rows = NULL;
static _Thread_local int g_p2b_nrows = 0, g_p2b_cap = 0;
static const double *p2b_row_lookup(double mu_a, int m, double wind,
                                    int sigma_type, double n_water,
                                    int q_convention, int ntot, int nphi) {
    for (int i = 0; i < g_p2b_nrows; ++i) {
        const p2b_row_t *r = &g_p2b_rows[i];
        if (r->m == m && r->ntot == ntot && r->nphi == nphi &&
            r->sigma_type == sigma_type && r->q_convention == q_convention &&
            r->mu_a == mu_a && r->wind == wind && r->n_water == n_water)
            return r->Tm;
    }
    return NULL;
}
static void p2b_row_insert(double mu_a, int m, double wind, int sigma_type,
                           double n_water, int q_convention, int ntot,
                           int nphi, const double *Tm) {
    if (g_p2b_nrows == g_p2b_cap) {
        int nc = g_p2b_cap ? 2 * g_p2b_cap : 64;
        p2b_row_t *nr = (p2b_row_t*)realloc(g_p2b_rows,
                                            (size_t)nc * sizeof(p2b_row_t));
        if (!nr) return;               /* cache full: silently skip */
        g_p2b_rows = nr; g_p2b_cap = nc;
    }
    double *cp = (double*)malloc((size_t)ntot * 9 * sizeof(double));
    if (!cp) return;
    memcpy(cp, Tm, (size_t)ntot * 9 * sizeof(double));
    p2b_row_t *r = &g_p2b_rows[g_p2b_nrows++];
    r->mu_a = mu_a; r->wind = wind; r->n_water = n_water;
    r->m = m; r->sigma_type = sigma_type; r->q_convention = q_convention;
    r->ntot = ntot; r->nphi = nphi; r->Tm = cp;
}

int rt_air_water_couple_water_to_atm(const double *I_water_per_m,
                                     const double *Q_water_per_m,
                                     const double *U_water_per_m,
                                     const double *mu_water_pos,
                                     const double *w_water_pos,
                                     int n_mu_water, int m_max,
                                     double n_water, int q_convention,
                                     double wind_speed, int sigma_type,
                                     const double *mu_atm_pos, int n_mu_atm,
                                     const double *w_atm_pos,
                                     double *I_air_per_m,
                                     double *Q_air_per_m,
                                     double *U_air_per_m) {
    if (!I_water_per_m || !Q_water_per_m || !U_water_per_m || !mu_water_pos ||
        !mu_atm_pos || !I_air_per_m || !Q_air_per_m || !U_air_per_m) return -1;
    if (n_water <= 1.0 || n_mu_water < 2 || n_mu_atm < 1 || m_max < 0) return -1;

    const double n2 = n_water * n_water;

    /* [P2 / WORKORDER-5, 2026-08-18] Milestone-2b rough water->air coupling
     * RESTORED.  The former node-sampled variant was disabled because its
     * near-delta transmission rows were under-resolved on the coarse water
     * quadrature; the restored path resolves the lobe on a fine internal
     * water-mu GL grid (Nfw = 192), symmetric to the forward FIX-B fine-grid
     * technique.  wind <= 0 keeps the flat analytic path below untouched
     * (bit-identical requirement, gate G1). */
    if (wind_speed > 0.0) {
        if (getenv("OCRT_DEBUG"))
            fprintf(stderr, "[P2b] water_to_atm rough branch: n_mu_atm=%d m_max=%d\n",
                    n_mu_atm, m_max);
        /* [P2b / WORKORDER-5 delta handling] The reverse transmission row is a
         * near-delta in mu_w centred at the reverse-Snell direction mu_w*(mu_a)
         * (SESSION 14 finding: width ~1e-3 at nadir, narrower than any global
         * quadrature).  Handling: per air direction, integrate on a composite
         * grid = [dense uniform window of NWIN points across mu_w* +/- WHALF]
         * + [global GL(NGLO) with window-interior weights zeroed].  The dense
         * spacing (~3e-4) resolves the lobe core; the global part carries the
         * broad wings incl. the TIR-cone leakage.  Deterministic construction,
         * no per-angle branching. */
        /* [P3c] Original single dense window restored (three-tier variant
         * introduced a composite-trapezoid defect, reverted).  Speed comes
         * from the row cache (24x cell dedup) and n_phi 128 (vs 2048). */
        const int    NWIN  = 1024;
        const int    NGLO  = 96;
        const double WHALF = 0.30;
        const int    NTOT  = NWIN + NGLO;
        const int    n_phi_q = 128; /* [P3] 64 broke G3 (phi-lobe under-resolved at vza>=40); 128 = converged floor (2048 identical) */
        double *xg = (double*)malloc((size_t)NGLO * sizeof(double));
        double *wg = (double*)malloc((size_t)NGLO * sizeof(double));
        double *xi = (double*)malloc((size_t)NTOT * sizeof(double));
        double *wi = (double*)malloc((size_t)NTOT * sizeof(double));
        double *Ff = (double*)malloc((size_t)NTOT * 3 * sizeof(double));
        double *Tm = (double*)malloc((size_t)NTOT * 9 * sizeof(double));
        int ok = (xg && wg && xi && wi && Ff && Tm) &&
                 (rt_quadrature_gauss_legendre_pos(NGLO, xg, wg) == 0);
        for (int j_a = 0; ok && j_a < n_mu_atm; ++j_a) {
            const double mu_a = mu_atm_pos[j_a];
            if (!(mu_a > 0.0) || mu_a > 1.0) {
                for (int m = 0; m <= m_max; ++m) {
                    I_air_per_m[(size_t)m * n_mu_atm + j_a] = 0.0;
                    Q_air_per_m[(size_t)m * n_mu_atm + j_a] = 0.0;
                    U_air_per_m[(size_t)m * n_mu_atm + j_a] = 0.0;
                }
                continue;
            }
            const double sin2_air = 1.0 - mu_a * mu_a;
            const double mu_ws    = sqrt(fmax(0.0, 1.0 - sin2_air / n2));
            double lo = mu_ws - WHALF; if (lo < 1e-6) lo = 1e-6;
            double hi = mu_ws + WHALF; if (hi > 1.0 - 1e-9) hi = 1.0 - 1e-9;
            const double h = (hi - lo) / (double)(NWIN - 1);
            for (int i = 0; i < NWIN; ++i) {
                xi[i] = lo + h * (double)i;
                wi[i] = ((i == 0 || i == NWIN - 1) ? 0.5 : 1.0) * h;
            }
            for (int i = 0; i < NGLO; ++i) {
                xi[NWIN + i] = xg[i];
                wi[NWIN + i] = (xg[i] > lo && xg[i] < hi) ? 0.0 : wg[i];
            }
            surface_fkc_set_build_mmax(m_max);   /* v1.11-speed S1: cap T_wa all-m build */
            for (int m = 0; m <= m_max; ++m) {
                const double *Trow;
                Trow = p2b_row_lookup(mu_a, m, wind_speed, sigma_type,
                                      n_water, q_convention, NTOT, n_phi_q);
                if (!Trow) {
                    const int krc = surface_T_wa_coxmunk_fourier_kernel(
                                                            &mu_a, 1, xi, NTOT,
                                                            m, n_phi_q,
                                                            wind_speed,
                                                            sigma_type,
                                                            n_water,
                                                            q_convention,
                                                            Tm);
                    if (krc != 0) {
                        if (krc == -3) {
                            free(Tm); free(Ff); free(xg); free(wg); free(xi); free(wi);
                            return -3;
                        }
                        ok = 0; break;
                    }
                    p2b_row_insert(mu_a, m, wind_speed, sigma_type, n_water,
                                   q_convention, NTOT, n_phi_q, Tm);
                    Trow = Tm;
                }
                const double az = (m == 0) ? (2.0 * M_PI) : M_PI;
                const double *I_w_m = &I_water_per_m[(size_t)m * (size_t)n_mu_water];
                const double *Q_w_m = &Q_water_per_m[(size_t)m * (size_t)n_mu_water];
                const double *U_w_m = &U_water_per_m[(size_t)m * (size_t)n_mu_water];
                double Ia = 0.0, Qa = 0.0, Ua = 0.0;
                for (int i = 0; i < NTOT; ++i) {
                    if (wi[i] == 0.0) continue;
                    const double wj = az * xi[i] * wi[i];
                    const double *T = Trow + (size_t)i * 9;
                    const double Iw = aw_interp_on_unsorted(mu_water_pos, I_w_m,
                                                            n_mu_water, xi[i]);
                    const double Qw = aw_interp_on_unsorted(mu_water_pos, Q_w_m,
                                                            n_mu_water, xi[i]);
                    const double Uw = aw_interp_on_unsorted(mu_water_pos, U_w_m,
                                                            n_mu_water, xi[i]);
                    Ia += wj * (T[0]*Iw + T[1]*Qw + T[2]*Uw);
                    Qa += wj * (T[3]*Iw + T[4]*Qw + T[5]*Uw);
                    Ua += wj * (T[6]*Iw + T[7]*Qw + T[8]*Uw);
                }
                I_air_per_m[(size_t)m * n_mu_atm + j_a] = Ia;
                Q_air_per_m[(size_t)m * n_mu_atm + j_a] = Qa;
                U_air_per_m[(size_t)m * n_mu_atm + j_a] = Ua;
            }
        }
        free(Tm); free(Ff);
        if (ok) {
            /* mb_re_escape retained per original OPTION-2 design (multi-
             * encounter add-on, both interface paths).  Gate measurement
             * 2026-08-18: on the weightless single-view path this term
             * contributes exactly zero (w_water NULL guard), so it does not
             * interact with the restored kernel there; node-grid double-
             * counting remains an open check (WORKORDER-5 follow-up). */
            mb_re_escape_add_m0(I_water_per_m, mu_water_pos, w_water_pos,
                                n_mu_water, n_water, wind_speed,
                                mu_atm_pos, n_mu_atm, w_atm_pos, I_air_per_m);
            free(xg); free(wg); free(xi); free(wi);
            return 0;
        }
        free(xg); free(wg); free(xi); free(wi);
        /* kernel/alloc failure: fall through to the flat analytic path */
    }

    for (int j_a = 0; j_a < n_mu_atm; ++j_a) {
        double mu_a = mu_atm_pos[j_a];
        /* Reverse Snell: sin²θ_w = sin²θ_a / n².  Every real air direction maps
         * into the water transmission cone, so μ_w is always well-defined
         * (sin²θ_a ≤ 1 ⇒ sin²θ_w ≤ 1/n² < 1). */
        double sin2_air   = 1.0 - mu_a * mu_a;
        double mu_w       = sqrt(fmax(0.0, 1.0 - sin2_air / n2));

        /* Water→air radiance Mueller (includes 1/n² upward enhancement; the
         * routine is TIR-safe and returns zero below the critical angle, which
         * cannot occur for a real air node here). */
        double M_T_wa[9];
        rt_air_water_T_wa(mu_w, n_water, q_convention, M_T_wa);

        for (int m = 0; m <= m_max; ++m) {
            const double *I_w_m = &I_water_per_m[(size_t)m * (size_t)n_mu_water];
            const double *Q_w_m = &Q_water_per_m[(size_t)m * (size_t)n_mu_water];
            const double *U_w_m = &U_water_per_m[(size_t)m * (size_t)n_mu_water];

            double I_w = aw_interp_on_unsorted(mu_water_pos, I_w_m, n_mu_water, mu_w);
            double Q_w = aw_interp_on_unsorted(mu_water_pos, Q_w_m, n_mu_water, mu_w);
            double U_w = aw_interp_on_unsorted(mu_water_pos, U_w_m, n_mu_water, mu_w);

            /* L_air = M_T_wa · L_water (mirror of the forward T_aw application). */
            double I_a = M_T_wa[0*3+0]*I_w + M_T_wa[0*3+1]*Q_w + M_T_wa[0*3+2]*U_w;
            double Q_a = M_T_wa[1*3+0]*I_w + M_T_wa[1*3+1]*Q_w + M_T_wa[1*3+2]*U_w;
            double U_a = M_T_wa[2*3+0]*I_w + M_T_wa[2*3+1]*Q_w + M_T_wa[2*3+2]*U_w;

            I_air_per_m[(size_t)m * (size_t)n_mu_atm + (size_t)j_a] = I_a;
            Q_air_per_m[(size_t)m * (size_t)n_mu_atm + (size_t)j_a] = Q_a;
            U_air_per_m[(size_t)m * (size_t)n_mu_atm + (size_t)j_a] = U_a;
        }
    }
    mb_re_escape_add_m0(I_water_per_m, mu_water_pos, w_water_pos,
                        n_mu_water, n_water, wind_speed,
                        mu_atm_pos, n_mu_atm, w_atm_pos, I_air_per_m);
    return 0;
}
