/* ============================================================================
 * rt_ipss.h — IPSS (improved pseudo-spherical shell) correction
 * ============================================================================
 * Zhai, P.-W., Hu, Y. (2022). An improved pseudo spherical shell algorithm for
 * vector radiative transfer.  JQSRT 282, 108132.  doi:10.1016/j.jqsrt.2022.108132
 *
 * WHY THIS REPLACED THE LEGACY PSSA PATH
 *   The legacy `--pssa` (rt_pssa.c; average-secant Chapman, Dahlback & Stamnes
 *   1991 / He et al. 2018) corrected ONLY the solar direct-beam attenuation and
 *   evaluated it on the NADIR COLUMN above the target level.  At slant view it
 *   therefore used the wrong solar-beam geometry: published error at
 *   SZA 84.26 deg is -16%..+5% over azimuth (paper Sec. 3), and the error grows
 *   with VZA.  Geostationary ocean-colour geometry (large fixed VZA x high SZA
 *   slots) sits exactly on that weakness.
 *
 *   v1.11 (2026-09-05): the legacy source files, the algorithm-selector field
 *   in rt_options_t and the --pssa-mode CLI option were all DELETED.  `--pssa`
 *   is a plain on/off switch and IPSS is the only spherical algorithm in the
 *   tree, so no run can enter the old path and rt_atm_t is never modified by a
 *   spherical correction.
 *
 *   IPSS instead (a) solves single scattering with the true spherical-shell
 *   geometry along the actual line of sight and (b) transports the
 *   multiple/single ratio from the plane-parallel solution:
 *
 *       I_ipss = (I_pp / I_1,pp) * I_1,ss ,  Q,U scaled by the same factor
 *              = kappa * I_pp,   kappa = I_1,ss / I_1,pp                (Eq. 7)
 *
 *   Because the scattering angle is constant along a straight line of sight and
 *   the sun beam is parallel, E0 and the phase function of ONE scatterer are
 *   identical in numerator and denominator.  They cancel in kappa only if the
 *   Rayleigh/aerosol mixing ratio is the same in every layer.  It is not in the
 *   production atmosphere (aerosol scale height 2 km under Rayleigh ~8 km), so
 *   the layer source must be weighted by the real single-scattering phase
 *   functions, w_j = w_ray_j P_R(Theta) + w_aer_j P_A(Theta) (v1.11.1; see the
 *   PHASE-WEIGHTED SOURCE note in rt_ipss_t).  Up to v1.11 the weight was the
 *   phase-free omega_j, which over-weighted the low aerosol layers by
 *   P_R/P_A ~ 8-10 at side/back scattering: -0.07 %p at low SZA and up to
 *   -1.9 %p at SZA 85 (555 nm, VZA 55, AOD 0.1) relative to the paper's ratio.
 *
 * THE SOS SOLVE STAYS PLANE-PARALLEL
 *   Eq. 7 requires a PLANE-PARALLEL I_pp, so nothing may perturb the SOS solve
 *   itself.  Since v1.11 that is structural rather than conditional: no
 *   spherical code touches rt_atm_t at all, and the correction is applied only
 *   to the finished Stokes vector.
 *
 * VZA ANCHOR (decides the sign AND size of the correction)
 *   A plane-parallel solve has no notion of WHERE its mu_v is the local zenith
 *   angle — the direction is the same at every altitude.  A spherical ray's
 *   local zenith angle is not: sin(theta(r)) * r = const, so theta falls with
 *   altitude (VZA 55 deg at the pixel is 53.755 deg at 100 km).  The anchor
 *   therefore has to be chosen, and kappa depends on it:
 *
 *     SURFACE anchor (default, and the only one used in production):
 *       vza is the zenith angle AT THE TARGET PIXEL — satellite L1B senz, and
 *       the same angle the plane-parallel solve used for mu_v.  The two share
 *       one geometry at the level where the atmosphere's mass sits, so the
 *       correction goes to zero in every plane-parallel limit.  Measured with
 *       an 8 km exponential profile, tau = 0.0935, sun at the zenith:
 *         VZA 55 -> +0.24 %, VZA 70 -> +0.83 % (and -> 0 as the scale height
 *         collapses: +0.016 % at 0.5 km, +0.005 % at 0.1 km).  Gate G-I3.
 *
 *     TOA anchor (retained only to reproduce the paper's Sec. 3 cases):
 *       vza is the zenith angle at 100 km, with the sensor on the nadir column.
 *       The ray then runs 1.3 deg steeper than mu_v through the dense lower
 *       atmosphere, so a spurious, SZA-INDEPENDENT offset appears — -2.9 % at
 *       VZA 55 and -11 % at VZA 70 even with the sun at the zenith — and it
 *       does NOT vanish as the atmosphere collapses onto the surface
 *       (-3.0 % at a 0.1 km scale height).  Wrong limit => wrong for OCRT.
 *
 *   A second reason: kappa's phase matrix cancels only if the plane-parallel
 *   scattering angle equals the true one.  Theta is fixed by the two space
 *   directions, and OCRT's --sza/--raa are pixel-referenced, so only the
 *   SURFACE anchor keeps the view angle in the same reference frame.
 *
 *   Side benefit: the SURFACE anchor has no TOA-horizon limit (79.86 deg for
 *   H = 100 km); every vza < 90 deg has a ground-connected path.
 *
 * IMPLEMENTATION NOTES (from the reference package, 2026-08-25)
 *   A. Paper Eq. (5c) prints the sine-rule radius as R_e; the exact form is
 *      (R_e + h_i).  Uniform-density validation hides this by telescoping, but
 *      a real (non-uniform) profile shifts tau_view by ~1.3%.  This module uses
 *      a parametric ray-shell tracer, which is exact by construction.
 *   B. Paper Eq. (6b) prints theta'_s = acos(+r_B.V0/|r_B|); only the LOCAL
 *      solar zenith angle acos(-r_B.V0/|r_B|) is self-consistent with (6c)-(6e).
 *      Not an issue here either: the tracer never uses the closed form.
 *   C. The solar ray MUST terminate at the first ground intersection
 *      (shadow test).  Omitting it double-counts the far-side atmosphere
 *      (exactly 2x error).  ipss_ray_segments() reports hit_ground.
 *
 * SCOPE (Phase I, per IPSS_IMPLEMENTATION_WORKORDER_v1.0 Sec. 3)
 *   Atmosphere-only / black-ocean TOA radiance.  The coupled ocean-atmosphere
 *   application of Eq. 7 is an open design decision and is refused fail-loud.
 * ========================================================================== */
#ifndef RT_IPSS_H
#define RT_IPSS_H

#include "rt_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Earth radius used by the paper's validation cases; kept identical to the
 * value the deleted legacy path used so historical geometry is comparable. */
#define RT_IPSS_EARTH_RADIUS_KM 6371.0

/* Gauss-Legendre nodes PER LAYER PANEL.  The line-of-sight integral is done
 * panel by panel (one panel per shell crossing) because beta jumps at every
 * shell boundary: a single global rule over the whole path wobbles ~1% between
 * 32 and 512 nodes on a realistic US62 grid, while the panel rule is converged
 * at 2 nodes (1 -> 16 moves the answer by <2e-9).  Gate G-I1b refines it. */
#define RT_IPSS_N_QUAD_DEFAULT 2

typedef enum {
    RT_IPSS_SUN_SPHERICAL = 0,  /* tau_sun at the true line-of-sight point (IPSS) */
    RT_IPSS_SUN_NADIR     = 1,  /* tau_sun on the nadir column at the same altitude
                                 * (== what legacy PSS effectively assumes) */
    RT_IPSS_SUN_PLANEPAR  = 2   /* tau_sun = tau(h)/mu_s (plane-parallel) */
} rt_ipss_sun_mode_t;

/* Where the input VZA is the local zenith angle.  See "VZA ANCHOR" above. */
typedef enum {
    RT_IPSS_VZA_ANCHOR_SURFACE = 0, /* at the target pixel (L1B senz) — production */
    RT_IPSS_VZA_ANCHOR_TOA     = 1  /* at TOA on the nadir column — paper Sec. 3 */
} rt_ipss_vza_anchor_t;

typedef struct {
    int     nlay;        /* number of layers */
    double *h_edge;      /* [nlay+1] ASCENDING altitude [km]; h_edge[0] = 0 (BOA) */
    double *beta;        /* [nlay] extinction [km^-1], index 0 = bottom layer */
    double *omega;       /* [nlay] single-scattering albedo of the layer */
    double  Re;          /* Earth radius [km] */
    double  H;           /* TOA altitude [km] = h_edge[nlay] */
    double  mu_sun;      /* cos(SZA) */
    int     n_quad;      /* line-of-sight quadrature nodes */
    /* Scratch owned by the object so the hot path allocates nothing.  One
     * rt_ipss_t belongs to one solve; the LUT grid loop is serial and the
     * batch driver gives every worker its own solve, so no sharing occurs. */
    double *scr_seg;     /* [nlay]   ray-segment lengths */
    double *scr_ts;      /* [2*nlay+4] shell-crossing parameters */
    double *tau_above;   /* [nlay+1] vertical optical depth above each edge */
    double *gl_x;        /* [n_quad] Gauss-Legendre nodes/weights, cached */
    double *gl_w;
    int     gl_n;        /* size the cache was built for */
    /* Line-of-sight crossing table, rebuilt once per viewing geometry so that
     * tau_view at a quadrature node is a binary search instead of a full trace. */
    double *los_t;       /* [2*nlay+8] crossing parameters, ascending */
    double *los_tau;     /* [2*nlay+8] cumulative tau at each crossing */
    int    *los_lay;     /* [2*nlay+8] layer index of the interval that starts there */
    int     los_n;
    /* PHASE-WEIGHTED SOURCE (v1.11.1, 2026-09-05).  Eq. (7) divides the real
     * single-scattering radiances, whose layer source is
     *   beta_j * [ w_ray_j * P_R(Theta) + w_aer_j * P_A(Theta) ]
     * (w = scattering fraction of the layer extinction).  Under a parallel
     * beam and a straight line of sight the scattering angle Theta is the same
     * at every point, so the phase functions cancel in kappa ONLY when the
     * Rayleigh/aerosol mixing ratio is height-independent.  In the production
     * atmosphere it is not (aerosol H = 2 km under Rayleigh ~8 km), and at
     * side/back-scattering angles P_A/P_R ~ 0.1 makes the Rayleigh layers
     * dominate the single scattering; the phase-free weight omega_j (used up to
     * v1.11) over-weighted the low aerosol layers.  RTSOS builds the same
     * ratio from its full phase matrices (rtsos_rao_dg.f90, PHASE_MATRIX_CALC
     * inside the spherical single-scatter recursion).  NULL w_ray/w_aer =
     * phase-free weights omega_j (single-component profiles, where P cancels
     * exactly; the G-I gates use these). */
    double *w_ray;       /* [nlay] Rayleigh scattering fraction of extinction */
    double *w_aer;       /* [nlay] aerosol  scattering fraction of extinction */
    double  depol;       /* Rayleigh depolarization factor */
    int     L_aer;       /* aerosol P11 Legendre order (<0: no moments) */
    double *betal_aer;   /* [L_aer+1] P11 moments, beta_0 = 1 (owned copy) */
    int     n_ang_aer;   /* tabulated aerosol P11 (preferred: Gibbs-free, what
                          * the production value kernel scatters with);
                          * pointers BORROWED from the atm value-kernel table */
    const double *theta_aer;   /* [n_ang_aer] ascending degrees */
    const double *p11_aer;     /* [n_ang_aer] normalized so (1/2) int P11 dmu = 1 */
    double *w_src;       /* [nlay] effective source weight for the current
                          * scattering angle (rt_ipss_set_source_weights) */
    double  cos_theta;   /* scattering angle the weights were built for */
    double  P_R, P_A;    /* the two phase-function values at that angle */
} rt_ipss_t;

/* Build w_src for one viewing geometry (no-op when w_ray/w_aer are NULL).
 * rt_ipss_kappa() calls it; call it explicitly before using the single-
 * scatter primitives on a phase-weighted profile with a new (vza, raa). */
void rt_ipss_set_source_weights(rt_ipss_t *G, double vza_deg, double raa_deg);
/* Normalized depolarized Rayleigh phase function and the aerosol P11 at a
 * given cos(scattering angle), as used by rt_ipss_set_source_weights(). */
double rt_ipss_phase_rayleigh(double depol, double cos_theta);
double rt_ipss_phase_aerosol(const rt_ipss_t *G, double cos_theta);
/* cos(Theta) for OCRT (sza, vza, raa): solar propagation V0 = (sin ts, 0, -cos ts),
 * view u = (sin tv cos fp, sin tv sin fp, cos tv), fp = 180 deg - raa. */
double rt_ipss_cos_scattering_angle(double mu_sun, double vza_deg, double raa_deg);

/* Build the IPSS geometry view of an already-populated rt_atm_t.
 * OCRT stores levels TOA-first (z[0] = TOA, z[nt] = 0) with cumulative tau
 * h[0] = 0 .. h[nt] = tau_total; this reverses them into ascending shells.
 * Returns 0, -1 on bad input, -2 if the altitude grid is not strictly
 * decreasing in k (zero-thickness layers cannot be ray-traced). */
int  rt_ipss_init_from_atm(rt_ipss_t *G, const rt_atm_t *atm, double Re_km);

/* Allocate the scratch/memo buffers and build the vertical-tau cache for a
 * struct whose nlay / h_edge / beta / omega / Re / H / mu_sun are already set.
 * rt_ipss_init_from_atm() calls it; tests and tools that construct a profile
 * by hand must call it too (the hot path assumes the scratch exists). */
int  rt_ipss_finalize(rt_ipss_t *G);
void rt_ipss_free(rt_ipss_t *G);

/* Eq. (7) correction factor for one viewing geometry.
 *   vza_deg : viewing zenith angle AT THE TARGET PIXEL (SURFACE anchor; the
 *             satellite L1B senz convention and the same angle the
 *             plane-parallel solve used for mu_v).
 *   raa_deg : OCRT public relative azimuth.  Internally mapped to the paper's
 *             frame by phi_paper = 180 deg - raa (verified against the OCRT
 *             direct-glint half-plane and the first-order phase dependence).
 * Outputs kappa = I_1,ss / I_1,pp plus both integrals for diagnostics.
 * *guard is set to 1 when I_1,pp underflows the guard threshold and kappa is
 * forced to 1 (no correction) instead of dividing by ~0. */
int rt_ipss_kappa(const rt_ipss_t *G, double vza_deg, double raa_deg,
                  double *kappa, double *I1_ss, double *I1_pp, int *guard);

/* Single-scatter line integral along the spherical line of sight with the
 * requested solar-attenuation convention.  Phase function and E0 are omitted
 * (they cancel in every ratio this module forms).  Exposed for gate G-I1.
 *   rt_ipss_single_scatter()      : TOA-anchored vza (paper Sec. 3 geometry).
 *   rt_ipss_single_scatter_surf() : pixel-anchored vza (production; what
 *                                   rt_ipss_kappa() uses). */
double rt_ipss_single_scatter(const rt_ipss_t *G, double vza_deg, double raa_deg,
                              rt_ipss_sun_mode_t mode);
double rt_ipss_single_scatter_surf(const rt_ipss_t *G, double vza_deg,
                                   double raa_deg, rt_ipss_sun_mode_t mode);

/* Fully plane-parallel single-scatter integral (same normalization). */
double rt_ipss_single_scatter_pp(const rt_ipss_t *G, double vza_deg);

/* Exact spherical slant optical depth of the solar beam from TOA down to the
 * ground point directly below the sensor; exp(-tau) is the direct-beam
 * transmittance that the legacy path approximated with the Chapman secant.
 * Returns a negative value if the point is not illuminated. */
double rt_ipss_tau_dir_dn_boa(const rt_ipss_t *G);

/* Geometry primitives (exposed for the G-I0 gate against the Python oracle). */
int    rt_ipss_ray_segments(const rt_ipss_t *G, const double r0[3],
                            const double u[3], double *seg, int *hit_ground);
double rt_ipss_tau_view_toa(const rt_ipss_t *G, double vza_deg);
double rt_ipss_tau_sun_at(const rt_ipss_t *G, const double rB[3]);

#ifdef __cplusplus
}
#endif
#endif /* RT_IPSS_H */
