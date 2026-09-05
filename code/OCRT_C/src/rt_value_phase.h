#ifndef OCRT_RT_VALUE_PHASE_H
#define OCRT_RT_VALUE_PHASE_H

#include "rt_kernel.h"
#include "rt_types.h"

/*
 * Direct particle phase table used by the value-space Fourier/SOS path.
 *
 * Stage-3B contract:
 *   - theta_deg is strictly ascending and spans 0..180 deg;
 *   - P11/P12/P33 are normalized together so 0.5*int(P11 dmu)=1;
 *   - interpolation is piecewise linear in the native scattering angle theta;
 *   - one common bracket/weight is used for P11/P12/P33, preserving
 *       P11>=0, |P12|<=P11 and |P33|<=P11 between valid input nodes;
 *   - FR631 tables use the O(1) affine-index fast path.  Other monotone
 *     legacy tables use a binary-search theta-linear fallback.
 */
typedef struct {
    int n;
    int is_fr631;
    double *theta_deg;   /* ascending, [n] */
    double *p11;         /* normalized, [n] */
    double *p12;         /* normalized, [n] */
    double *p33;         /* normalized, [n] */
    double norm_before;  /* 0.5*int(P11 sin(theta)dtheta) before normalization */
    double g_asym;       /* exact integral of the piecewise-linear normalized table */
    double bb_b_ratio;   /* backward / total scattering integral */
} rt_value_phase_interp_t;

int rt_value_phase_interp_build(rt_value_phase_interp_t *out,
                                const double *theta_deg,
                                const double *p11,
                                const double *p12,
                                const double *p33,
                                int n,
                                int norm_n_mu /* retained for source compatibility; ignored */);
void rt_value_phase_interp_free(rt_value_phase_interp_t *tab);

/* OSOAA-compatible real-angle hydrosol forward-peak truncation.
 *
 * The input table is already normalized and represented piecewise-linearly
 * in theta.  P11 is replaced below theta2=acos(mu2) by the log10-linear
 * continuation defined by the exact, theta-linearly interpolated values at
 * theta1=acos(mu1) and theta2.  P12/P33 are multiplied by the same local
 * P11 ratio, so polarization ratios and pointwise physical bounds are
 * preserved.  The returned table is normalized again for direct Fourier/SOS
 * use.  A_out is OSOAA A_TRONCA=2*f, where f is the removed forward fraction;
 * the matching transport scattering coefficient is b_eff=b*(1-A/2).
 *
 * If A is non-positive or below threshold_A, the operation is an exact
 * no-op (A_out=0 and out is a normalized copy of src).
 */
int rt_value_phase_interp_loglinear_truncate(
    const rt_value_phase_interp_t *src,
    double mu1, double mu2, double threshold_A,
    rt_value_phase_interp_t *out,
    double *A_out);

/* Evaluate with mu=cos(theta). */
void rt_value_phase_interp_eval(const rt_value_phase_interp_t *tab,
                                double mu,
                                double *p11,
                                double *p12,
                                double *p33);

/* Evaluate directly in degrees. */
void rt_value_phase_interp_eval_theta(const rt_value_phase_interp_t *tab,
                                      double theta_deg,
                                      double *p11,
                                      double *p12,
                                      double *p33);

/* Evaluate from a precomputed lower index and common linear weight.  This is
 * valid for an ascending table and is used by geometry-map caches. */
void rt_value_phase_interp_eval_cached(const rt_value_phase_interp_t *tab,
                                       int lower,
                                       double weight,
                                       double *p11,
                                       double *p12,
                                       double *p33);

/* Direct azimuth-space vector phase kernel.  Output layout is
 * [m][j][k+n_mu], j=0..n_mu, k=-n_mu..+n_mu for each of the six tables.
 * All Fourier modes are accumulated in one (j,k,phi) sweep so the phase
 * interpolation is performed only once per angle sample. */

/* Scalar P11 value kernel, all Fourier modes in one azimuth sweep. */
int rt_value_phase_fourier_scalar_allm(
    const rt_atm_t *atm, int m_count, const rt_value_phase_interp_t *phase,
    int nphi, double *pfm);

/* Recompute an inclusive incoming-row interval [j_begin,j_end] in existing
 * all-mode arrays.  Used by the water beam cache: GL rows are immutable,
 * while row 0 follows the current incident-beam cosine. */
int rt_aerosol_value_phase_fourier_pol_rows(
    const rt_atm_t *atm, int m_count, const rt_value_phase_interp_t *phase,
    int nphi, int j_begin, int j_end,
    double *pfm, double *gr, double *gt, double *arr,
    double *art, double *att);

int rt_aerosol_value_phase_fourier_pol_allm(
    const rt_atm_t *atm,
    int m_count,
    const rt_value_phase_interp_t *phase,
    int nphi,
    double *pfm,
    double *gr,
    double *gt,
    double *arr,
    double *art,
    double *att);

/* Recompute only entries coupled to the beam slot rm[0]: incoming row j=0
 * and outgoing column k=0, for all Fourier modes.  The GL/view block is
 * beam-independent and may be held in the dedicated value-kernel cache. */
int rt_aerosol_value_phase_fourier_pol_solar_allm(
    const rt_atm_t *atm,
    int m_count,
    const rt_value_phase_interp_t *phase,
    int nphi,
    double *pfm,
    double *gr,
    double *gt,
    double *arr,
    double *art,
    double *att);

/* Exact response of the direct particle phase matrix to an incoming Stokes-Q
 * solar component.  Output layout is [m][k+n_mu] for I<-Q, Q<-Q and U<-Q.
 * This column is not completely represented by the standard six diffuse
 * kernel tables. */
int rt_aerosol_value_phase_fourier_pol_beamq_allm(
    const rt_atm_t *atm,
    int m_count,
    const rt_value_phase_interp_t *phase,
    int nphi,
    double *i_from_q,
    double *q_from_q,
    double *u_from_q);

/* Convenience copy of one mode into an allocated Legendre workspace. */
int rt_aerosol_value_phase_fourier_pol(
    rt_legendre_workspace_t *ws,
    const rt_atm_t *atm,
    int m,
    const rt_value_phase_interp_t *phase,
    int nphi);

#endif
