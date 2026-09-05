/* shared/mat3.h
 *
 * 3x3 Mueller matrix operations and Hovenier (1983) rotation between
 * meridian and scattering-plane bases.
 *
 * Convention: row-major flat array of 9 doubles, M[i*3+j].
 * Stokes 3-vector: [I, Q, U] (no V circular polarization).
 *
 * compute_hovenier_rotation produces (cos, sin) of rotation angles
 * i1, i2 such that:
 *   Z_meridian = L(i2) × P_scattering × L(i1)
 * where L is the rotation matrix (build_rotation_L).
 *
 * Both vrt_solver and ocean_solver include this file.
 */
#ifndef OCRT_SHARED_MAT3_H
#define OCRT_SHARED_MAT3_H

#include <math.h>

static inline void mat3_zero(double* M) {
    for (int i = 0; i < 9; ++i) M[i] = 0.0;
}

static inline void mat3_mul(const double* A, const double* B, double* out) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double s = 0.0;
            for (int k = 0; k < 3; ++k) {
                s += A[i*3+k] * B[k*3+j];
            }
            out[i*3+j] = s;
        }
    }
}

/* Matrix-vector product: out[3] = M[3x3] * v[3] (vrt only) */
static inline void mat3_vec(const double* M, const double* v, double* out) {
    out[0] = M[0]*v[0] + M[1]*v[1] + M[2]*v[2];
    out[1] = M[3]*v[0] + M[4]*v[1] + M[5]*v[2];
    out[2] = M[6]*v[0] + M[7]*v[1] + M[8]*v[2];
}

/* Build rotation matrix L(i) = diag(1, cos(2i), -sin(2i); 0, sin(2i), cos(2i))
 * for Stokes 3-vector basis rotation. Uses cos(i), sin(i) inputs. */
static inline void build_rotation_L(double ci, double si, double* L) {
    double c2 = ci * ci - si * si;       /* cos(2i) */
    double s2 = 2.0 * ci * si;            /* sin(2i) */
    mat3_zero(L);
    L[0*3+0] = 1.0;
    L[1*3+1] =  c2;
    L[1*3+2] = -s2;
    L[2*3+1] =  s2;
    L[2*3+2] =  c2;
}

/* Hovenier (1983) rotation angles i1, i2 between meridian and scattering plane.
 * Inputs:
 *   mu_out, mu_in: cosines of out/in zenith
 *   s1, s2: sines (= sqrt(1 - mu²))
 *   cT, sT: cosine and sine of scattering angle
 *   sin_phi: sin of relative azimuth (signed)
 * Outputs: cos, sin of i1 and i2 (normalized, numerically stable at sT→0).
 */
static inline void compute_hovenier_rotation(double mu_out, double mu_in,
                                                 double s1, double s2,
                                                 double cT, double sT,
                                                 double sin_phi,
                                                 double* ci1, double* si1,
                                                 double* ci2, double* si2)
{
    double denom1 = fmax(1e-12, s2 * sT);
    double denom2 = fmax(1e-12, s1 * sT);
    double sT_safe = fmax(1e-12, sT);
    double c1 = (mu_out - mu_in * cT) / denom1;
    double s1r = s1 * sin_phi / sT_safe;
    double n1 = sqrt(c1 * c1 + s1r * s1r);
    if (n1 < 1e-12) { *ci1 = 1.0; *si1 = 0.0; }
    else            { *ci1 = c1 / n1; *si1 = s1r / n1; }
    double c2 = (mu_in - mu_out * cT) / denom2;
    double s2r = s2 * sin_phi / sT_safe;
    double n2 = sqrt(c2 * c2 + s2r * s2r);
    if (n2 < 1e-12) { *ci2 = 1.0; *si2 = 0.0; }
    else            { *ci2 = c2 / n2; *si2 = s2r / n2; }
}

/* WORKORDER_1 (2026-08-18): exact analytic ±1 pole limits, cos_phi-aware.
 * Identical to compute_hovenier_rotation on every non-degenerate node (the
 * generic branch below is a verbatim copy), plus the exact limits at the
 * exact poles.  Pole nodes are inserted into rm[] as exact ±1.0, so
 * s = sqrt(fmax(0, 1-mu*mu)) == 0.0 exactly and the branch tests are exact.
 *
 * Derivation (outgoing pole, s1==0, mu_out=±1): cT = mu_out*mu_in exactly and
 * sT = s2; expanding the generic c2 numerator/denominator to first order in
 * s1 gives  lim cos(i2) = -mu_out*cos_phi,  lim sin(i2) = sin_phi  (unit norm
 * automatic).  The incoming pole (s2==0, mu_in=±1) is symmetric:
 * lim cos(i1) = -mu_in*cos_phi, lim sin(i1) = sin_phi.  Consequently at
 * mu_out=+1 the Z rows carry cos(2*i2)=... = pure cos/sin(2*phi) structure,
 * so I keeps only m=0 and Q/U keep only m=2 after Fourier projection
 * (discrete orthogonality exact on the uniform nphi grid).  At a double pole
 * P12(0 or 180 deg)=0 removes any Q/U ambiguity. */
static inline void compute_hovenier_rotation_cp(double mu_out, double mu_in,
                                                 double s1, double s2,
                                                 double cT, double sT,
                                                 double cos_phi, double sin_phi,
                                                 double* ci1, double* si1,
                                                 double* ci2, double* si2)
{
    if (s2 == 0.0) { *ci1 = -mu_in * cos_phi; *si1 = sin_phi; }
    else {
        double denom1 = fmax(1e-12, s2 * sT);
        double sT_safe = fmax(1e-12, sT);
        double c1 = (mu_out - mu_in * cT) / denom1;
        double s1r = s1 * sin_phi / sT_safe;
        double n1 = sqrt(c1 * c1 + s1r * s1r);
        if (n1 < 1e-12) { *ci1 = 1.0; *si1 = 0.0; }
        else            { *ci1 = c1 / n1; *si1 = s1r / n1; }
    }
    if (s1 == 0.0) { *ci2 = -mu_out * cos_phi; *si2 = sin_phi; }
    else {
        double denom2 = fmax(1e-12, s1 * sT);
        double sT_safe = fmax(1e-12, sT);
        double c2 = (mu_in - mu_out * cT) / denom2;
        double s2r = s2 * sin_phi / sT_safe;
        double n2 = sqrt(c2 * c2 + s2r * s2r);
        if (n2 < 1e-12) { *ci2 = 1.0; *si2 = 0.0; }
        else            { *ci2 = c2 / n2; *si2 = s2r / n2; }
    }
}

#endif /* OCRT_SHARED_MAT3_H */
