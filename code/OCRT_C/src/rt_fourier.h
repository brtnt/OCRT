#ifndef OCRT_RT_FOURIER_H
#define OCRT_RT_FOURIER_H

/* OCRT azimuthal Fourier convention.
 *
 * The public OCRT RAA is defined operationally in rt_raa_convention.h
 * (direct-glint branch at RAA=180 deg).  The stored atmospheric modal
 * coefficients use the photon-propagation azimuth, which differs from the
 * public reconstruction argument by pi.  For a real Stokes field the negative
 * Fourier orders are the complex conjugates of the positive orders, so the
 * real expansion can be folded onto m >= 0:
 *
 *   X(phi) = X_0 + 2 sum_{m=1}^M X_m cos[m(phi + pi)]  for I and Q,
 *   U(phi) =       2 sum_{m=1}^M U_m sin[m(phi + pi)]  for U.
 *
 * When a particular operator is explicitly represented in a coordinate whose
 * azimuth is shifted by pi, its modal coefficient acquires exp(i m pi)=(-1)^m.
 * This is a coordinate-conversion helper, not a universal reflection factor.
 * In particular, water-side specular internal reflection preserves the shared
 * physical azimuth and must not use this sign.
 */

double rt_fourier_reconstruct_cos(const double *mode_values,
                                  int m_max,
                                  double public_relative_azimuth_rad);

double rt_fourier_reconstruct_sin(const double *mode_values,
                                  int m_max,
                                  double public_relative_azimuth_rad);

double rt_fourier_pi_shift_sign(int m);

#endif /* OCRT_RT_FOURIER_H */
