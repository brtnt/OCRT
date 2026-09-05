#ifndef OCRT_V2_RT_QUADRATURE_H
#define OCRT_V2_RT_QUADRATURE_H

/* Single-Gauss quadrature for plane-parallel radiative transfer.
 *
 * Computes a (2 * n_mu)-point Gauss-Legendre rule on [-1, 1] and returns
 * the n_mu positive nodes/weights only. The n_mu negative nodes are the
 * mirror -mu[i] with the same weights w[i] (caller fills these into the
 * full rm/gb arrays in rt_atm).
 *
 * Output (caller-allocated, length n_mu):
 *   mu[i] in (0, 1), ascending      mu[0] smallest .. mu[n_mu-1] largest
 *   w[i]  positive,  Σ w[i] = 1     (half of full [-1,1] sum, by symmetry)
 *
 * Polynomial integration (the full 2*n_mu rule on [-1, 1] is exact for
 * degree ≤ 4*n_mu - 1):
 *   Σ_pos w[i] * mu[i]^k = 1/(k+1)   for even k = 0, 2, ..., 4*n_mu - 2
 * Odd-k sums on the positive half do NOT match 1/(k+1) — that integral
 * is recovered only when both halves are summed (and is zero by symmetry).
 *
 * Returns 0 on success, -1 on bad args, -2 if Newton iteration failed.
 */
int rt_quadrature_gauss_legendre_pos(int n_mu, double *mu, double *w);

#endif /* OCRT_V2_RT_QUADRATURE_H */
