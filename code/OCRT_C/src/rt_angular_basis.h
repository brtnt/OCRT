#ifndef OCRT_RT_ANGULAR_BASIS_H
#define OCRT_RT_ANGULAR_BASIS_H

/* Build the scalar and real spin-2 angular bases used by the OCRT
 * azimuthal-Fourier representation.  Direction-indexed rows are offset so
 * that the second index is valid on [-n_mu,+n_mu]. */
int rt_angular_basis_build_scalar(int n_mu,
                                  int l_max,
                                  int fourier_order,
                                  const double *mu,
                                  double **scalar);

int rt_angular_basis_build_spin2(int n_mu,
                                 int l_max,
                                 int fourier_order,
                                 const double *mu,
                                 double **linear_r,
                                 double **linear_t);

/* Evaluate the same normalized scalar/spin-2 basis at one arbitrary signed
 * direction cosine.  These routines are used by PSSA reflected-beam sources,
 * whose local incident direction varies with altitude and is not a quadrature
 * node.  Output arrays have length l_max+1. */
int rt_angular_basis_eval_scalar(int l_max,
                                 int fourier_order,
                                 double mu,
                                 double *scalar);
int rt_angular_basis_eval_spin2(int l_max,
                                int fourier_order,
                                double mu,
                                double *linear_r,
                                double *linear_t);

/* Contract precomputed bases with beta/gamma/alpha/zeta phase moments. */
void rt_angular_kernel_scalar(int n_mu,
                              int l_max,
                              int fourier_order,
                              double **scalar,
                              const double *beta,
                              double **out);

void rt_angular_kernel_cross(int n_mu,
                             int l_max,
                             int fourier_order,
                             double **scalar,
                             double **linear_r,
                             double **linear_t,
                             const double *gamma,
                             double **out_r,
                             double **out_t);

void rt_angular_kernel_linear(int n_mu,
                              int l_max,
                              int fourier_order,
                              double **linear_r,
                              double **linear_t,
                              const double *alpha,
                              const double *zeta,
                              double **out_rr,
                              double **out_rt,
                              double **out_tt);

#endif
