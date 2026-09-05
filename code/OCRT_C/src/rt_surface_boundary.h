#ifndef OCRT_RT_SURFACE_BOUNDARY_H
#define OCRT_RT_SURFACE_BOUNDARY_H

#include "rt_kernel.h"
#include "rt_solver.h"

/* Fourier-mode representation of the lower atmospheric reflection operator.
 * The kernel is immutable after initialization and can be reused by every SOS
 * scattering order for the same angular grid and Fourier order. */
typedef struct {
    int n_mu;
    int fourier_order;
    int is_flat;
    double *mu_positive;
    double *kernel_3x3;  /* [outgoing][incoming][3][3] */
} rt_surface_boundary_t;

int rt_surface_boundary_init(rt_surface_boundary_t *boundary,
                             const rt_atm_t *atm,
                             int fourier_order,
                             rt_surface_kind_t surface_kind,
                             double n_water,
                             double wind_speed,
                             int sigma_type,
                             int q_convention);

void rt_surface_boundary_free(rt_surface_boundary_t *boundary);

/* Apply the lower-boundary reflection operator to the downward Stokes field at
 * the final atmospheric level. */
int rt_surface_boundary_apply(const rt_surface_boundary_t *boundary,
                              const rt_atm_t *atm,
                              double n_water,
                              int q_convention,
                              const double *I_field,
                              const double *Q_field,
                              const double *U_field,
                              double *bc_I,
                              double *bc_Q,
                              double *bc_U);

/* Build all first-order lower-boundary sources used to seed atmosphere-surface
 * multiple scattering.  The caller supplies zero-initialized order1 and
 * direct_seed arrays.  direct_seed contains only the rough-surface direct
 * component that may be removed from a glint-decoupled final product. */
int rt_surface_build_initial_sources(const rt_surface_boundary_t *boundary,
                                     const rt_atm_t *atm,
                                     const rt_legendre_workspace_t *ws,
                                     int fourier_order,
                                     const rt_solver_sos_options_t *opts,
                                     rt_surface_kind_t surface_kind,
                                     double n_water,
                                     double wind_speed,
                                     int sigma_type,
                                     int q_convention,
                                     double *total_I,
                                     double *total_Q,
                                     double *total_U,
                                     double *order1_I,
                                     double *order1_Q,
                                     double *order1_U,
                                     double *direct_seed_I,
                                     double *direct_seed_Q,
                                     double *direct_seed_U);

void rt_surface_remove_direct_seed(const rt_atm_t *atm,
                                   const double *direct_seed_I,
                                   const double *direct_seed_Q,
                                   const double *direct_seed_U,
                                   double *total_I,
                                   double *total_Q,
                                   double *total_U);

#endif /* OCRT_RT_SURFACE_BOUNDARY_H */
