#ifndef OCRT_RT_SOS_OPERATOR_H
#define OCRT_RT_SOS_OPERATOR_H

#include "rt_kernel.h"
#include "rt_types.h"

/* Prepare the immutable direction-pair operator for one azimuthal Fourier
 * order.  The packed operator and reusable scratch arena are owned by ws and
 * remain valid until any angular kernel in ws is rebuilt. */
int rt_sos_operator_prepare(const rt_atm_t *atm,
                            int fourier_order,
                            rt_legendre_workspace_t *ws);

/* Apply the vector scattering operator to one SOS radiance order.  The solar
 * slot (signed direction index zero) is left as NaN; all physical ordinates
 * are written. */
int rt_sos_operator_apply_vector(const rt_atm_t *atm,
                                 int fourier_order,
                                 const rt_legendre_workspace_t *ws,
                                 const double *I_previous,
                                 const double *Q_previous,
                                 const double *U_previous,
                                 double *J_I,
                                 double *J_Q,
                                 double *J_U);

/* Scalar counterpart retained for the low-level diagnostic scalar API. */
int rt_sos_operator_apply_scalar(const rt_atm_t *atm,
                                 int fourier_order,
                                 const rt_legendre_workspace_t *ws,
                                 const double *I_previous,
                                 double *J_I);

#endif /* OCRT_RT_SOS_OPERATOR_H */
