/* rt_eap_coated_mie.h -- internal only.  Aden-Kerker stratified sphere. */
#include <complex.h>
#ifndef OCRT_RT_EAP_COATED_MIE_H
#define OCRT_RT_EAP_COATED_MIE_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
enum {
    RT_EAP_MIE_OK       =  0,
    RT_EAP_MIE_EINVAL   = -1,
    RT_EAP_MIE_ERANGE   = -2,   /* size parameter beyond the supported range */
    RT_EAP_MIE_ENOMEM   = -3,
    RT_EAP_MIE_ENUMERIC = -4,
    RT_EAP_MIE_ENOCONV  = -5
};
typedef struct {
    double qext, qsca, gqsc, qbs;
    int    n_terms;             /* series order actually used */
} rt_eap_mie_eff_t;

/* With numang > 0, mu[numang] must lie in [0, 1] and the routine returns,
 * for each mu, both the forward-hemisphere angle (k = 0) and its supplement
 * (k = 1). Output layout is out[j * 2 + k], matching the reference routine.
 * With numang == 0, mu and all four angular-output pointers may be NULL; only
 * the efficiency structure is computed. This is used by the .mie writer's
 * spectral table and does not change the phase-calculation path. */
int rt_eap_coated_mie(double rcore, double rshell, double wvno,
                      double _Complex rindsh, double _Complex rindco,
                      const double *mu, size_t numang,
                      double *m1, double *m2, double *s21, double *d21,
                      rt_eap_mie_eff_t *eff);
#ifdef __cplusplus
}
#endif
#endif
