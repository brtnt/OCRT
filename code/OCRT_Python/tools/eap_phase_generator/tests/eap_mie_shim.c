/* ctypes 시험용 껍데기: 복소수 인자를 실수 쌍으로 받는다. */
#include <complex.h>
#include "../src/internal/rt_eap_coated_mie.h"
int rt_eap_coated_mie_ri(double rcore, double rshell, double wvno,
                         double nsh, double ksh, double nco, double kco,
                         const double *mu, size_t numang,
                         double *m1, double *m2, double *s21, double *d21,
                         rt_eap_mie_eff_t *eff)
{
    return rt_eap_coated_mie(rcore, rshell, wvno,
                             nsh + ksh * _Complex_I, nco + kco * _Complex_I,
                             mu, numang, m1, m2, s21, d21, eff);
}
