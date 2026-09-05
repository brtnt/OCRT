/* Gate E1: 입력 검증과 오류코드. */
#include <stdio.h>
#include <math.h>
#include "../src/rt_eap_mie_phase.h"
static int fails = 0;
static void ck(const char *what, int got, int want)
{
    if (got != want) { printf("  FAIL %-30s got %d want %d\n", what, got, want); ++fails; }
    else               printf("  ok   %-30s (%d)\n", what, got);
}
int main(void)
{
    enum { NT = 5, NW = 2 };
    double th[NT] = {0.0, 45.0, 90.0, 135.0, 180.0};
    double wl[NW] = {443.0, 675.0};
    double a[NW * NT], b[NW * NT], c[NW * NT];
    ck("null wavelength",  rt_eap_mie_phase_compute(0, NULL, NW, th, NT, a, b, c), RT_EAP_PHASE_EINVAL);
    ck("null output",      rt_eap_mie_phase_compute(0, wl, NW, th, NT, NULL, b, c), RT_EAP_PHASE_EINVAL);
    ck("n_wavelength = 0", rt_eap_mie_phase_compute(0, wl, 0, th, NT, a, b, c), RT_EAP_PHASE_EINVAL);
    ck("n_theta < 3",      rt_eap_mie_phase_compute(0, wl, NW, th, 2, a, b, c), RT_EAP_PHASE_EINVAL);
    ck("species >= 17",    rt_eap_mie_phase_compute(17, wl, NW, th, NT, a, b, c), RT_EAP_PHASE_ESPECIES);
    double bad[NW] = {443.0, NAN};
    ck("NaN wavelength",   rt_eap_mie_phase_compute(0, bad, NW, th, NT, a, b, c), RT_EAP_PHASE_ERANGE);
    double oor[1] = {200.0};
    ck("wavelength < range", rt_eap_mie_phase_compute(0, oor, 1, th, NT, a, b, c), RT_EAP_PHASE_ERANGE);
    double nm[NT] = {0.0, 90.0, 45.0, 135.0, 180.0};
    ck("non-monotonic angles", rt_eap_mie_phase_compute(0, wl, NW, nm, NT, a, b, c), RT_EAP_PHASE_EINVAL);
    double ne[NT] = {5.0, 45.0, 90.0, 135.0, 180.0};
    ck("no 0 deg endpoint", rt_eap_mie_phase_compute(0, wl, NW, ne, NT, a, b, c), RT_EAP_PHASE_ERANGE);
    double ne2[NT] = {0.0, 45.0, 90.0, 135.0, 170.0};
    ck("no 180 deg endpoint", rt_eap_mie_phase_compute(0, wl, NW, ne2, NT, a, b, c), RT_EAP_PHASE_ERANGE);
    ck("valid request",    rt_eap_mie_phase_compute(0, wl, NW, th, NT, a, b, c), RT_EAP_PHASE_OK);
    printf(fails ? "\nFAILED (%d)\n" : "\nall passed\n", fails);
    return fails ? 1 : 0;
}
