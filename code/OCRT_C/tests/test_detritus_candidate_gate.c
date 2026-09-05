#include <math.h>
#include <stdio.h>
#include "rt_iop_organic.h"

int main(int argc, char **argv)
{
    if (argc != 2) return 2;
    if (rt_iop_organic_init_selected(argv[1], ORGANIC_PHYTO_MICRO, 0) != 0)
        return 3;
    if (fabs(rt_iop_organic_detritus_wavelength_min_nm() - 350.0) > 1e-12 ||
        fabs(rt_iop_organic_detritus_wavelength_max_nm() - 850.0) > 1e-12)
        return 4;
    if (rt_iop_organic_detritus_wavelength_supported(330.0) ||
        !rt_iop_organic_detritus_wavelength_supported(350.0) ||
        !rt_iop_organic_detritus_wavelength_supported(850.0) ||
        rt_iop_organic_detritus_wavelength_supported(1100.0))
        return 5;
    rt_iop_t x = {0};
    if (rt_iop_detritus_eval(330.0, 1.0, 0.0,
                             ORGANIC_DETRITUS_SLOPE_DEFAULT, &x) != -7)
        return 6;
    if (rt_iop_detritus_eval(443.0, 1.0, 0.0,
                             ORGANIC_DETRITUS_SLOPE_DEFAULT, &x) != 0 ||
        !(x.b > 0.0 && x.bb > 0.0))
        return 7;
    if (rt_iop_detritus_eval(1100.0, 1.0, 0.0,
                             ORGANIC_DETRITUS_SLOPE_DEFAULT, &x) != -7)
        return 8;
    rt_iop_organic_free();
    puts("PASS: validated detritus active 350-850; rejected 330-1100 candidate is not runtime-active");
    return 0;
}
