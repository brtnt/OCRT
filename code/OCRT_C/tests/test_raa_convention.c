#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "rt_raa_convention.h"
#include "rt_fourier.h"

static void require_close(const char *name, double got, double expected, double tol)
{
    if (!isfinite(got) || fabs(got - expected) > tol) {
        fprintf(stderr, "FAIL %s: got %.17g expected %.17g tol %.3g\n",
                name, got, expected, tol);
        exit(1);
    }
}

int main(void)
{
    const double pi = RT_F_SOLAR_PI;
    const double mu = cos(pi / 6.0); /* 30 deg */

    require_close("norm -180", rt_norm_raa_0_360(-180.0), 180.0, 1e-14);
    require_close("norm 360", rt_norm_raa_0_360(360.0), 0.0, 1e-14);

    /* Atmospheric public branches, SZA=VZA=30 deg. */
    require_close("atm cosTheta RAA0",
                  rt_atm_scatter_cos_from_public_raa(mu, mu, 0.0),
                  -0.5, 1e-14);
    require_close("atm cosTheta RAA180",
                  rt_atm_scatter_cos_from_public_raa(mu, mu, 180.0),
                  -1.0, 1e-14);

    /* Local water single-scatter geometry deliberately uses pi-RAA. */
    require_close("water cosTheta RAA0",
                  rt_water_scatter_cos_from_public_raa(mu, mu, 0.0),
                  -1.0, 1e-14);
    require_close("water cosTheta RAA180",
                  rt_water_scatter_cos_from_public_raa(mu, mu, 180.0),
                  -0.5, 1e-14);
    require_close("water phi RAA0", rt_raa_to_water_scatter_phi(0.0), pi, 1e-14);
    require_close("water phi RAA180", rt_raa_to_water_scatter_phi(180.0), 0.0, 1e-14);


    /* Reported atmospheric and water Fourier outputs share the public RAA.
     * The synthetic modes make the 0/180 branch swap and U sign observable. */
    {
        const double iq_modes[3] = {1.0, 0.25, -0.10};
        const double u_modes[3]  = {0.0, 0.20, -0.05};
        const double p0   = rt_raa_to_atm_fourier_phi(0.0);
        const double p45  = rt_raa_to_atm_fourier_phi(45.0);
        const double p180 = rt_raa_to_atm_fourier_phi(180.0);

        require_close("public phi RAA0", p0, 0.0, 1e-14);
        require_close("public phi RAA180", p180, pi, 1e-14);
        require_close("IQ reconstruct RAA0",
                      rt_fourier_reconstruct_cos(iq_modes, 2, p0),
                      1.0 + 2.0*0.25*cos(pi) - 2.0*0.10*cos(2.0*pi),
                      1e-14);
        require_close("IQ reconstruct RAA180",
                      rt_fourier_reconstruct_cos(iq_modes, 2, p180),
                      1.0 + 2.0*0.25*cos(2.0*pi) - 2.0*0.10*cos(4.0*pi),
                      1e-14);
        require_close("U reconstruct RAA45",
                      rt_fourier_reconstruct_sin(u_modes, 2, p45),
                      2.0*0.20*sin(5.0*pi/4.0)
                        - 2.0*0.05*sin(5.0*pi/2.0),
                      1e-14);
    }

    /* Surface local azimuth must put the direct-glint branch at RAA=180. */
    require_close("surface phi RAA0", rt_raa_to_surface_view_phi(0.0), -pi, 1e-14);
    require_close("surface phi RAA180", rt_raa_to_surface_view_phi(180.0), 0.0, 1e-14);

    puts("PASS: OCRT public RAA convention helpers");
    return 0;
}
