#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "shared/mie_io.h"

static double fr631_desc(int i)
{
    if (i <= 320) return 180.0 - 0.5 * i;
    if (i <= 470) return 19.9 - 0.1 * (i - 321);
    if (i <= 550) return 4.95 - 0.05 * (i - 471);
    if (i <= 590) return 0.98 - 0.02 * (i - 551);
    return 0.195 - 0.005 * (i - 591);
}

int main(void)
{
    double *fr = calloc(631, sizeof(double));
    double *legacy = calloc(361, sizeof(double));
    double custom[5] = {180.0, 20.0, 1.0, 0.1, 0.0};
    if (!fr || !legacy) return 2;
    for (int i = 0; i < 631; ++i) fr[i] = fr631_desc(i);
    for (int i = 0; i < 361; ++i) legacy[i] = 180.0 - 0.5 * i;

    if (mie_classify_angle_grid(fr, 631) != MIE_ANGLE_GRID_FR631) return 10;
    if (mie_classify_angle_grid(legacy, 361) != MIE_ANGLE_GRID_LEGACY_361_UNIFORM) return 11;
    if (mie_classify_angle_grid(custom, 5) != MIE_ANGLE_GRID_NONCANONICAL) return 12;

    for (int i = 0; i < 631 / 2; ++i) {
        const double t = fr[i]; fr[i] = fr[630 - i]; fr[630 - i] = t;
    }
    if (mie_classify_angle_grid(fr, 631) != MIE_ANGLE_GRID_FR631) return 13;

    free(fr); free(legacy);
    puts("PASS mie angle-grid policy");
    return 0;
}
