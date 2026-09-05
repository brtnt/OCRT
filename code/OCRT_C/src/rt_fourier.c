#include "rt_fourier.h"

#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double rt_fourier_reconstruct_cos(const double *mode_values,
                                  int m_max,
                                  double public_relative_azimuth_rad)
{
    if (!mode_values || m_max < 0) return 0.0;

    const double phi = public_relative_azimuth_rad + M_PI;
    double value = mode_values[0];
    for (int m = 1; m <= m_max; ++m)
        value += 2.0 * mode_values[m] * cos((double)m * phi);
    return value;
}

double rt_fourier_reconstruct_sin(const double *mode_values,
                                  int m_max,
                                  double public_relative_azimuth_rad)
{
    if (!mode_values || m_max < 1) return 0.0;

    const double phi = public_relative_azimuth_rad + M_PI;
    double value = 0.0;
    for (int m = 1; m <= m_max; ++m)
        value += 2.0 * mode_values[m] * sin((double)m * phi);
    return value;
}

double rt_fourier_pi_shift_sign(int m)
{
    return (m & 1) ? -1.0 : 1.0;
}
