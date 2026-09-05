/*
 * OCRT atmosphere/surface boundary operator.
 *
 * The implementation follows the Fourier convolution of a Mueller reflection
 * operator with the downward Stokes field.  Flat Fresnel reflection is a
 * specular delta represented by a diagonal discrete operator; a rough ocean
 * uses the azimuthal Fourier coefficients of the Cox-Munk microfacet BRDF.
 * Direct lower-boundary sources are prepared once before the SOS order loop.
 */

#include "rt_surface_boundary.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "shared/surface.h"
#include "rt_angular_basis.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double moment_gr_pair(const rt_legendre_workspace_t *ws,
                             const double *gamma,
                             int m,
                             int first,
                             int second)
{
    if (!ws || !gamma || m < 0 || m > ws->l_max) return 0.0;
    double value = 0.0;
    for (int l = m; l <= ws->l_max; ++l) {
        value += ws->plm[l][first] * ws->rrl[l][second] * gamma[l];
    }
    return value;
}

static double moment_art_pair(const rt_legendre_workspace_t *ws,
                              const double *alpha,
                              const double *zeta,
                              int m,
                              int first,
                              int second)
{
    if (!ws || !alpha || !zeta || m < 0 || m > ws->l_max) return 0.0;
    double value = 0.0;
    for (int l = m; l <= ws->l_max; ++l) {
        value += ws->rtl[l][first] * ws->rrl[l][second] * alpha[l] +
                 ws->rrl[l][first] * ws->rtl[l][second] * zeta[l];
    }
    return value;
}


typedef struct {
    double P;
    double R;
    double T;
    double RT;
    double RR;
    double TR;
} rt_reflected_phase_kernel_t;

int rt_surface_boundary_init(rt_surface_boundary_t *boundary,
                             const rt_atm_t *atm,
                             int m,
                             rt_surface_kind_t surface_kind,
                             double n_water,
                             double wind_speed,
                             int sigma_type,
                             int q_convention)
{
    if (!boundary || !atm || m < 0 || atm->n_mu < 1 || !(n_water > 0.0)) {
        return -1;
    }
    memset(boundary, 0, sizeof(*boundary));

    const int n_mu = atm->n_mu;
    boundary->n_mu = n_mu;
    boundary->fourier_order = m;
    boundary->is_flat = (surface_kind == RT_SURFACE_FLAT) ||
                        (surface_kind == RT_SURFACE_BLACK_FRESNEL_OCEAN && wind_speed <= 0.0);
    boundary->mu_positive = (double *)malloc((size_t)n_mu * sizeof(double));
    boundary->kernel_3x3 =
        (double *)calloc((size_t)n_mu * (size_t)n_mu * 9u, sizeof(double));
    if (!boundary->mu_positive || !boundary->kernel_3x3) {
        rt_surface_boundary_free(boundary);
        return -1;
    }
    for (int j = 1; j <= n_mu; ++j) boundary->mu_positive[j - 1] = atm->rm[j];

    if (boundary->is_flat) {
        for (int outgoing = 0; outgoing < n_mu; ++outgoing) {
            const double weight_mu = atm->gb[outgoing + 1];
            if (weight_mu <= 0.0) continue;
            double fresnel[9];
            surface_flat_fresnel_matrix_raw(boundary->mu_positive[outgoing],
                                            n_water, q_convention, fresnel);
            const double discrete_delta =
                1.0 / (2.0 * M_PI * boundary->mu_positive[outgoing] * weight_mu);
            double *entry = boundary->kernel_3x3 +
                ((size_t)outgoing * (size_t)n_mu + (size_t)outgoing) * 9u;
            for (int element = 0; element < 9; ++element) {
                entry[element] = discrete_delta * fresnel[element];
            }
        }
        return 0;
    }

    int n_phi = 1024;
    const char *override = getenv("OCRT_COXMUNK_FOURIER_NPHI");
    if (override && override[0]) {
        const int requested = atoi(override);
        if (requested >= 64 && requested <= 4096) n_phi = requested;
    }
    const int status = surface_coxmunk_fourier_kernel(
        boundary->mu_positive, n_mu,
        boundary->mu_positive, n_mu,
        m, n_phi, wind_speed, sigma_type, n_water, q_convention,
        boundary->kernel_3x3);
    if (status != 0) {
        rt_surface_boundary_free(boundary);
        return status;
    }
    return 0;
}

void rt_surface_boundary_free(rt_surface_boundary_t *boundary)
{
    if (!boundary) return;
    free(boundary->mu_positive);
    free(boundary->kernel_3x3);
    memset(boundary, 0, sizeof(*boundary));
}

int rt_surface_boundary_apply(const rt_surface_boundary_t *boundary,
                              const rt_atm_t *atm,
                              double n_water,
                              int q_convention,
                              const double *I_field,
                              const double *Q_field,
                              const double *U_field,
                              double *bc_I,
                              double *bc_Q,
                              double *bc_U)
{
    if (!boundary || !atm || !I_field || !Q_field || !U_field ||
        !bc_I || !bc_Q || !bc_U || boundary->n_mu != atm->n_mu) {
        return -1;
    }
    const int n_mu = atm->n_mu;
    const int directions = 2 * n_mu + 1;
    const int bottom = atm->n_layers;
    const double Fourier_weight = boundary->is_flat
        ? (2.0 * RT_F_SOLAR_PI)
        : ((boundary->fourier_order == 0)
            ? (2.0 * RT_F_SOLAR_PI) : RT_F_SOLAR_PI);

    for (int outgoing = 1; outgoing <= n_mu; ++outgoing) {
        double I_sum = 0.0;
        double Q_sum = 0.0;
        double U_sum = 0.0;

        if (boundary->is_flat && atm->gb[outgoing] == 0.0) {
            double fresnel[9];
            surface_flat_fresnel_matrix_raw(
                boundary->mu_positive[outgoing - 1], n_water,
                q_convention, fresnel);
            const size_t down = (size_t)bottom * (size_t)directions +
                                (size_t)(n_mu - outgoing);
            const double incident_I = I_field[down];
            const double incident_Q = Q_field[down];
            const double incident_U = U_field[down];
            I_sum = fresnel[0] * incident_I + fresnel[1] * incident_Q +
                    fresnel[2] * incident_U;
            Q_sum = fresnel[3] * incident_I + fresnel[4] * incident_Q +
                    fresnel[5] * incident_U;
            U_sum = fresnel[6] * incident_I + fresnel[7] * incident_Q +
                    fresnel[8] * incident_U;
        } else {
            for (int incoming = 1; incoming <= n_mu; ++incoming) {
                const double quadrature_weight = atm->gb[incoming];
                if (quadrature_weight == 0.0) continue;
                const size_t down = (size_t)bottom * (size_t)directions +
                                    (size_t)(n_mu - incoming);
                const double incident_I = I_field[down];
                const double incident_Q = Q_field[down];
                const double incident_U = U_field[down];
                const double projected_weight = Fourier_weight *
                    boundary->mu_positive[incoming - 1] * quadrature_weight;
                const double *R = boundary->kernel_3x3 +
                    ((size_t)(outgoing - 1) * (size_t)n_mu +
                     (size_t)(incoming - 1)) * 9u;
                I_sum += projected_weight *
                    (R[0] * incident_I + R[1] * incident_Q + R[2] * incident_U);
                Q_sum += projected_weight *
                    (R[3] * incident_I + R[4] * incident_Q + R[5] * incident_U);
                U_sum += projected_weight *
                    (R[6] * incident_I + R[7] * incident_Q + R[8] * incident_U);
            }
        }
        bc_I[outgoing - 1] = I_sum;
        bc_Q[outgoing - 1] = Q_sum;
        bc_U[outgoing - 1] = U_sum;
    }
    return 0;
}

static int add_flat_reflected_beam_scattering(
    const rt_surface_boundary_t *boundary,
    const rt_atm_t *atm,
    const rt_legendre_workspace_t *ws,
    int m,
    double n_water,
    int q_convention,
    double *total_I, double *total_Q, double *total_U,
    double *order1_I, double *order1_Q, double *order1_U)
{
    const int n_layers = atm->n_layers;
    const int n_mu = atm->n_mu;
    const int directions = 2 * n_mu + 1;
    const size_t field_count =
        (size_t)(n_layers + 1) * (size_t)directions;
    const double mu_sun = atm->mu_sun;

    /* Fresnel reflection of the downward direct beam.  This equivalent
     * algebra deliberately uses the same amplitude ordering as the historical
     * OCRT public outputs. */
    const double n2 = n_water * n_water;
    const double transmitted_numerator = sqrt(n2 - 1.0 + mu_sun * mu_sun);
    const double rp = (n2 * mu_sun - transmitted_numerator) /
                      (n2 * mu_sun + transmitted_numerator);
    const double rs = (mu_sun - transmitted_numerator) /
                      (mu_sun + transmitted_numerator);
    const double reflected_I_pp = (rp * rp + rs * rs) / 2.0;
    const double reflected_Q_pp = (q_convention == 1)
        ? (rp * rp - rs * rs) / 2.0
        : (rs * rs - rp * rp) / 2.0;

    const double tau_total = atm->h[n_layers];
    const double surface_beam = 0.5 * exp(-2.0 * tau_total / mu_sun);
    const double beam_P = atm->xpl[0];
    const double beam_R = atm->xrl[0];
    const double beta0 = (m == 0) ? atm->beta0 : 0.0;
    const double beta2 = atm->beta2;
    const double gamma2 = atm->gamma2;
    const double alpha2 = atm->alpha2;
    const int molecular_active = (m <= 2);
    const double *gamma_particle =
        (atm->aerosol_active && atm->gammal_aer) ? atm->gammal_aer : NULL;
    const double *alpha_particle =
        (atm->aerosol_active && atm->alphal_aer) ? atm->alphal_aer : NULL;
    const double *zeta_particle =
        (atm->aerosol_active && atm->zetal_aer) ? atm->zetal_aer : NULL;

    double *source_I = (double *)calloc(field_count, sizeof(double));
    double *source_Q = (double *)calloc(field_count, sizeof(double));
    double *source_U = (double *)calloc(field_count, sizeof(double));
    double *field_I = (double *)calloc(field_count, sizeof(double));
    double *field_Q = (double *)calloc(field_count, sizeof(double));
    double *field_U = (double *)calloc(field_count, sizeof(double));
    if (!source_I || !source_Q || !source_U ||
        !field_I || !field_Q || !field_U) {
        free(source_I); free(source_Q); free(source_U);
        free(field_I); free(field_Q); free(field_U);
        return -1;
    }


    for (int layer = 0; layer < n_layers; ++layer) {
        const double molecular_top = atm->ydel[layer];
        const double particle_top = atm->xdel[layer];
        const double molecular_bottom = atm->ydel[layer + 1];
        const double particle_bottom = atm->xdel[layer + 1];
        double beam_top = surface_beam * exp(atm->h[layer] / mu_sun);
        double beam_bottom = surface_beam * exp(+atm->h[layer + 1] / mu_sun);
        double reflected_I_top = reflected_I_pp, reflected_Q_top = reflected_Q_pp;
        double reflected_I_bottom = reflected_I_pp, reflected_Q_bottom = reflected_Q_pp;

        for (int direction = 1; direction <= n_mu; ++direction) {
            rt_reflected_phase_kernel_t top_kernel;
            rt_reflected_phase_kernel_t bottom_kernel;

            {
                const double p_up = atm->xpl[direction];
                const double p_down = atm->xpl[-direction];
                const double r_up = atm->xrl[direction];
                const double r_down = atm->xrl[-direction];
                const double t_up = atm->xtl[direction];
                const double t_down = atm->xtl[-direction];

                const double particle_P_down = ws->phase_fourier_m[0][-direction];
                const double particle_R_down = ws->gr_pol[0][-direction];
                const double particle_T_down = ws->gt_pol[0][-direction];
                const double particle_R_transpose_down =
                    moment_gr_pair(ws, gamma_particle, m, -direction, 0);
                const double particle_RR_down = ws->arr_pol[0][-direction];
                const double particle_TR_down =
                    moment_art_pair(ws, alpha_particle, zeta_particle,
                                    m, -direction, 0);

                const double particle_P_up = ws->phase_fourier_m[0][direction];
                const double particle_R_up = ws->gr_pol[0][direction];
                const double particle_T_up = ws->gt_pol[0][direction];
                const double particle_R_transpose_up =
                    moment_gr_pair(ws, gamma_particle, m, direction, 0);
                const double particle_RR_up = ws->arr_pol[0][direction];
                const double particle_TR_up =
                    moment_art_pair(ws, alpha_particle, zeta_particle,
                                    m, direction, 0);

                const double molecular_P_down = molecular_active
                    ? beta0 + beta2 * p_down * beam_P : 0.0;
                const double molecular_R_down = molecular_active
                    ? gamma2 * r_down * beam_P : 0.0;
                const double molecular_T_down = molecular_active
                    ? gamma2 * t_down * beam_P : 0.0;
                const double molecular_R_transpose_down = molecular_active
                    ? gamma2 * beam_R * p_down : 0.0;
                const double molecular_RR_down = molecular_active
                    ? alpha2 * beam_R * r_down : 0.0;
                const double molecular_TR_down = molecular_active
                    ? alpha2 * t_down * beam_R : 0.0;

                const double molecular_P_up = molecular_active
                    ? beta0 + beta2 * p_up * beam_P : 0.0;
                const double molecular_R_up = molecular_active
                    ? gamma2 * r_up * beam_P : 0.0;
                const double molecular_T_up = molecular_active
                    ? gamma2 * t_up * beam_P : 0.0;
                const double molecular_R_transpose_up = molecular_active
                    ? gamma2 * beam_R * p_up : 0.0;
                const double molecular_RR_up = molecular_active
                    ? alpha2 * beam_R * r_up : 0.0;
                const double molecular_TR_up = molecular_active
                    ? alpha2 * t_up * beam_R : 0.0;

                top_kernel.P = particle_P_down * particle_top +
                               molecular_P_down * molecular_top;
                top_kernel.R = particle_R_down * particle_top +
                               molecular_R_down * molecular_top;
                top_kernel.T = particle_T_down * particle_top +
                               molecular_T_down * molecular_top;
                top_kernel.RT = particle_R_transpose_down * particle_top +
                                molecular_R_transpose_down * molecular_top;
                top_kernel.RR = particle_RR_down * particle_top +
                                molecular_RR_down * molecular_top;
                top_kernel.TR = particle_TR_down * particle_top +
                                molecular_TR_down * molecular_top;

                bottom_kernel.P = particle_P_up * particle_bottom +
                                  molecular_P_up * molecular_bottom;
                bottom_kernel.R = particle_R_up * particle_bottom +
                                  molecular_R_up * molecular_bottom;
                bottom_kernel.T = particle_T_up * particle_bottom +
                                  molecular_T_up * molecular_bottom;
                bottom_kernel.RT =
                    particle_R_transpose_up * particle_bottom +
                    molecular_R_transpose_up * molecular_bottom;
                bottom_kernel.RR = particle_RR_up * particle_bottom +
                                   molecular_RR_up * molecular_bottom;
                bottom_kernel.TR = particle_TR_up * particle_bottom +
                                   molecular_TR_up * molecular_bottom;
            }

            const size_t up = (size_t)layer * (size_t)directions +
                              (size_t)(n_mu + direction);
            source_I[up] = beam_top *
                (reflected_I_top * top_kernel.P + reflected_Q_top * top_kernel.RT);
            source_Q[up] = beam_top *
                (reflected_I_top * top_kernel.R + reflected_Q_top * top_kernel.RR);
            source_U[up] = beam_top *
                (reflected_I_top * top_kernel.T + reflected_Q_top * top_kernel.TR);

            const size_t down = (size_t)(layer + 1) * (size_t)directions +
                                (size_t)(n_mu - direction);
            source_I[down] = beam_bottom *
                (reflected_I_bottom * bottom_kernel.P + reflected_Q_bottom * bottom_kernel.RT);
            source_Q[down] = beam_bottom *
                (reflected_I_bottom * bottom_kernel.R + reflected_Q_bottom * bottom_kernel.RR);
            source_U[down] = beam_bottom *
                (reflected_I_bottom * bottom_kernel.T + reflected_Q_bottom * bottom_kernel.TR);
        }
    }
    for (int level = 0; level <= n_layers; ++level) {
        const size_t solar = (size_t)level * (size_t)directions + (size_t)n_mu;
        source_I[solar] = NAN;
        source_Q[solar] = NAN;
        source_U[solar] = NAN;
    }

    int status = 0;
    if (rt_solver_integrate(atm, m, source_I,
                            RT_INTEGRATION_METHOD_LINEAR, field_I) != 0 ||
        rt_solver_integrate(atm, m, source_Q,
                            RT_INTEGRATION_METHOD_LINEAR, field_Q) != 0 ||
        rt_solver_integrate(atm, m, source_U,
                            RT_INTEGRATION_METHOD_LINEAR, field_U) != 0) {
        status = -2;
    } else {
        for (int level = 0; level <= n_layers; ++level) {
            for (int direction = 1; direction <= n_mu; ++direction) {
                const size_t up = (size_t)level * (size_t)directions +
                                  (size_t)(n_mu + direction);
                const size_t down = (size_t)level * (size_t)directions +
                                    (size_t)(n_mu - direction);
                total_I[up] += field_I[up];
                total_Q[up] += field_Q[up];
                total_U[up] += field_U[up];
                total_I[down] += field_I[down];
                total_Q[down] += field_Q[down];
                total_U[down] += field_U[down];
                order1_I[up] += field_I[up];
                order1_Q[up] += field_Q[up];
                order1_U[up] += field_U[up];
                order1_I[down] += field_I[down];
                order1_Q[down] += field_Q[down];
                order1_U[down] += field_U[down];
            }
        }
    }

    free(source_I); free(source_Q); free(source_U);
    free(field_I); free(field_Q); free(field_U);
    (void)boundary;
    return status;
}

static int add_rough_direct_beam(
    const rt_surface_boundary_t *boundary,
    const rt_atm_t *atm,
    int m,
    double n_water,
    double wind_speed,
    int sigma_type,
    int q_convention,
    double *total_I, double *total_Q, double *total_U,
    double *order1_I, double *order1_Q, double *order1_U,
    double *direct_I, double *direct_Q, double *direct_U)
{
    const int n_mu = atm->n_mu;
    const int n_layers = atm->n_layers;
    const int directions = 2 * n_mu + 1;
    const double mu_sun = atm->mu_sun;
    const double tau_total = atm->h[n_layers];
    const double down_transmission = exp(-tau_total / mu_sun);
    const double Fourier_scale =
        ((m == 0) ? RT_F_SOLAR_PI : 0.5 * RT_F_SOLAR_PI) * mu_sun;

    int n_phi = 1024;
    const char *override = getenv("OCRT_COXMUNK_FOURIER_NPHI");
    if (override && override[0]) {
        const int requested = atoi(override);
        if (requested >= 64 && requested <= 4096) n_phi = requested;
    }

    double solar_mu[1] = {mu_sun};
    double *solar_kernel = (double *)calloc((size_t)n_mu * 9u, sizeof(double));
    if (!solar_kernel) return -1;
    const int status = surface_coxmunk_fourier_kernel(
        boundary->mu_positive, n_mu, solar_mu, 1, m, n_phi,
        wind_speed, sigma_type, n_water, q_convention, solar_kernel);
    if (status != 0) {
        free(solar_kernel);
        return status;
    }

    for (int level = 0; level <= n_layers; ++level) {
        const double upward_depth = tau_total - atm->h[level];
        for (int outgoing = 0; outgoing < n_mu; ++outgoing) {
            const double transmission = down_transmission *
                exp(-upward_depth / boundary->mu_positive[outgoing]);
            const double factor = Fourier_scale * transmission;
            const double I_value = solar_kernel[(size_t)outgoing * 9u + 0] * factor;
            const double Q_value = solar_kernel[(size_t)outgoing * 9u + 3] * factor;
            const double U_value = solar_kernel[(size_t)outgoing * 9u + 6] * factor;
            const size_t up = (size_t)level * (size_t)directions +
                              (size_t)(n_mu + 1 + outgoing);
            total_I[up] += I_value;
            total_Q[up] += Q_value;
            total_U[up] += U_value;
            order1_I[up] += I_value;
            order1_Q[up] += Q_value;
            order1_U[up] += U_value;
            direct_I[up] = I_value;
            direct_Q[up] = Q_value;
            direct_U[up] = U_value;
        }
    }
    free(solar_kernel);
    return 0;
}

static void add_external_bottom_source(
    const rt_surface_boundary_t *boundary,
    const rt_atm_t *atm,
    const rt_solver_sos_options_t *opts,
    double *total_I, double *total_Q, double *total_U,
    double *order1_I, double *order1_Q, double *order1_U)
{
    if (!opts->ext_bottom_I || !opts->ext_bottom_Q || !opts->ext_bottom_U) {
        return;
    }
    const int n_mu = atm->n_mu;
    const int n_layers = atm->n_layers;
    const int directions = 2 * n_mu + 1;
    const double tau_total = atm->h[n_layers];
    for (int level = 0; level <= n_layers; ++level) {
        const double upward_depth = tau_total - atm->h[level];
        for (int outgoing = 0; outgoing < n_mu; ++outgoing) {
            const double attenuation =
                exp(-upward_depth / boundary->mu_positive[outgoing]);
            const double I_value = opts->ext_bottom_I[outgoing] * attenuation;
            const double Q_value = opts->ext_bottom_Q[outgoing] * attenuation;
            const double U_value = opts->ext_bottom_U[outgoing] * attenuation;
            const size_t up = (size_t)level * (size_t)directions +
                              (size_t)(n_mu + 1 + outgoing);
            total_I[up] += I_value;
            total_Q[up] += Q_value;
            total_U[up] += U_value;
            order1_I[up] += I_value;
            order1_Q[up] += Q_value;
            order1_U[up] += U_value;
        }
    }
}

int rt_surface_build_initial_sources(const rt_surface_boundary_t *boundary,
                                     const rt_atm_t *atm,
                                     const rt_legendre_workspace_t *ws,
                                     int m,
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
                                     double *direct_seed_U)
{
    if (!boundary || !atm || !ws || !opts ||
        !total_I || !total_Q || !total_U ||
        !order1_I || !order1_Q || !order1_U) return -1;

    if (!opts->bottom_source_only && boundary->is_flat) {
        const int status = add_flat_reflected_beam_scattering(
            boundary, atm, ws, m, n_water, q_convention,
            total_I, total_Q, total_U, order1_I, order1_Q, order1_U);
        if (status != 0) return status;
    }

    if (!opts->bottom_source_only &&
        surface_kind == RT_SURFACE_BLACK_FRESNEL_OCEAN && wind_speed > 0.0) {
        if (!direct_seed_I || !direct_seed_Q || !direct_seed_U) return -1;
        const int status = add_rough_direct_beam(
            boundary, atm, m, n_water, wind_speed, sigma_type, q_convention,
            total_I, total_Q, total_U,
            order1_I, order1_Q, order1_U,
            direct_seed_I, direct_seed_Q, direct_seed_U);
        if (status != 0) return status;
    }

    add_external_bottom_source(boundary, atm, opts,
                               total_I, total_Q, total_U,
                               order1_I, order1_Q, order1_U);
    return 0;
}

void rt_surface_remove_direct_seed(const rt_atm_t *atm,
                                   const double *direct_seed_I,
                                   const double *direct_seed_Q,
                                   const double *direct_seed_U,
                                   double *total_I,
                                   double *total_Q,
                                   double *total_U)
{
    if (!atm || !direct_seed_I || !direct_seed_Q || !direct_seed_U ||
        !total_I || !total_Q || !total_U) return;
    const int n_mu = atm->n_mu;
    const int directions = 2 * n_mu + 1;
    for (int level = 0; level <= atm->n_layers; ++level) {
        for (int outgoing = 0; outgoing < n_mu; ++outgoing) {
            const size_t up = (size_t)level * (size_t)directions +
                              (size_t)(n_mu + 1 + outgoing);
            total_I[up] -= direct_seed_I[up];
            total_Q[up] -= direct_seed_Q[up];
            total_U[up] -= direct_seed_U[up];
        }
    }
}
