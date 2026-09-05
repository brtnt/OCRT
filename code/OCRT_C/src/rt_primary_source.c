/*
 * OCRT direct-beam single-scattering source for one azimuthal Fourier mode.
 *
 * Molecular and particulate phase operators are mixed by the layer scattering
 * fractions xdel/ydel.  The stored direct-beam attenuation already carries
 * the OCRT source normalization, so this module performs no reference-code
 * compatibility adjustment.
 */

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "rt_solver.h"

int rt_solver_primary_source(const rt_atm_t *atm,
                             int m,
                             const rt_legendre_workspace_t *ws,
                             double *source_first)
{
    if (!atm || !ws || !source_first || m < 0 || m > ws->l_max ||
        atm->n_mu != ws->n_mu) return -1;

    const int n_layers = atm->n_layers;
    const int n_mu = atm->n_mu;
    const int directions = 2 * n_mu + 1;
    const double beta0 = (m == 0) ? atm->beta0 : 0.0;
    const double beta2 = atm->beta2;
    const double beam_scalar_l2 = atm->xpl[0];
    const double beam_q = atm->beam_q;
    const int molecular_active = (m <= 2);
    const double gamma2 = atm->gamma2;
    const double beam_r_l2 = atm->xrl ? atm->xrl[0] : 0.0;

    double *beam_q_to_i_particle = NULL;
    if (beam_q != 0.0) {
        beam_q_to_i_particle =
            (double *)calloc((size_t)directions, sizeof(double));
        if (!beam_q_to_i_particle) return -1;

        if (ws->beam_q_direct_valid && ws->beam_q_to_i_direct) {
            memcpy(beam_q_to_i_particle, ws->beam_q_to_i_direct,
                   (size_t)directions * sizeof(double));
        } else if (atm->gammal_aer && ws->plm && ws->rrl) {
            for (int direction = -n_mu; direction <= n_mu; ++direction) {
                double value = 0.0;
                for (int l = m; l <= ws->l_max; ++l) {
                    value += ws->plm[l][direction] * ws->rrl[l][0] *
                             atm->gammal_aer[l];
                }
                beam_q_to_i_particle[direction + n_mu] = value;
            }
        }
    }

    for (int level = 0; level <= n_layers; ++level) {
        const double beam = atm->ch[level];
        const double particle_fraction = atm->xdel[level];
        const double molecular_fraction = atm->ydel[level];
        double *row = source_first +
                      (size_t)level * (size_t)directions;

        for (int direction = -n_mu; direction <= n_mu; ++direction) {
            const double molecular_phase =
                beta0 + beta2 * atm->xpl[direction] * beam_scalar_l2;
            const double particle_phase =
                ws->phase_fourier_m[0][direction];
            double value = beam *
                (particle_fraction * particle_phase +
                 molecular_fraction * molecular_phase);

            if (beam_q != 0.0) {
                const double molecular_q_to_i = molecular_active
                    ? gamma2 * atm->xpl[direction] * beam_r_l2 : 0.0;
                value += beam * beam_q *
                    (particle_fraction *
                         beam_q_to_i_particle[direction + n_mu] +
                     molecular_fraction * molecular_q_to_i);
            }
            row[direction + n_mu] = value;
        }
        row[n_mu] = NAN;
    }

    free(beam_q_to_i_particle);
    return 0;
}

int rt_solver_primary_source_pol(const rt_atm_t *atm,
                                 int m,
                                 const rt_legendre_workspace_t *ws,
                                 double *source_q,
                                 double *source_u)
{
    if (!atm || !ws || !source_q || !source_u || m < 0 ||
        m > ws->l_max || atm->n_mu != ws->n_mu ||
        !ws->gr_pol || !ws->gt_pol) return -1;

    const int n_layers = atm->n_layers;
    const int n_mu = atm->n_mu;
    const int directions = 2 * n_mu + 1;
    const double gamma2 = atm->gamma2;
    const double beam_scalar_l2 = atm->xpl[0];
    const int molecular_active = (m <= 2);
    const double beam_q = atm->beam_q;
    const double alpha2 = atm->alpha2;
    const double beam_r_l2 = atm->xrl ? atm->xrl[0] : 0.0;

    double *beam_q_to_q_particle = NULL;
    double *beam_q_to_u_particle = NULL;
    if (beam_q != 0.0) {
        beam_q_to_q_particle =
            (double *)calloc((size_t)directions, sizeof(double));
        beam_q_to_u_particle =
            (double *)calloc((size_t)directions, sizeof(double));
        if (!beam_q_to_q_particle || !beam_q_to_u_particle) {
            free(beam_q_to_q_particle);
            free(beam_q_to_u_particle);
            return -1;
        }

        if (ws->beam_q_direct_valid && ws->beam_q_to_q_direct &&
            ws->beam_q_to_u_direct) {
            memcpy(beam_q_to_q_particle, ws->beam_q_to_q_direct,
                   (size_t)directions * sizeof(double));
            memcpy(beam_q_to_u_particle, ws->beam_q_to_u_direct,
                   (size_t)directions * sizeof(double));
        } else if (atm->alphal_aer && atm->zetal_aer && ws->rrl && ws->rtl) {
            for (int direction = -n_mu; direction <= n_mu; ++direction) {
                double q_value = 0.0;
                double u_value = 0.0;
                for (int l = m; l <= ws->l_max; ++l) {
                    const double rj = ws->rrl[l][direction];
                    const double tj = ws->rtl[l][direction];
                    const double r0 = ws->rrl[l][0];
                    const double t0 = ws->rtl[l][0];
                    const double alpha = atm->alphal_aer[l];
                    const double zeta = atm->zetal_aer[l];
                    q_value += t0 * tj * zeta + r0 * rj * alpha;
                    u_value += tj * r0 * alpha + rj * t0 * zeta;
                }
                beam_q_to_q_particle[direction + n_mu] = q_value;
                beam_q_to_u_particle[direction + n_mu] = u_value;
            }
        }
    }

    for (int level = 0; level <= n_layers; ++level) {
        const double beam = atm->ch[level];
        const double particle_fraction = atm->xdel[level];
        const double molecular_fraction = atm->ydel[level];
        double *q_row = source_q +
                        (size_t)level * (size_t)directions;
        double *u_row = source_u +
                        (size_t)level * (size_t)directions;

        for (int direction = -n_mu; direction <= n_mu; ++direction) {
            const double molecular_i_to_q = molecular_active
                ? gamma2 * atm->xrl[direction] * beam_scalar_l2 : 0.0;
            const double molecular_i_to_u = molecular_active
                ? gamma2 * atm->xtl[direction] * beam_scalar_l2 : 0.0;
            const double particle_i_to_q = ws->gr_pol[0][direction];
            const double particle_i_to_u = ws->gt_pol[0][direction];

            double q_value = beam *
                (particle_fraction * particle_i_to_q +
                 molecular_fraction * molecular_i_to_q);
            double u_value = -beam *
                (particle_fraction * particle_i_to_u +
                 molecular_fraction * molecular_i_to_u);

            if (beam_q != 0.0) {
                const double molecular_q_to_q = molecular_active
                    ? alpha2 * atm->xrl[direction] * beam_r_l2 : 0.0;
                const double molecular_q_to_u = molecular_active
                    ? alpha2 * atm->xtl[direction] * beam_r_l2 : 0.0;
                q_value += beam * beam_q *
                    (particle_fraction *
                         beam_q_to_q_particle[direction + n_mu] +
                     molecular_fraction * molecular_q_to_q);
                u_value += -beam * beam_q *
                    (particle_fraction *
                         beam_q_to_u_particle[direction + n_mu] +
                     molecular_fraction * molecular_q_to_u);
            }

            q_row[direction + n_mu] = q_value;
            u_row[direction + n_mu] = u_value;
        }
        q_row[n_mu] = NAN;
        u_row[n_mu] = NAN;
    }

    free(beam_q_to_q_particle);
    free(beam_q_to_u_particle);
    return 0;
}
