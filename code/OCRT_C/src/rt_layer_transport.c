/*
 * OCRT plane-parallel formal solution on a layered optical-depth grid.
 *
 * The source is represented by its values at layer interfaces.  Production
 * mode integrates the linear interpolant analytically; the constant mode is
 * retained only as a numerical diagnostic.  Upward and downward sweeps share
 * one signed-direction expression.
 */

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "rt_solver.h"

static int transport_cache_enabled_(void)
{
    /* Startup option: snapshot once per worker.  getenv() itself can serialize
     * in some C runtimes, and this function is reached for every Stokes/formal
     * integration at every scattering order. */
    static _Thread_local int enabled = -1;
    if (enabled < 0) {
        const char *off = getenv("OCRT_TRANSPORT_CACHE_OFF");
        enabled = !(off && off[0] && strcmp(off, "0") != 0 &&
                    strcmp(off, "off") != 0 && strcmp(off, "false") != 0);
    }
    return enabled;
}

/* Build the immutable Beer-Lambert transmission stencil once per finalized
 * medium.  h[] and rm[] snapshots make the lazy cache safe for view-node and
 * gas-absorption paths that modify an already allocated rt_atm_t before the
 * first (or a later) formal-solution call.  A failed allocation merely falls
 * back to the historical per-call exp() path. */
static int transport_cache_prepare_(rt_atm_t *atm)
{
    if (!atm || !atm->h || !atm->rm || atm->n_layers < 1 || atm->n_mu < 1)
        return -1;
    if (!transport_cache_enabled_()) return 1;

    const int n_layers = atm->n_layers;
    const int n_mu = atm->n_mu;
    const int directions = 2 * n_mu + 1;
    const size_t h_bytes = (size_t)(n_layers + 1) * sizeof(double);
    const size_t rm_bytes = (size_t)directions * sizeof(double);
    const size_t t_count = (size_t)directions * (size_t)n_layers;

    if (atm->transport_cache_ready && atm->transport_T &&
        atm->transport_h_key && atm->transport_rm_key &&
        memcmp(atm->transport_h_key, atm->h, h_bytes) == 0 &&
        memcmp(atm->transport_rm_key, &atm->rm[-n_mu], rm_bytes) == 0) {
        ++atm->transport_cache_hits;
        return 0;
    }

    if (!atm->transport_T)
        atm->transport_T = (double *)malloc(t_count * sizeof(double));
    if (!atm->transport_h_key)
        atm->transport_h_key = (double *)malloc(h_bytes);
    if (!atm->transport_rm_key)
        atm->transport_rm_key = (double *)malloc(rm_bytes);
    if (!atm->transport_T || !atm->transport_h_key ||
        !atm->transport_rm_key) {
        free(atm->transport_T); atm->transport_T = NULL;
        free(atm->transport_h_key); atm->transport_h_key = NULL;
        free(atm->transport_rm_key); atm->transport_rm_key = NULL;
        atm->transport_cache_ready = 0;
        return 1;
    }

    for (int direction = -n_mu; direction <= n_mu; ++direction) {
        if (direction == 0) continue;
        const double abs_mu = fabs(atm->rm[direction]);
        const size_t row = (size_t)(direction + n_mu) * (size_t)n_layers;
        for (int layer = 0; layer < n_layers; ++layer) {
            const double delta_tau = atm->h[layer + 1] - atm->h[layer];
            atm->transport_T[row + (size_t)layer] =
                exp(-delta_tau / abs_mu);
        }
    }
    memcpy(atm->transport_h_key, atm->h, h_bytes);
    memcpy(atm->transport_rm_key, &atm->rm[-n_mu], rm_bytes);
    atm->transport_cache_ready = 1;
    ++atm->transport_cache_builds;
    return 0;
}

static double linear_source_integral(double slope,
                                     double intercept,
                                     double signed_mu,
                                     double tau_out,
                                     double tau_in,
                                     double transmission)
{
    const double one_minus_transmission = 1.0 - transmission;
    return (one_minus_transmission * (intercept + slope * signed_mu) +
            slope * (tau_out - tau_in * transmission)) * 0.5;
}

static double constant_source_integral(double source,
                                       double transmission)
{
    return (1.0 - transmission) * source * 0.5;
}

int rt_solver_integrate_bcs(const rt_atm_t *atm,
                            int m,
                            const double *source_first,
                            const double *surface_bc,
                            const double *top_down_bc,
                            rt_integration_method_t method,
                            double *radiation)
{
    if (!atm || !source_first || !radiation || m < 0) return -1;
    if (method != RT_INTEGRATION_METHOD_LINEAR &&
        method != RT_INTEGRATION_METHOD_CONSTANT) return -1;

    const int n_layers = atm->n_layers;
    const int n_mu = atm->n_mu;
    const int directions = 2 * n_mu + 1;
    const double *tau = atm->h;
    rt_atm_t *mutable_atm = (rt_atm_t *)(void *)atm;
    const int transport_cache_status = transport_cache_prepare_(mutable_atm);
    const int use_transport_cache = (transport_cache_status == 0);

    for (size_t i = 0;
         i < (size_t)(n_layers + 1) * (size_t)directions;
         ++i) {
        radiation[i] = NAN;
    }

    for (int direction = -n_mu; direction <= n_mu; ++direction) {
        if (direction == 0) continue;

        const size_t column = (size_t)(direction + n_mu);
        const double signed_mu = atm->rm[direction];
        const double abs_mu = fabs(signed_mu);

        if (direction > 0) {
            double running = surface_bc ? surface_bc[direction - 1] : 0.0;
            radiation[(size_t)n_layers * (size_t)directions + column] = running;

            for (int level = n_layers - 1; level >= 0; --level) {
                const int upstream = level + 1;
                const double source_out =
                    source_first[(size_t)level * (size_t)directions + column];
                const double source_in =
                    source_first[(size_t)upstream * (size_t)directions + column];
                if (!isfinite(source_out) || !isfinite(source_in)) return -2;

                const double delta_tau = tau[upstream] - tau[level];
                const double transmission = use_transport_cache
                    ? atm->transport_T[
                          (size_t)(direction + n_mu) * (size_t)n_layers +
                          (size_t)level]
                    : exp(-delta_tau / abs_mu);
                double contribution;
                if (method == RT_INTEGRATION_METHOD_LINEAR) {
                    const double slope = (source_in - source_out) / delta_tau;
                    const double intercept = source_out - slope * tau[level];
                    contribution = linear_source_integral(
                        slope, intercept, signed_mu,
                        tau[level], tau[upstream], transmission);
                } else {
                    contribution = constant_source_integral(source_out,
                                                            transmission);
                }
                running = transmission * running + contribution;
                radiation[(size_t)level * (size_t)directions + column] = running;
            }
        } else {
            double running = top_down_bc ? top_down_bc[-direction - 1] : 0.0;
            radiation[column] = running;

            for (int level = 1; level <= n_layers; ++level) {
                const int upstream = level - 1;
                const double source_out =
                    source_first[(size_t)level * (size_t)directions + column];
                const double source_in =
                    source_first[(size_t)upstream * (size_t)directions + column];
                if (!isfinite(source_out) || !isfinite(source_in)) return -2;

                const double delta_tau = tau[level] - tau[upstream];
                const double transmission = use_transport_cache
                    ? atm->transport_T[
                          (size_t)(direction + n_mu) * (size_t)n_layers +
                          (size_t)upstream]
                    : exp(-delta_tau / abs_mu);
                double contribution;
                if (method == RT_INTEGRATION_METHOD_LINEAR) {
                    const double slope = (source_out - source_in) / delta_tau;
                    const double intercept = source_out - slope * tau[level];
                    contribution = linear_source_integral(
                        slope, intercept, signed_mu,
                        tau[level], tau[upstream], transmission);
                } else {
                    contribution = constant_source_integral(source_out,
                                                            transmission);
                }
                running = transmission * running + contribution;
                radiation[(size_t)level * (size_t)directions + column] = running;
            }
        }
    }

    return 0;
}

int rt_solver_integrate_with_surface_bc(const rt_atm_t *atm,
                                         int m,
                                         const double *source_first,
                                         const double *surface_bc,
                                         rt_integration_method_t method,
                                         double *radiation)
{
    return rt_solver_integrate_bcs(atm, m, source_first, surface_bc, NULL,
                                   method, radiation);
}

int rt_solver_integrate(const rt_atm_t *atm,
                        int m,
                        const double *source_first,
                        rt_integration_method_t method,
                        double *radiation)
{
    return rt_solver_integrate_with_surface_bc(atm, m, source_first, NULL,
                                                method, radiation);
}
