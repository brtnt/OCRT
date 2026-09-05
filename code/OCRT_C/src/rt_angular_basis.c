/*
 * OCRT angular basis and phase-kernel contractions.
 *
 * The scalar basis is the factorial-normalized associated Legendre function
 * without a Condon-Shortley phase.  The linear-polarization basis is the real
 * spin-2 pair used by the reciprocal three-Stokes expansion.  Both are built
 * from their defining three-term recurrences; the phase kernels are direct
 * contractions with the beta/gamma/alpha/zeta moment arrays.
 */

#include <math.h>
#include <stddef.h>
#include <string.h>

#include "rt_angular_basis.h"
#include "rt_kernel.h"

static double diagonal_normalization(int order)
{
    double value = 1.0;
    for (int i = 1; i <= order; ++i) {
        value *= sqrt((double)(i + order) / (double)i) * 0.5;
    }
    return value;
}

int rt_angular_basis_eval_scalar(int l_max,
                                 int fourier_order,
                                 double mu,
                                 double *scalar)
{
    if (!scalar || fourier_order < 0 || l_max < fourier_order || l_max < 2 ||
        !(mu >= -1.0 && mu <= 1.0)) return -1;

    for (int l = 0; l <= l_max; ++l) scalar[l] = 0.0;

    if (fourier_order == 0) {
        scalar[0] = 1.0;
        scalar[1] = mu;
        scalar[2] = 0.5 * (3.0 * mu * mu - 1.0);
    } else if (fourier_order == 1) {
        scalar[1] = sqrt(0.5 * (1.0 - mu * mu));
        scalar[2] = mu * scalar[1] * sqrt(3.0);
    } else {
        const double norm = diagonal_normalization(fourier_order);
        scalar[fourier_order - 1] = 0.0;
        scalar[fourier_order] = norm *
            pow(1.0 - mu * mu, 0.5 * (double)fourier_order);
    }

    const int first_degree = (fourier_order < 2) ? 2 : fourier_order;
    for (int l = first_degree; l < l_max; ++l) {
        const double forward = (2.0 * l + 1.0) /
            sqrt((double)(l + fourier_order + 1) *
                 (double)(l - fourier_order + 1));
        const double backward =
            sqrt((double)(l + fourier_order) *
                 (double)(l - fourier_order)) /
            (2.0 * l + 1.0);
        scalar[l + 1] = forward *
            (mu * scalar[l] - backward * scalar[l - 1]);
    }
    return 0;
}

int rt_angular_basis_eval_spin2(int l_max,
                                int fourier_order,
                                double mu,
                                double *linear_r,
                                double *linear_t)
{
    if (!linear_r || !linear_t || fourier_order < 0 ||
        l_max < fourier_order || l_max < 2 ||
        !(mu >= -1.0 && mu <= 1.0)) return -1;

    for (int l = 0; l <= l_max; ++l) {
        linear_r[l] = 0.0;
        linear_t[l] = 0.0;
    }

    if (fourier_order == 0) {
        linear_r[2] = 3.0 * (1.0 - mu * mu) / (2.0 * sqrt(6.0));
    } else if (fourier_order == 1) {
        const double root = sqrt(1.0 - mu * mu);
        linear_r[2] = -mu * root * 0.5;
        linear_t[2] = -root * 0.5;
    } else {
        const double scalar_norm = diagonal_normalization(fourier_order);
        const double spin_norm = scalar_norm *
            sqrt((double)fourier_order / (double)(fourier_order + 1)) *
            sqrt((double)(fourier_order - 1) /
                 (double)(fourier_order + 2));
        const double one_minus_mu2 = 1.0 - mu * mu;
        const double envelope =
            pow(one_minus_mu2, 0.5 * (double)fourier_order - 1.0);
        linear_r[fourier_order] =
            spin_norm * (1.0 + mu * mu) * envelope;
        linear_t[fourier_order] =
            spin_norm * 2.0 * mu * envelope;
    }

    const int first_degree = (fourier_order < 2) ? 2 : fourier_order;
    for (int l = first_degree; l < l_max; ++l) {
        const double scale = (double)(l + 1) * (double)(2 * l + 1) /
            sqrt((double)(l + 3) * (double)(l - 1) *
                 (double)(l + fourier_order + 1) *
                 (double)(l - fourier_order + 1));
        const double degree_coupling =
            sqrt((double)(l + 2) * (double)(l - 2) *
                 (double)(l + fourier_order) *
                 (double)(l - fourier_order)) /
            ((double)l * (double)(2 * l + 1));
        const double spin_coupling =
            2.0 * (double)fourier_order /
            ((double)l * (double)(l + 1));
        linear_r[l + 1] = scale *
            (mu * linear_r[l] - spin_coupling * linear_t[l] -
             degree_coupling * linear_r[l - 1]);
        linear_t[l + 1] = scale *
            (mu * linear_t[l] - spin_coupling * linear_r[l] -
             degree_coupling * linear_t[l - 1]);
    }
    return 0;
}

int rt_angular_basis_build_scalar(int n_mu,
                                  int l_max,
                                  int fourier_order,
                                  const double *mu,
                                  double **scalar)
{
    if (n_mu < 1 || !mu || !scalar || fourier_order < 0 ||
        l_max < fourier_order || l_max < 2) {
        return -1;
    }

    for (int l = 0; l < fourier_order; ++l) {
        for (int j = -n_mu; j <= n_mu; ++j) scalar[l][j] = 0.0;
    }

    if (fourier_order == 0) {
        for (int j = -n_mu; j <= n_mu; ++j) {
            const double x = mu[j];
            scalar[0][j] = 1.0;
            scalar[1][j] = x;
            scalar[2][j] = 0.5 * (3.0 * x * x - 1.0);
        }
    } else if (fourier_order == 1) {
        const double root3 = sqrt(3.0);
        for (int j = -n_mu; j <= n_mu; ++j) {
            const double x = mu[j];
            scalar[1][j] = sqrt(0.5 * (1.0 - x * x));
            scalar[2][j] = x * scalar[1][j] * root3;
        }
    } else {
        const double norm = diagonal_normalization(fourier_order);
        for (int j = -n_mu; j <= n_mu; ++j) {
            const double x = mu[j];
            const double one_minus_x2 = 1.0 - x * x;
            scalar[fourier_order - 1][j] = 0.0;
            scalar[fourier_order][j] =
                norm * pow(one_minus_x2, 0.5 * (double)fourier_order);
        }
    }

    const int first_degree = (fourier_order < 2) ? 2 : fourier_order;
    for (int l = first_degree; l < l_max; ++l) {
        const double forward = (2.0 * l + 1.0) /
            sqrt((double)(l + fourier_order + 1) *
                 (double)(l - fourier_order + 1));
        const double backward =
            sqrt((double)(l + fourier_order) *
                 (double)(l - fourier_order)) /
            (2.0 * l + 1.0);
        const double *previous = scalar[l - 1];
        const double *current = scalar[l];
        double *next = scalar[l + 1];

        for (int j = -n_mu; j <= n_mu; ++j) {
            next[j] = forward *
                (mu[j] * current[j] - backward * previous[j]);
        }
    }

    return 0;
}

int rt_angular_basis_build_spin2(int n_mu,
                                 int l_max,
                                 int fourier_order,
                                 const double *mu,
                                 double **linear_r,
                                 double **linear_t)
{
    if (n_mu < 1 || !mu || !linear_r || !linear_t || fourier_order < 0 ||
        l_max < fourier_order || l_max < 2) {
        return -1;
    }

    for (int l = 0; l <= l_max; ++l) {
        if (l < 2 || l < fourier_order) {
            for (int j = -n_mu; j <= n_mu; ++j) {
                linear_r[l][j] = 0.0;
                linear_t[l][j] = 0.0;
            }
        }
    }

    if (fourier_order == 0) {
        for (int j = -n_mu; j <= n_mu; ++j) {
            const double x = mu[j];
            const double one_minus_x2 = 1.0 - x * x;
            linear_r[2][j] = 3.0 * one_minus_x2 / (2.0 * sqrt(6.0));
            linear_t[2][j] = 0.0;
        }
    } else if (fourier_order == 1) {
        for (int j = -n_mu; j <= n_mu; ++j) {
            const double x = mu[j];
            const double root = sqrt(1.0 - x * x);
            linear_r[2][j] = -x * root * 0.5;
            linear_t[2][j] = -root * 0.5;
        }
    } else {
        const double scalar_norm = diagonal_normalization(fourier_order);
        const double spin_norm = scalar_norm *
            sqrt((double)fourier_order / (double)(fourier_order + 1)) *
            sqrt((double)(fourier_order - 1) /
                 (double)(fourier_order + 2));

        for (int j = -n_mu; j <= n_mu; ++j) {
            const double x = mu[j];
            const double one_minus_x2 = 1.0 - x * x;
            const double envelope =
                pow(one_minus_x2, 0.5 * (double)fourier_order - 1.0);
            linear_r[fourier_order][j] =
                spin_norm * (1.0 + x * x) * envelope;
            linear_t[fourier_order][j] =
                spin_norm * 2.0 * x * envelope;
        }
    }

    const int first_degree = (fourier_order < 2) ? 2 : fourier_order;
    for (int l = first_degree; l < l_max; ++l) {
        const double scale = (double)(l + 1) * (double)(2 * l + 1) /
            sqrt((double)(l + 3) * (double)(l - 1) *
                 (double)(l + fourier_order + 1) *
                 (double)(l - fourier_order + 1));
        const double degree_coupling =
            sqrt((double)(l + 2) * (double)(l - 2) *
                 (double)(l + fourier_order) *
                 (double)(l - fourier_order)) /
            ((double)l * (double)(2 * l + 1));
        const double spin_coupling =
            2.0 * (double)fourier_order /
            ((double)l * (double)(l + 1));

        const double *r_previous = linear_r[l - 1];
        const double *r_current = linear_r[l];
        double *r_next = linear_r[l + 1];
        const double *t_previous = linear_t[l - 1];
        const double *t_current = linear_t[l];
        double *t_next = linear_t[l + 1];

        for (int j = -n_mu; j <= n_mu; ++j) {
            const double x = mu[j];
            r_next[j] = scale *
                (x * r_current[j] - spin_coupling * t_current[j] -
                 degree_coupling * r_previous[j]);
            t_next[j] = scale *
                (x * t_current[j] - spin_coupling * r_current[j] -
                 degree_coupling * t_previous[j]);
        }
    }

    return 0;
}

void rt_angular_kernel_scalar(int n_mu,
                              int l_max,
                              int fourier_order,
                              double **scalar,
                              const double *beta,
                              double **out)
{
    for (int j = 0; j <= n_mu; ++j) {
        for (int k = -n_mu; k <= n_mu; ++k) {
            double value = 0.0;
            for (int l = fourier_order; l <= l_max; ++l) {
                value += scalar[l][j] * scalar[l][k] * beta[l];
            }
            out[j][k] = value;
        }
    }
}

void rt_angular_kernel_cross(int n_mu,
                             int l_max,
                             int fourier_order,
                             double **scalar,
                             double **linear_r,
                             double **linear_t,
                             const double *gamma,
                             double **out_r,
                             double **out_t)
{
    for (int j = 0; j <= n_mu; ++j) {
        for (int k = -n_mu; k <= n_mu; ++k) {
            double r_value = 0.0;
            double t_value = 0.0;
            for (int l = fourier_order; l <= l_max; ++l) {
                const double scalar_out = scalar[l][j];
                const double coupling = gamma[l];
                r_value += scalar_out * linear_r[l][k] * coupling;
                t_value += scalar_out * linear_t[l][k] * coupling;
            }
            out_r[j][k] = r_value;
            out_t[j][k] = t_value;
        }
    }
}

void rt_angular_kernel_linear(int n_mu,
                              int l_max,
                              int fourier_order,
                              double **linear_r,
                              double **linear_t,
                              const double *alpha,
                              const double *zeta,
                              double **out_rr,
                              double **out_rt,
                              double **out_tt)
{
    for (int j = 0; j <= n_mu; ++j) {
        for (int k = -n_mu; k <= n_mu; ++k) {
            double rr_value = 0.0;
            double rt_value = 0.0;
            double tt_value = 0.0;
            for (int l = fourier_order; l <= l_max; ++l) {
                const double rj = linear_r[l][j];
                const double tj = linear_t[l][j];
                const double rk = linear_r[l][k];
                const double tk = linear_t[l][k];
                const double a = alpha[l];
                const double z = zeta[l];
                const double tt = tj * tk;
                const double rr = rj * rk;

                tt_value += tt * a + rr * z;
                rr_value += tt * z + rr * a;
                rt_value += tj * rk * a + rj * tk * z;
            }
            out_rr[j][k] = rr_value;
            out_rt[j][k] = rt_value;
            out_tt[j][k] = tt_value;
        }
    }
}


static void invalidate_operator_pack(rt_legendre_workspace_t *ws)
{
    ws->sos_kernel_pack_valid = 0;
    ws->sos_kernel_pack_m = -1;
}

int rt_legendre_compute(rt_legendre_workspace_t *ws,
                        rt_atm_t *atm,
                        int m)
{
    if (!ws || !atm || atm->n_mu != ws->n_mu) return -1;
    invalidate_operator_pack(ws);
    if (rt_angular_basis_build_scalar(ws->n_mu, ws->l_max, m,
                                      atm->rm, ws->plm) != 0) {
        return -1;
    }
    for (int j = -ws->n_mu; j <= ws->n_mu; ++j) {
        atm->xpl[j] = ws->plm[2][j];
    }
    return 0;
}

int rt_legendre_compute_pol(rt_legendre_workspace_t *ws,
                            rt_atm_t *atm,
                            int m)
{
    if (rt_legendre_compute(ws, atm, m) != 0) return -1;
    if (rt_angular_basis_build_spin2(ws->n_mu, ws->l_max, m,
                                     atm->rm, ws->rrl, ws->rtl) != 0) {
        return -1;
    }
    for (int j = -ws->n_mu; j <= ws->n_mu; ++j) {
        atm->xrl[j] = ws->rrl[2][j];
        atm->xtl[j] = ws->rtl[2][j];
    }
    return 0;
}

int rt_kernel_phase_fourier(rt_legendre_workspace_t *ws,
                            int m,
                            const double *betal)
{
    if (!ws || m < 0 || m > ws->l_max) return -1;
    invalidate_operator_pack(ws);
    const size_t count = (size_t)(ws->n_mu + 1) *
                         (size_t)(2 * ws->n_mu + 1);
    if (!betal) {
        memset(ws->pfm_storage, 0, count * sizeof(double));
        return 0;
    }
    rt_angular_kernel_scalar(ws->n_mu, ws->l_max, m, ws->plm,
                             betal, ws->phase_fourier_m);
    return 0;
}

int rt_kernel_phase_fourier_pol(rt_legendre_workspace_t *ws,
                                int m,
                                const double *gammal)
{
    if (!ws || m < 0 || m > ws->l_max) return -1;
    invalidate_operator_pack(ws);
    const size_t count = (size_t)(ws->n_mu + 1) *
                         (size_t)(2 * ws->n_mu + 1);
    if (!gammal) {
        memset(ws->gr_storage, 0, count * sizeof(double));
        memset(ws->gt_storage, 0, count * sizeof(double));
        return 0;
    }
    rt_angular_kernel_cross(ws->n_mu, ws->l_max, m,
                            ws->plm, ws->rrl, ws->rtl, gammal,
                            ws->gr_pol, ws->gt_pol);
    return 0;
}

int rt_kernel_phase_fourier_aerosol_full(rt_legendre_workspace_t *ws,
                                         int m,
                                         const double *alphal,
                                         const double *zetal)
{
    if (!ws || m < 0 || m > ws->l_max) return -1;
    invalidate_operator_pack(ws);
    const size_t count = (size_t)(ws->n_mu + 1) *
                         (size_t)(2 * ws->n_mu + 1);
    if (!alphal || !zetal) {
        memset(ws->arr_storage, 0, count * sizeof(double));
        memset(ws->art_storage, 0, count * sizeof(double));
        memset(ws->att_storage, 0, count * sizeof(double));
        return 0;
    }
    rt_angular_kernel_linear(ws->n_mu, ws->l_max, m,
                             ws->rrl, ws->rtl, alphal, zetal,
                             ws->arr_pol, ws->art_pol, ws->att_pol);
    return 0;
}
