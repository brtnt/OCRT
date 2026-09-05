/*
 * OCRT successive-scattering source operator.
 *
 * The phase operator is separated into particle and molecular direction-pair
 * coefficients.  Those coefficients depend on Fourier order and quadrature,
 * but not on scattering order, and are therefore packed once.  Each SOS step
 * then performs only weighted I/Q/U contractions and no allocation when the
 * prepared path is used.
 */

#include "rt_sos_operator.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Exact optimization gate.  Requested/view-only ordinates carry a zero
 * quadrature weight and therefore contribute identically zero to the source
 * contraction.  The environment switch exists only for byte-regression and
 * performance diagnosis; production keeps the skip enabled. */
static int ocrt_sos_zero_col_skip_enabled(void)
{
    const char *off = getenv("OCRT_SOS_ZERO_COL_SKIP_OFF");
    return !(off && off[0] && strcmp(off, "0") != 0);
}

/* Eighteen coefficients per component.  The particle block occupies 0..17;
 * the molecular block uses the same ordering at 18..35. */
enum {
    K_II_PP = 0, K_II_PM,
    K_IU_PP, K_IU_PM, K_IU_TP, K_IU_TM,
    K_IQ_PP, K_IQ_PM, K_IQ_TP, K_IQ_TM,
    K_QQ_PP, K_QQ_PM,
    K_QU_PP, K_QU_PM, K_QU_TP, K_QU_TM,
    K_UU_PP, K_UU_PM,
    K_COMPONENT_COUNT
};

static size_t pair_width(const rt_atm_t *atm)
{
    return (size_t)(atm->n_mu + 1);
}

static size_t packed_row_size(const rt_atm_t *atm)
{
    return (size_t)(2 * K_COMPONENT_COUNT) * pair_width(atm);
}

static void fill_pair_operator(const rt_atm_t *atm,
                               int m,
                               const rt_legendre_workspace_t *ws,
                               int outgoing,
                               double *packed)
{
    const int n_mu = atm->n_mu;
    const size_t width = pair_width(atm);
    const int molecular_active = (m <= 2);
    const double beta0 = (m == 0) ? atm->beta0 : 0.0;
    const double beta2 = atm->beta2;
    const double gamma2 = atm->gamma2;
    const double alpha2 = atm->alpha2;

    const double p_out = atm->xpl[+outgoing];
    const double r_out = atm->xrl[+outgoing];
    const double t_out = atm->xtl[+outgoing];
    const double p_out_mirror = atm->xpl[-outgoing];
    const double r_out_mirror = atm->xrl[-outgoing];
    const double t_out_mirror = atm->xtl[-outgoing];

    for (int incoming = 1; incoming <= n_mu; ++incoming) {
        const double p_in = atm->xpl[+incoming];
        const double r_in = atm->xrl[+incoming];
        const double t_in = atm->xtl[+incoming];
        const double r_in_mirror = atm->xrl[-incoming];
        const double t_in_mirror = atm->xtl[-incoming];

        packed[K_II_PP * width + incoming] =
            ws->phase_fourier_m[+incoming][+outgoing];
        packed[K_II_PM * width + incoming] =
            ws->phase_fourier_m[+incoming][-outgoing];
        packed[K_IU_PP * width + incoming] = ws->gt_pol[+incoming][+outgoing];
        packed[K_IU_PM * width + incoming] = ws->gt_pol[+incoming][-outgoing];
        packed[K_IU_TP * width + incoming] = ws->gt_pol[+outgoing][+incoming];
        packed[K_IU_TM * width + incoming] = ws->gt_pol[+outgoing][-incoming];
        packed[K_IQ_PP * width + incoming] = ws->gr_pol[+incoming][+outgoing];
        packed[K_IQ_PM * width + incoming] = ws->gr_pol[+incoming][-outgoing];
        packed[K_IQ_TP * width + incoming] = ws->gr_pol[+outgoing][+incoming];
        packed[K_IQ_TM * width + incoming] = ws->gr_pol[+outgoing][-incoming];
        packed[K_QQ_PP * width + incoming] = ws->arr_pol[+incoming][+outgoing];
        packed[K_QQ_PM * width + incoming] = ws->arr_pol[+incoming][-outgoing];
        packed[K_QU_PP * width + incoming] = ws->art_pol[+incoming][+outgoing];
        packed[K_QU_PM * width + incoming] = ws->art_pol[+incoming][-outgoing];
        packed[K_QU_TP * width + incoming] = ws->art_pol[+outgoing][+incoming];
        packed[K_QU_TM * width + incoming] = ws->art_pol[+outgoing][-incoming];
        packed[K_UU_PP * width + incoming] = ws->att_pol[+incoming][+outgoing];
        packed[K_UU_PM * width + incoming] = ws->att_pol[+incoming][-outgoing];

        double *molecular = packed + (size_t)K_COMPONENT_COUNT * width;
        if (!molecular_active) {
            for (int channel = 0; channel < K_COMPONENT_COUNT; ++channel) {
                molecular[(size_t)channel * width + incoming] = 0.0;
            }
            continue;
        }

        molecular[K_II_PP * width + incoming] = beta0 + beta2 * p_in * p_out;
        molecular[K_II_PM * width + incoming] = beta0 + beta2 * p_in * p_out_mirror;
        molecular[K_IU_PP * width + incoming] = gamma2 * p_in * t_out;
        molecular[K_IU_PM * width + incoming] = gamma2 * p_in * t_out_mirror;
        molecular[K_IU_TP * width + incoming] = gamma2 * p_out * t_in;
        molecular[K_IU_TM * width + incoming] = gamma2 * p_out * t_in_mirror;
        molecular[K_IQ_PP * width + incoming] = gamma2 * p_in * r_out;
        molecular[K_IQ_PM * width + incoming] = gamma2 * p_in * r_out_mirror;
        molecular[K_IQ_TP * width + incoming] = gamma2 * p_out * r_in;
        molecular[K_IQ_TM * width + incoming] = gamma2 * p_out * r_in_mirror;
        molecular[K_QQ_PP * width + incoming] = alpha2 * r_in * r_out;
        molecular[K_QQ_PM * width + incoming] = alpha2 * r_in * r_out_mirror;
        molecular[K_QU_PP * width + incoming] = alpha2 * t_in * r_out;
        molecular[K_QU_PM * width + incoming] = alpha2 * t_in * r_out_mirror;
        molecular[K_QU_TP * width + incoming] = alpha2 * t_out * r_in;
        molecular[K_QU_TM * width + incoming] = alpha2 * t_out * r_in_mirror;
        molecular[K_UU_PP * width + incoming] = alpha2 * t_in * t_out;
        molecular[K_UU_PM * width + incoming] = alpha2 * t_in * t_out_mirror;
    }
}

int rt_sos_operator_prepare(const rt_atm_t *atm,
                            int m,
                            rt_legendre_workspace_t *ws)
{
    if (!atm || !ws || atm->n_mu != ws->n_mu || m < 0 ||
        !ws->phase_fourier_m || !ws->gr_pol || !ws->gt_pol ||
        !ws->arr_pol || !ws->art_pol || !ws->att_pol) {
        return -1;
    }
    if (ws->sos_kernel_pack_valid && ws->sos_kernel_pack_m == m) {
        return 0;
    }

    const size_t levels = (size_t)(atm->n_layers + 1);
    const size_t pack_per_outgoing = packed_row_size(atm);
    const size_t pack_count = (size_t)atm->n_mu * pack_per_outgoing;
    const size_t transpose_count = (size_t)6 * (size_t)atm->n_mu * levels;
    const size_t accumulator_count = (size_t)6 * levels;
    const size_t arena_count = transpose_count + accumulator_count;

    if (ws->sos_kernel_pack_capacity < pack_count) {
        double *replacement = (double *)realloc(ws->sos_kernel_pack,
                                                 pack_count * sizeof(double));
        if (!replacement) return -3;
        ws->sos_kernel_pack = replacement;
        ws->sos_kernel_pack_capacity = pack_count;
    }
    if (ws->sos_source_arena_capacity < arena_count) {
        double *replacement = (double *)realloc(ws->sos_source_arena,
                                                 arena_count * sizeof(double));
        if (!replacement) return -3;
        ws->sos_source_arena = replacement;
        ws->sos_source_arena_capacity = arena_count;
    }

    for (int outgoing = 1; outgoing <= atm->n_mu; ++outgoing) {
        double *row = ws->sos_kernel_pack +
                      (size_t)(outgoing - 1) * pack_per_outgoing;
        fill_pair_operator(atm, m, ws, outgoing, row);
    }

    ws->sos_kernel_pack_m = m;
    ws->sos_kernel_pack_valid = 1;

    const char *trace = getenv("OCRT_DUMP_SOS_PHASE_CACHE");
    if (trace && trace[0] && strcmp(trace, "0") != 0) {
        fprintf(stderr,
                "SOS_PHASE_CACHE_BUILD m=%d n_mu=%d n_layers=%d "
                "pack_bytes=%zu scratch_bytes=%zu\n",
                m, atm->n_mu, atm->n_layers,
                pack_count * sizeof(double),
                arena_count * sizeof(double));
    }
    return 0;
}

static int validate_vector_field(const rt_atm_t *atm,
                                 const double *I,
                                 const double *Q,
                                 const double *U)
{
    const int n_mu = atm->n_mu;
    const int directions = 2 * n_mu + 1;
    for (int level = 0; level <= atm->n_layers; ++level) {
        const size_t row = (size_t)level * (size_t)directions;
        for (int incoming = 1; incoming <= n_mu; ++incoming) {
            const size_t up = row + (size_t)(n_mu + incoming);
            const size_t down = row + (size_t)(n_mu - incoming);
            if (!isfinite(I[up]) || !isfinite(I[down]) ||
                !isfinite(Q[up]) || !isfinite(Q[down]) ||
                !isfinite(U[up]) || !isfinite(U[down])) {
                return -2;
            }
        }
    }
    return 0;
}

int rt_sos_operator_apply_vector(const rt_atm_t *atm,
                                 int m,
                                 const rt_legendre_workspace_t *ws,
                                 const double *restrict I_previous,
                                 const double *restrict Q_previous,
                                 const double *restrict U_previous,
                                 double *restrict J_I,
                                 double *restrict J_Q,
                                 double *restrict J_U)
{
    if (!atm || !ws || !I_previous || !Q_previous || !U_previous ||
        !J_I || !J_Q || !J_U || atm->n_mu != ws->n_mu || m < 0) {
        return -1;
    }

    const int n_mu = atm->n_mu;
    const int directions = 2 * n_mu + 1;
    const size_t levels = (size_t)(atm->n_layers + 1);
    const size_t field_count = levels * (size_t)directions;
    for (size_t index = 0; index < field_count; ++index) {
        J_I[index] = NAN;
        J_Q[index] = NAN;
        J_U[index] = NAN;
    }

    const int finite_status = validate_vector_field(atm, I_previous, Q_previous,
                                                     U_previous);
    if (finite_status != 0) return finite_status;

    const size_t width = pair_width(atm);
    const size_t pack_per_outgoing = packed_row_size(atm);
    const size_t transpose_count = (size_t)6 * (size_t)n_mu * levels;
    const size_t accumulator_count = (size_t)6 * levels;
    const int prepared = ws->sos_kernel_pack_valid &&
                         ws->sos_kernel_pack_m == m && ws->sos_kernel_pack &&
                         ws->sos_kernel_pack_capacity >=
                             (size_t)n_mu * pack_per_outgoing &&
                         ws->sos_source_arena &&
                         ws->sos_source_arena_capacity >=
                             transpose_count + accumulator_count;

    double *owned = NULL;
    double *arena = prepared ? ws->sos_source_arena :
        (owned = (double *)malloc((transpose_count + accumulator_count +
                                  pack_per_outgoing) * sizeof(double)));
    if (!arena) return -3;

    double *transpose = arena;
    double *accumulator = arena + transpose_count;
    double *fallback_pack = prepared ? NULL : accumulator + accumulator_count;

#define TRANSPOSE_PLANE(component, incoming) \
    (transpose + ((size_t)(component) * (size_t)n_mu + \
                  (size_t)((incoming) - 1)) * levels)

    for (int incoming = 1; incoming <= n_mu; ++incoming) {
        const size_t up_column = (size_t)(n_mu + incoming);
        const size_t down_column = (size_t)(n_mu - incoming);
        double *I_up = TRANSPOSE_PLANE(0, incoming);
        double *I_down = TRANSPOSE_PLANE(1, incoming);
        double *Q_up = TRANSPOSE_PLANE(2, incoming);
        double *Q_down = TRANSPOSE_PLANE(3, incoming);
        double *U_up = TRANSPOSE_PLANE(4, incoming);
        double *U_down = TRANSPOSE_PLANE(5, incoming);

        for (size_t level = 0; level < levels; ++level) {
            const size_t row = level * (size_t)directions;
            I_up[level] = I_previous[row + up_column];
            I_down[level] = I_previous[row + down_column];
            Q_up[level] = Q_previous[row + up_column];
            Q_down[level] = Q_previous[row + down_column];
            U_up[level] = U_previous[row + up_column];
            U_down[level] = U_previous[row + down_column];
        }
    }

    const double *restrict particle_weight = atm->xdel;
    const double *restrict molecular_weight = atm->ydel;
    const int skip_zero_columns = ocrt_sos_zero_col_skip_enabled();

    for (int outgoing = 1; outgoing <= n_mu; ++outgoing) {
        const double *packed = prepared ?
            ws->sos_kernel_pack + (size_t)(outgoing - 1) * pack_per_outgoing :
            fallback_pack;
        if (!prepared) fill_pair_operator(atm, m, ws, outgoing, fallback_pack);

        double *I_out_up = accumulator + 0 * levels;
        double *I_out_down = accumulator + 1 * levels;
        double *Q_out_up = accumulator + 2 * levels;
        double *Q_out_down = accumulator + 3 * levels;
        double *U_out_up = accumulator + 4 * levels;
        double *U_out_down = accumulator + 5 * levels;
        memset(accumulator, 0, accumulator_count * sizeof(double));

        for (int incoming = 1; incoming <= n_mu; ++incoming) {
            const double quadrature_weight = atm->gb[+incoming];
            if (skip_zero_columns && quadrature_weight == 0.0) continue;
            const double *particle = packed;
            const double *molecular = packed + (size_t)K_COMPONENT_COUNT * width;

#define COEF(block, channel) ((block)[(size_t)(channel) * width + (size_t)incoming])
            const double p0 = COEF(particle, K_II_PP);
            const double p1 = COEF(particle, K_II_PM);
            const double p2 = COEF(particle, K_IU_PP);
            const double p3 = COEF(particle, K_IU_PM);
            const double p4 = COEF(particle, K_IU_TP);
            const double p5 = COEF(particle, K_IU_TM);
            const double p6 = COEF(particle, K_IQ_PP);
            const double p7 = COEF(particle, K_IQ_PM);
            const double p8 = COEF(particle, K_IQ_TP);
            const double p9 = COEF(particle, K_IQ_TM);
            const double p10 = COEF(particle, K_QQ_PP);
            const double p11 = COEF(particle, K_QQ_PM);
            const double p12 = COEF(particle, K_QU_PP);
            const double p13 = COEF(particle, K_QU_PM);
            const double p14 = COEF(particle, K_QU_TP);
            const double p15 = COEF(particle, K_QU_TM);
            const double p16 = COEF(particle, K_UU_PP);
            const double p17 = COEF(particle, K_UU_PM);

            const double r0 = COEF(molecular, K_II_PP);
            const double r1 = COEF(molecular, K_II_PM);
            const double r2 = COEF(molecular, K_IU_PP);
            const double r3 = COEF(molecular, K_IU_PM);
            const double r4 = COEF(molecular, K_IU_TP);
            const double r5 = COEF(molecular, K_IU_TM);
            const double r6 = COEF(molecular, K_IQ_PP);
            const double r7 = COEF(molecular, K_IQ_PM);
            const double r8 = COEF(molecular, K_IQ_TP);
            const double r9 = COEF(molecular, K_IQ_TM);
            const double r10 = COEF(molecular, K_QQ_PP);
            const double r11 = COEF(molecular, K_QQ_PM);
            const double r12 = COEF(molecular, K_QU_PP);
            const double r13 = COEF(molecular, K_QU_PM);
            const double r14 = COEF(molecular, K_QU_TP);
            const double r15 = COEF(molecular, K_QU_TM);
            const double r16 = COEF(molecular, K_UU_PP);
            const double r17 = COEF(molecular, K_UU_PM);
#undef COEF

            const double *restrict I_up = TRANSPOSE_PLANE(0, incoming);
            const double *restrict I_down = TRANSPOSE_PLANE(1, incoming);
            const double *restrict Q_up = TRANSPOSE_PLANE(2, incoming);
            const double *restrict Q_down = TRANSPOSE_PLANE(3, incoming);
            const double *restrict U_up = TRANSPOSE_PLANE(4, incoming);
            const double *restrict U_down = TRANSPOSE_PLANE(5, incoming);

#pragma omp simd
            for (size_t level = 0; level < levels; ++level) {
                const double x = particle_weight[level];
                const double y = molecular_weight[level];
                const double ii_same = p0 * x + y * r0;
                const double ii_cross = p1 * x + y * r1;
                const double iu_same = p2 * x + y * r2;
                const double iu_cross = p3 * x + y * r3;
                const double iu_transpose_same = p4 * x + y * r4;
                const double iu_transpose_cross = p5 * x + y * r5;
                const double iq_same = p6 * x + y * r6;
                const double iq_cross = p7 * x + y * r7;
                const double iq_transpose_same = p8 * x + y * r8;
                const double iq_transpose_cross = p9 * x + y * r9;
                const double qq_same = p10 * x + y * r10;
                const double qq_cross = p11 * x + y * r11;
                const double qu_same = p12 * x + y * r12;
                const double qu_cross = p13 * x + y * r13;
                const double qu_transpose_same = p14 * x + y * r14;
                const double qu_transpose_cross = p15 * x + y * r15;
                const double uu_same = p16 * x + y * r16;
                const double uu_cross = p17 * x + y * r17;

                const double i_up = I_up[level];
                const double i_down = I_down[level];
                const double q_up = Q_up[level];
                const double q_down = Q_down[level];
                const double u_up = U_up[level];
                const double u_down = U_down[level];

                I_out_up[level] += quadrature_weight *
                    (i_up * ii_same + q_up * iq_transpose_same - u_up * iu_transpose_same) +
                    quadrature_weight *
                    (i_down * ii_cross + q_down * iq_transpose_cross - u_down * iu_transpose_cross);
                I_out_down[level] += quadrature_weight *
                    (i_up * ii_cross + q_up * iq_transpose_cross + u_up * iu_transpose_cross) +
                    quadrature_weight *
                    (i_down * ii_same + q_down * iq_transpose_same + u_down * iu_transpose_same);
                Q_out_up[level] += quadrature_weight *
                    (i_up * iq_same + q_up * qq_same - u_up * qu_same) +
                    quadrature_weight *
                    (i_down * iq_cross + q_down * qq_cross + u_down * qu_cross);
                Q_out_down[level] += quadrature_weight *
                    (i_up * iq_cross + q_up * qq_cross - u_up * qu_cross) +
                    quadrature_weight *
                    (i_down * iq_same + q_down * qq_same + u_down * qu_same);
                U_out_up[level] += -quadrature_weight *
                    (i_up * iu_same + q_up * qu_transpose_same - u_up * uu_same) +
                    -quadrature_weight *
                    (-i_down * iu_cross + q_down * qu_transpose_cross - u_down * uu_cross);
                U_out_down[level] += -quadrature_weight *
                    (i_up * iu_cross - q_up * qu_transpose_cross - u_up * uu_cross) +
                    -quadrature_weight *
                    (-i_down * iu_same - q_down * qu_transpose_same - u_down * uu_same);
            }
        }

        for (size_t level = 0; level < levels; ++level) {
            const size_t row = level * (size_t)directions;
            J_I[row + (size_t)(n_mu + outgoing)] = I_out_up[level];
            J_I[row + (size_t)(n_mu - outgoing)] = I_out_down[level];
            J_Q[row + (size_t)(n_mu + outgoing)] = Q_out_up[level];
            J_Q[row + (size_t)(n_mu - outgoing)] = Q_out_down[level];
            J_U[row + (size_t)(n_mu + outgoing)] = U_out_up[level];
            J_U[row + (size_t)(n_mu - outgoing)] = U_out_down[level];
        }
    }

#undef TRANSPOSE_PLANE
    free(owned);
    return 0;
}

int rt_sos_operator_apply_scalar(const rt_atm_t *atm,
                                 int m,
                                 const rt_legendre_workspace_t *ws,
                                 const double *I_previous,
                                 double *J_I)
{
    if (!atm || !ws || !I_previous || !J_I || atm->n_mu != ws->n_mu || m < 0) {
        return -1;
    }

    const int n_mu = atm->n_mu;
    const int directions = 2 * n_mu + 1;
    const size_t levels = (size_t)(atm->n_layers + 1);
    const double beta0 = (m == 0) ? atm->beta0 : 0.0;
    const double beta2 = atm->beta2;
    const int skip_zero_columns = ocrt_sos_zero_col_skip_enabled();

    for (size_t index = 0; index < levels * (size_t)directions; ++index) {
        J_I[index] = NAN;
    }

    for (size_t level = 0; level < levels; ++level) {
        const size_t row = level * (size_t)directions;
        const double x = atm->xdel[level];
        const double y = atm->ydel[level];
        for (int outgoing = 1; outgoing <= n_mu; ++outgoing) {
            double up_sum = 0.0;
            double down_sum = 0.0;
            for (int incoming = 1; incoming <= n_mu; ++incoming) {
                const size_t up = row + (size_t)(n_mu + incoming);
                const size_t down = row + (size_t)(n_mu - incoming);
                if (!isfinite(I_previous[up]) || !isfinite(I_previous[down])) return -2;

                const double weight = atm->gb[incoming];
                if (skip_zero_columns && weight == 0.0) continue;
                const double molecular_same = beta0 +
                    beta2 * atm->xpl[incoming] * atm->xpl[outgoing];
                const double molecular_cross = beta0 +
                    beta2 * atm->xpl[incoming] * atm->xpl[-outgoing];
                const double particle_same = ws->phase_fourier_m[incoming][outgoing];
                const double particle_cross = ws->phase_fourier_m[incoming][-outgoing];
                const double same = x * particle_same + y * molecular_same;
                const double cross = x * particle_cross + y * molecular_cross;

                up_sum += weight *
                    (I_previous[up] * same + I_previous[down] * cross);
                down_sum += weight *
                    (I_previous[up] * cross + I_previous[down] * same);
            }
            J_I[row + (size_t)(n_mu + outgoing)] = up_sum;
            J_I[row + (size_t)(n_mu - outgoing)] = down_sum;
        }
    }
    return 0;
}
