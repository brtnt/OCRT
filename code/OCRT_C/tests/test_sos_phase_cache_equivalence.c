#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rt_atm.h"
#include "rt_kernel.h"
#include "rt_solver.h"
#include "rt_sos_operator.h"


static uint64_t dbits(double x) {
    uint64_t u;
    memcpy(&u, &x, sizeof u);
    return u;
}

static int arrays_equal_bits(const double *a, const double *b, size_t n,
                             const char *name, int m, int seed) {
    for (size_t i = 0; i < n; ++i) {
        if (dbits(a[i]) != dbits(b[i])) {
            fprintf(stderr,
                    "MISMATCH %s m=%d seed=%d i=%zu a=%.17g b=%.17g bits=%016llx/%016llx\n",
                    name, m, seed, i, a[i], b[i],
                    (unsigned long long)dbits(a[i]),
                    (unsigned long long)dbits(b[i]));
            return 0;
        }
    }
    return 1;
}

int main(void) {
    enum { N_MU = 8, N_LAYERS = 13, LMAX = 12, MMAX = 6, NSEED = 4 };
    rt_atm_t atm = {0};
    rt_legendre_workspace_t ws = {0};
    double beta[LMAX + 1], gamma[LMAX + 1], alpha[LMAX + 1], zeta[LMAX + 1];

    for (int l = 0; l <= LMAX; ++l) {
        const double d = (double)l;
        beta[l]  = (l == 0) ? 1.0 : (2.0*d + 1.0) * pow(0.73, d);
        gamma[l] = (l < 2) ? 0.0 : 0.13 * pow(-0.61, d - 2.0);
        alpha[l] = (l < 2) ? 0.0 : 0.17 * pow( 0.57, d - 2.0);
        zeta[l]  = (l < 2) ? 0.0 : 0.11 * pow(-0.49, d - 2.0);
    }

    if (rt_atm_alloc(&atm, N_LAYERS, N_MU) != 0) return 2;
    if (rt_atm_build_aerosol_rayleigh(&atm,
            0.17, 0.0279, 0.43, 0.92, LMAX,
            beta, gamma, alpha, zeta,
            cos(37.0 * 3.14159265358979323846 / 180.0),
            RT_RAYLEIGH_MODEL_BODHAINE_1999, 2.0) != 0) return 3;
    if (rt_legendre_workspace_alloc(&ws, N_MU, LMAX) != 0) return 4;

    const int dirs = 2 * N_MU + 1;
    const size_t n = (size_t)(N_LAYERS + 1) * (size_t)dirs;
    double *buf = calloc((size_t)12 * n, sizeof(double));
    if (!buf) return 5;
    double *I = buf + 0*n, *Q = buf + 1*n, *U = buf + 2*n;
    double *faI = buf + 3*n, *faQ = buf + 4*n, *faU = buf + 5*n;
    double *caI = buf + 6*n, *caQ = buf + 7*n, *caU = buf + 8*n;
    double *c2I = buf + 9*n, *c2Q = buf + 10*n, *c2U = buf + 11*n;

    double *arena_ptr = NULL;
    double *pack_ptr = NULL;
    int comparisons = 0;

    for (int m = 0; m <= MMAX; ++m) {
        if (rt_legendre_compute_pol(&ws, &atm, m) != 0 ||
            rt_kernel_phase_fourier(&ws, m, atm.betal_aer) != 0 ||
            rt_kernel_phase_fourier_pol(&ws, m, atm.gammal_aer) != 0 ||
            rt_kernel_phase_fourier_aerosol_full(&ws, m,
                                                  atm.alphal_aer,
                                                  atm.zetal_aer) != 0) {
            fprintf(stderr, "kernel build failed at m=%d\n", m);
            return 6;
        }

        for (int seed = 0; seed < NSEED; ++seed) {
            for (int k = 0; k <= N_LAYERS; ++k) {
                for (int j = -N_MU; j <= N_MU; ++j) {
                    const size_t idx = (size_t)k * (size_t)dirs + (size_t)(j + N_MU);
                    const double x = (double)(1 + seed) * 0.013
                                   + (double)(k + 1) * 0.0017
                                   + (double)(j + N_MU + 1) * 0.00031;
                    I[idx] = 0.21 + sin(x * (double)(m + 1));
                    Q[idx] = 0.03 * cos((x + 0.2) * (double)(m + 2));
                    U[idx] = 0.02 * sin((x - 0.1) * (double)(m + 3));
                }
            }

            ws.sos_kernel_pack_valid = 0;
            ws.sos_kernel_pack_m = -1;
            if (rt_sos_operator_apply_vector(&atm, m, &ws, I, Q, U, faI, faQ, faU) != 0)
                return 7;
            if (rt_sos_operator_prepare(&atm, m, &ws) != 0)
                return 8;
            if (rt_sos_operator_apply_vector(&atm, m, &ws, I, Q, U, caI, caQ, caU) != 0)
                return 9;
            if (rt_sos_operator_apply_vector(&atm, m, &ws, I, Q, U, c2I, c2Q, c2U) != 0)
                return 10;

            if (!arrays_equal_bits(faI, caI, n, "I fallback/cache", m, seed) ||
                !arrays_equal_bits(faQ, caQ, n, "Q fallback/cache", m, seed) ||
                !arrays_equal_bits(faU, caU, n, "U fallback/cache", m, seed) ||
                !arrays_equal_bits(caI, c2I, n, "I cache/reuse", m, seed) ||
                !arrays_equal_bits(caQ, c2Q, n, "Q cache/reuse", m, seed) ||
                !arrays_equal_bits(caU, c2U, n, "U cache/reuse", m, seed))
                return 11;
            comparisons += 6;

            if (!arena_ptr) arena_ptr = ws.sos_source_arena;
            if (!pack_ptr) pack_ptr = ws.sos_kernel_pack;
            if (arena_ptr != ws.sos_source_arena || pack_ptr != ws.sos_kernel_pack) {
                fprintf(stderr, "cache allocation moved after initial reserve\n");
                return 12;
            }
        }
    }

    printf("PASS direct-cache-equivalence modes=%d seeds=%d array_compares=%d values_per_array=%zu\n",
           MMAX + 1, NSEED, comparisons, n);
    free(buf);
    rt_legendre_workspace_free(&ws);
    rt_atm_free(&atm);
    return 0;
}
