/* Workspace allocation and immutable Fourier-kernel cache for OCRT. */

#include "rt_kernel.h"

#include <stdlib.h>
#include <string.h>

int rt_legendre_workspace_alloc(rt_legendre_workspace_t *ws,
                                int n_mu,
                                int l_max)
{
    if (!ws || n_mu < 1 || l_max < 0) return -1;

    memset(ws, 0, sizeof(*ws));
    ws->n_mu = n_mu;
    ws->l_max = l_max;
    ws->sos_kernel_pack_m = -1;

    const int directions = 2 * n_mu + 1;
    const size_t basis_count = (size_t)(l_max + 1) * (size_t)directions;
    const size_t kernel_count = (size_t)(n_mu + 1) * (size_t)directions;

    ws->plm_storage = (double *)calloc(basis_count, sizeof(double));
    ws->plm = (double **)calloc((size_t)(l_max + 1), sizeof(double *));
    ws->rrl_storage = (double *)calloc(basis_count, sizeof(double));
    ws->rrl = (double **)calloc((size_t)(l_max + 1), sizeof(double *));
    ws->rtl_storage = (double *)calloc(basis_count, sizeof(double));
    ws->rtl = (double **)calloc((size_t)(l_max + 1), sizeof(double *));

    ws->pfm_storage = (double *)calloc(kernel_count, sizeof(double));
    ws->phase_fourier_m = (double **)calloc((size_t)(n_mu + 1), sizeof(double *));
    ws->gr_storage = (double *)calloc(kernel_count, sizeof(double));
    ws->gr_pol = (double **)calloc((size_t)(n_mu + 1), sizeof(double *));
    ws->gt_storage = (double *)calloc(kernel_count, sizeof(double));
    ws->gt_pol = (double **)calloc((size_t)(n_mu + 1), sizeof(double *));
    ws->arr_storage = (double *)calloc(kernel_count, sizeof(double));
    ws->arr_pol = (double **)calloc((size_t)(n_mu + 1), sizeof(double *));
    ws->art_storage = (double *)calloc(kernel_count, sizeof(double));
    ws->art_pol = (double **)calloc((size_t)(n_mu + 1), sizeof(double *));
    ws->att_storage = (double *)calloc(kernel_count, sizeof(double));
    ws->att_pol = (double **)calloc((size_t)(n_mu + 1), sizeof(double *));

    ws->beam_q_to_i_direct = (double *)calloc((size_t)directions, sizeof(double));
    ws->beam_q_to_q_direct = (double *)calloc((size_t)directions, sizeof(double));
    ws->beam_q_to_u_direct = (double *)calloc((size_t)directions, sizeof(double));
    ws->beam_q_direct_valid = 0;

    if (!ws->plm_storage || !ws->plm || !ws->rrl_storage || !ws->rrl ||
        !ws->rtl_storage || !ws->rtl || !ws->pfm_storage ||
        !ws->phase_fourier_m || !ws->gr_storage || !ws->gr_pol ||
        !ws->gt_storage || !ws->gt_pol || !ws->arr_storage || !ws->arr_pol ||
        !ws->art_storage || !ws->art_pol || !ws->att_storage || !ws->att_pol ||
        !ws->beam_q_to_i_direct || !ws->beam_q_to_q_direct ||
        !ws->beam_q_to_u_direct) {
        rt_legendre_workspace_free(ws);
        return -1;
    }

    for (int l = 0; l <= l_max; ++l) {
        const size_t row = (size_t)l * (size_t)directions + (size_t)n_mu;
        ws->plm[l] = ws->plm_storage + row;
        ws->rrl[l] = ws->rrl_storage + row;
        ws->rtl[l] = ws->rtl_storage + row;
    }
    for (int j = 0; j <= n_mu; ++j) {
        const size_t row = (size_t)j * (size_t)directions + (size_t)n_mu;
        ws->phase_fourier_m[j] = ws->pfm_storage + row;
        ws->gr_pol[j] = ws->gr_storage + row;
        ws->gt_pol[j] = ws->gt_storage + row;
        ws->arr_pol[j] = ws->arr_storage + row;
        ws->art_pol[j] = ws->art_storage + row;
        ws->att_pol[j] = ws->att_storage + row;
    }
    return 0;
}

void rt_legendre_workspace_free(rt_legendre_workspace_t *ws)
{
    if (!ws) return;
    free(ws->plm_storage);
    free(ws->plm);
    free(ws->rrl_storage);
    free(ws->rrl);
    free(ws->rtl_storage);
    free(ws->rtl);
    free(ws->pfm_storage);
    free(ws->phase_fourier_m);
    free(ws->gr_storage);
    free(ws->gr_pol);
    free(ws->gt_storage);
    free(ws->gt_pol);
    free(ws->arr_storage);
    free(ws->arr_pol);
    free(ws->art_storage);
    free(ws->art_pol);
    free(ws->att_storage);
    free(ws->att_pol);
    free(ws->beam_q_to_i_direct);
    free(ws->beam_q_to_q_direct);
    free(ws->beam_q_to_u_direct);
    free(ws->sos_source_arena);
    free(ws->sos_kernel_pack);
    memset(ws, 0, sizeof(*ws));
    ws->sos_kernel_pack_m = -1;
}

#ifdef OCRT_FAST_KERNELS

typedef struct {
    int valid;
    int n_mu;
    int l_max;
    int m_count;
    unsigned long long key;
    unsigned char *have;
    double *slab;
} rt_moment_kernel_cache_t;

static rt_moment_kernel_cache_t g_moment_cache = {0};
#pragma omp threadprivate(g_moment_cache)

static unsigned long long hash_bytes(unsigned long long hash,
                                     const void *data,
                                     size_t size)
{
    const unsigned char *bytes = (const unsigned char *)data;
    if (!data) {
        hash ^= 0x9e3779b97f4a7c15ULL;
        return hash * 0x100000001b3ULL;
    }
    for (size_t i = 0; i < size; ++i) {
        hash ^= bytes[i];
        hash *= 0x100000001b3ULL;
    }
    return hash;
}

static double *cache_block(int m, int component)
{
    const size_t block = (size_t)(g_moment_cache.n_mu + 1) *
                         (size_t)(2 * g_moment_cache.n_mu + 1);
    return g_moment_cache.slab +
           ((size_t)m * 6u + (size_t)component) * block;
}

int rt_mkc_begin(const rt_atm_t *atm,
                 const rt_legendre_workspace_t *ws,
                 int m_count_hint)
{
    if (!atm || !ws || m_count_hint < 1) {
        return -1;
    }

    const int n_mu = ws->n_mu;
    const int l_max = ws->l_max;
    const int directions = 2 * n_mu + 1;
    const size_t moments = (size_t)(l_max + 1);
    unsigned long long key = 0xcbf29ce484222325ULL;

    key = hash_bytes(key, &n_mu, sizeof(n_mu));
    key = hash_bytes(key, &l_max, sizeof(l_max));
    key = hash_bytes(key, &atm->rm[-n_mu],
                     (size_t)directions * sizeof(double));
    key = hash_bytes(key, atm->betal_aer,
                     atm->betal_aer ? moments * sizeof(double) : 0);
    key = hash_bytes(key, atm->gammal_aer,
                     atm->gammal_aer ? moments * sizeof(double) : 0);
    key = hash_bytes(key, atm->alphal_aer,
                     atm->alphal_aer ? moments * sizeof(double) : 0);
    key = hash_bytes(key, atm->zetal_aer,
                     atm->zetal_aer ? moments * sizeof(double) : 0);

    if (g_moment_cache.valid && g_moment_cache.key == key &&
        g_moment_cache.n_mu == n_mu && g_moment_cache.l_max == l_max &&
        g_moment_cache.m_count >= m_count_hint) {
        return 1;
    }

    free(g_moment_cache.have);
    free(g_moment_cache.slab);
    memset(&g_moment_cache, 0, sizeof(g_moment_cache));

    const size_t block = (size_t)(n_mu + 1) * (size_t)directions;
    g_moment_cache.have = (unsigned char *)calloc((size_t)m_count_hint, 1);
    g_moment_cache.slab = (double *)malloc((size_t)m_count_hint * 6u *
                                           block * sizeof(double));
    if (!g_moment_cache.have || !g_moment_cache.slab) {
        free(g_moment_cache.have);
        free(g_moment_cache.slab);
        memset(&g_moment_cache, 0, sizeof(g_moment_cache));
        return -1;
    }

    g_moment_cache.valid = 1;
    g_moment_cache.key = key;
    g_moment_cache.n_mu = n_mu;
    g_moment_cache.l_max = l_max;
    g_moment_cache.m_count = m_count_hint;
    return 0;
}

int rt_mkc_load(int m, rt_legendre_workspace_t *ws)
{
    if (!g_moment_cache.valid || !ws || m < 0 ||
        m >= g_moment_cache.m_count || g_moment_cache.have[m] != 3 ||
        ws->n_mu != g_moment_cache.n_mu) {
        return 0;
    }

    const size_t bytes = (size_t)(g_moment_cache.n_mu + 1) *
                         (size_t)(2 * g_moment_cache.n_mu + 1) * sizeof(double);
    memcpy(ws->pfm_storage, cache_block(m, 0), bytes);
    memcpy(ws->gr_storage, cache_block(m, 1), bytes);
    memcpy(ws->gt_storage, cache_block(m, 2), bytes);
    memcpy(ws->arr_storage, cache_block(m, 3), bytes);
    memcpy(ws->art_storage, cache_block(m, 4), bytes);
    memcpy(ws->att_storage, cache_block(m, 5), bytes);
    return 1;
}

void rt_mkc_store_pfm(int m, const rt_legendre_workspace_t *ws)
{
    if (!g_moment_cache.valid || !ws || m < 0 ||
        m >= g_moment_cache.m_count || ws->n_mu != g_moment_cache.n_mu) {
        return;
    }
    const size_t bytes = (size_t)(g_moment_cache.n_mu + 1) *
                         (size_t)(2 * g_moment_cache.n_mu + 1) * sizeof(double);
    memcpy(cache_block(m, 0), ws->pfm_storage, bytes);
    g_moment_cache.have[m] |= 1;
}

void rt_mkc_store_pol(int m, const rt_legendre_workspace_t *ws)
{
    if (!g_moment_cache.valid || !ws || m < 0 ||
        m >= g_moment_cache.m_count || ws->n_mu != g_moment_cache.n_mu) {
        return;
    }
    const size_t bytes = (size_t)(g_moment_cache.n_mu + 1) *
                         (size_t)(2 * g_moment_cache.n_mu + 1) * sizeof(double);
    memcpy(cache_block(m, 1), ws->gr_storage, bytes);
    memcpy(cache_block(m, 2), ws->gt_storage, bytes);
    memcpy(cache_block(m, 3), ws->arr_storage, bytes);
    memcpy(cache_block(m, 4), ws->art_storage, bytes);
    memcpy(cache_block(m, 5), ws->att_storage, bytes);
    g_moment_cache.have[m] |= 2;
}

#endif /* OCRT_FAST_KERNELS */
