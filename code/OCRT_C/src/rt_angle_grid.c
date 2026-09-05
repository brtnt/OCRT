#include "rt_angle_grid.h"

#include <math.h>
#include <stddef.h>
#include <string.h>

#include "rt_quadrature.h"

static int lower_bound_mu(const double *values, int count, double target)
{
    int lo = 0;
    int hi = count;
    while (lo < hi) {
        const int mid = lo + (hi - lo) / 2;
        if (values[mid] < target) lo = mid + 1;
        else                      hi = mid;
    }
    return lo;
}

int rt_uangles_init(rt_uangles_t *table, int n_gauss)
{
    if (table == NULL || n_gauss < 2 || n_gauss > RT_UANG_MAX) return -1;

    memset(table, 0, sizeof(*table));
    if (rt_quadrature_gauss_legendre_pos(n_gauss,
                                         table->mu,
                                         table->w) != 0) {
        return -1;
    }

    table->n_gauss = n_gauss;
    table->n_total = n_gauss;
    return 0;
}

int rt_uangles_add(rt_uangles_t *table, double mu_new, int *index_out)
{
    if (table == NULL || !isfinite(mu_new) || mu_new < 0.0 || mu_new > 1.0) {
        return -1;
    }
    if (table->n_total < 0 || table->n_total >= RT_UANG_MAX) return -1;

    const int pos = lower_bound_mu(table->mu, table->n_total, mu_new);

    /* Only the immediate neighbours of the insertion point can satisfy the
     * tolerance because the table is sorted.  Prefer the lower-index match,
     * matching the stable behaviour expected by existing index maps. */
    if (pos > 0 && fabs(table->mu[pos - 1] - mu_new) < RT_ANGLE_DEDUP_TOL) {
        if (index_out != NULL) *index_out = pos - 1;
        return 0;
    }
    if (pos < table->n_total &&
        fabs(table->mu[pos] - mu_new) < RT_ANGLE_DEDUP_TOL) {
        if (index_out != NULL) *index_out = pos;
        return 0;
    }

    const size_t tail = (size_t)(table->n_total - pos);
    if (tail > 0) {
        memmove(&table->mu[pos + 1], &table->mu[pos], tail * sizeof(table->mu[0]));
        memmove(&table->w[pos + 1], &table->w[pos], tail * sizeof(table->w[0]));
        memmove(&table->is_special[pos + 1], &table->is_special[pos],
                tail * sizeof(table->is_special[0]));
    }

    table->mu[pos] = mu_new;
    table->w[pos] = 0.0;
    table->is_special[pos] = 1;
    table->n_total += 1;
    if (index_out != NULL) *index_out = pos;
    return 0;
}
