#ifndef OCRT_RT_ANGLE_GRID_H
#define OCRT_RT_ANGLE_GRID_H

/* Shared positive-cosine table used by atmosphere and water solvers.
 *
 * The table begins with an ascending positive Gauss-Legendre rule.  Output
 * or boundary directions may then be inserted as zero-weight nodes.  A new
 * value within RT_ANGLE_DEDUP_TOL of an existing node reuses that node.
 */

#define RT_UANG_MAX 640
#define RT_ANGLE_DEDUP_TOL 1.0e-12

typedef struct {
    int n_total;
    int n_gauss;
    double mu[RT_UANG_MAX];
    double w[RT_UANG_MAX];
    unsigned char is_special[RT_UANG_MAX];
} rt_uangles_t;

int rt_uangles_init(rt_uangles_t *table, int n_gauss);
int rt_uangles_add(rt_uangles_t *table, double mu_new, int *index_out);

#endif /* OCRT_RT_ANGLE_GRID_H */
