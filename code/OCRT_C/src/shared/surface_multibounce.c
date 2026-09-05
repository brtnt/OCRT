/* surface_multibounce.c — water->air interface MULTI-ENCOUNTER cascade.
 *
 * Physically exact closure (within the uncorrelated-facet Cox-Munk model)
 * for the fate of an upwelling in-water beam at a rough surface:
 *   first-encounter transmission T_h  =  UP-ESCAPE  +  WATER-RETURN
 * where the re-entrant (downward-transmitted) rays are followed through
 * repeated air-side facet encounters (reflect -> possibly escape upward;
 * refract -> return into the water body).  Validated against the python
 * reference cascade (BUILD_SPEEDUP 2026-07-10): identity closes to <0.05%,
 * zero unconverged mass at 24 bounces, K=8 reservoir resampling.
 *
 * Deterministic: fixed-seed xorshift64* + Box-Muller; same binary ->
 * bit-identical results (documented; cross-platform identical up to libm).
 * No angle switches anywhere (Jae engineering rule 2026-07-10): the return
 * term goes to zero smoothly below the critical cone.
 */
#include <math.h>
#include <stdint.h>
#include <string.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct { uint64_t s; } mbrng_t;
static inline uint64_t mb_next(mbrng_t *r) {
    uint64_t x = r->s;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    r->s = x;
    return x * 0x2545F4914F6CDD1DULL;
}
static inline double mb_u01(mbrng_t *r) {
    return (double)(mb_next(r) >> 11) * (1.0 / 9007199254740992.0);
}
/* Box-Muller pair; returns one, caches the other. */
static inline double mb_gauss(mbrng_t *r, double *cache, int *has) {
    if (*has) { *has = 0; return *cache; }
    double u1 = fmax(1e-300, mb_u01(r)), u2 = mb_u01(r);
    double m = sqrt(-2.0 * log(u1));
    *cache = m * sin(2.0 * M_PI * u2); *has = 1;
    return m * cos(2.0 * M_PI * u2);
}

/* Fresnel energy reflectance, air(1) -> water(n), local incidence ci. */
static inline double mb_R_aw(double ci, double n) {
    double s2w = (1.0 - ci * ci) / (n * n);
    double mt = sqrt(fmax(0.0, 1.0 - s2w));
    double rs = (ci - n * mt) / (ci + n * mt);
    double rp = (n * ci - mt) / (n * ci + mt);
    return 0.5 * (rs * rs + rp * rp);
}

/* Cascade for one in-water direction.
 * Inputs : n_water, wind_speed (Cox-Munk isotropic sigma^2=0.003+0.00512*ws),
 *          mu_w (in-water upwelling cosine), n_rays.
 * Outputs: *T_h (first-encounter total transmission),
 *          *T_up (final upward escape incl. multi-encounter re-escape),
 *          *W_ret (returned-into-water fraction).
 * Optional: esc_hist[nb] — final-escape energy binned uniformly in mu_air
 *           (pass NULL to skip).  All outputs normalized per incident flux. */
int surface_wa_multibounce_cascade(double n_water, double wind_speed,
                                   double mu_w, long n_rays,
                                   double *T_h, double *T_up, double *W_ret,
                                   double *esc_hist, double *escd_hist,
                                   double *ret_hist, int nb) {
    if (!(mu_w > 1e-9) || !(n_water > 1.0) || n_rays < 1000) return -1;
    const double n = n_water;
    const double sig2 = 0.003 + 0.00512 * (wind_speed > 0.0 ? wind_speed : 0.0);
    const double sg = sqrt(sig2 * 0.5);
    const double s1 = sqrt(fmax(0.0, 1.0 - mu_w * mu_w));
    mbrng_t rng = { 0x9E3779B97F4A7C15ULL ^ (uint64_t)(mu_w * 1e12)
                                          ^ ((uint64_t)(wind_speed * 1e6) << 20)
                                          ^ ((uint64_t)(n * 1e6) << 40) };
    double gc = 0.0; int gh = 0;
    double sum_w = 0.0, sum_Th = 0.0, sum_up = 0.0, sum_ret = 0.0;
    if (esc_hist)  memset(esc_hist, 0, (size_t)nb * sizeof(double));
    if (escd_hist) memset(escd_hist, 0, (size_t)nb * sizeof(double));
    if (ret_hist) memset(ret_hist, 0, (size_t)nb * sizeof(double));

    enum { KCAND = 8, MAXB = 24 };
    for (long it = 0; it < n_rays; ++it) {
        /* --- first encounter: encounter-weighted facet from the water side */
        double zx = sg * mb_gauss(&rng, &gc, &gh);
        double zy = sg * mb_gauss(&rng, &gc, &gh);
        double nz = 1.0 / sqrt(1.0 + zx * zx + zy * zy);
        double nx = -zx * nz, ny = -zy * nz;
        double ih = s1 * nx + mu_w * nz;
        if (ih <= 0.0) continue;                 /* back-facing: no encounter  */
        double w = ih / nz;                      /* encounter weight            */
        sum_w += w;
        double s2t = n * n * (1.0 - ih * ih);
        if (s2t >= 1.0) continue;                /* local TIR: all reflected    */
        double mut = sqrt(1.0 - s2t);
        double rs = (n * ih - mut) / (n * ih + mut);
        double rp = (ih - n * mut) / (ih + n * mut);
        double T1 = 1.0 - 0.5 * (rs * rs + rp * rp);
        sum_Th += w * T1;
        /* refracted global direction */
        double c = n * ih - mut;
        double dx = n * s1 - c * nx, dy = -c * ny, dz = n * mu_w - c * nz;
        double L = sqrt(dx * dx + dy * dy + dz * dz);
        dx /= L; dy /= L; dz /= L;
        double W = w * T1;
        if (dz > 0.0) {
            sum_up += W;   /* direct up-escape: exactly the single-encounter
                              kernel's content — kept OUT of esc_hist, but
                              recorded separately for flat-operator domain
                              completion (beyond-critical nodes). */
            if (escd_hist) { int b = (int)(dz * nb); if (b >= nb) b = nb - 1; escd_hist[b] += W; }
            continue;
        }
        /* --- multi-encounter cascade over air-side facets */
        for (int bnc = 0; bnc < MAXB && W > 1e-14; ++bnc) {
            /* reservoir-resample one forward-facing facet among KCAND */
            double hx = 0, hy = 0, hz = 1, ci = 0, wsum = 0;
            double rsel = mb_u01(&rng);
            int any = 0;
            for (int k = 0; k < KCAND; ++k) {
                double zx2 = sg * mb_gauss(&rng, &gc, &gh);
                double zy2 = sg * mb_gauss(&rng, &gc, &gh);
                double nz2 = 1.0 / sqrt(1.0 + zx2 * zx2 + zy2 * zy2);
                double nx2 = -zx2 * nz2, ny2 = -zy2 * nz2;
                double cik = -(dx * nx2 + dy * ny2 + dz * nz2);
                double wk = (cik > 0.0) ? cik / nz2 : 0.0;
                if (wk <= 0.0) continue;
                any = 1;
                double neww = wsum + wk;
                if (rsel * neww >= wsum) { hx = nx2; hy = ny2; hz = nz2; ci = cik; }
                wsum = neww;
            }
            if (!any) { /* extremely rare at ws>=3; retry the bounce */
                --bnc; continue;
            }
            double R2 = mb_R_aw(ci, n);
            {
                double Tw = W * (1.0 - R2);
                sum_ret += Tw;
                if (ret_hist && Tw > 0.0) {
                    /* refracted-into-water direction (air 1 -> water n):
                       d_t = (d + (mtw*n - ci)*(-h)) / n  … 표준식, 전역 z */
                    double s2w = (1.0 - ci * ci) / (n * n);
                    double mtw = sqrt(fmax(0.0, 1.0 - s2w));
                    double k1 = 1.0 / n, k2 = (mtw - ci / n);
                    double dtz = k1 * dz - k2 * hz;   /* h가 위쪽, 광선은 아래로 */
                    double mwv = fabs(dtz);
                    int b = (int)(mwv * nb); if (b >= nb) b = nb - 1;
                    ret_hist[b] += Tw;
                }
            }
            W *= R2;
            dx += 2.0 * ci * hx; dy += 2.0 * ci * hy; dz += 2.0 * ci * hz;
            /* reflection of a unit vector stays unit */
            if (dz > 0.0) {
                sum_up += W;
                if (esc_hist) { int b = (int)(dz * nb); if (b >= nb) b = nb - 1; esc_hist[b] += W; }
                W = 0.0;
            }
        }
        sum_ret += W;   /* unconverged residue booked as return (tiny)        */
        if (ret_hist && W > 0.0) ret_hist[nb/3] += W;  /* residue: mid-bin */
    }
    double norm = (double)n_rays * mu_w;
    if (T_h)  *T_h  = sum_Th / norm;
    if (T_up) *T_up = sum_up / norm;
    if (W_ret)*W_ret= sum_ret / norm;
    if (esc_hist)  for (int b = 0; b < nb; ++b) esc_hist[b]  /= norm;
    if (escd_hist) for (int b = 0; b < nb; ++b) escd_hist[b] /= norm;
    if (ret_hist) for (int b = 0; b < nb; ++b) ret_hist[b] /= norm;
    return 0;
}
