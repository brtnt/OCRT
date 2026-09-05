/* rt_water_iop.c — Pure water IOP module implementation
 *
 * See rt_water_iop.h for detailed source/method documentation.
 * 2026-05-22 KST (Phase A: pure water only; seawater constituents in Phase B).
 */
#include "rt_water_iop.h"
#include "rt_spectral_contract.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ===== LUT load / free ===== */

int rt_water_iop_lut_load(const char* filepath, rt_water_iop_lut_t* lut) {
    if (!filepath || !lut) return -1;
    lut->n = 0;
    lut->lambda_nm = NULL;
    lut->aw_m_inv  = NULL;
    lut->bw_m_inv  = NULL;

    FILE* f = fopen(filepath, "r");
    if (!f) {
        fprintf(stderr, "rt_water_iop_lut_load: cannot open '%s'\n", filepath);
        return -2;
    }

    /* Two-pass: count data rows, then allocate, then re-parse */
    char line[1024];
    int in_data = 0, n_rows = 0;
    while (fgets(line, sizeof line, f)) {
        if (!in_data) {
            if (strncmp(line, "/end_header", 11) == 0) in_data = 1;
            continue;
        }
        /* Data row: at least 3 tokens */
        double a, b, c;
        if (sscanf(line, "%lf %lf %lf", &a, &b, &c) == 3) n_rows++;
    }
    if (n_rows < 2) {
        fclose(f);
        fprintf(stderr, "rt_water_iop_lut_load: too few data rows in '%s' (%d)\n",
                filepath, n_rows);
        return -3;
    }

    lut->n         = n_rows;
    lut->lambda_nm = (double*)malloc(n_rows * sizeof(double));
    lut->aw_m_inv  = (double*)malloc(n_rows * sizeof(double));
    lut->bw_m_inv  = (double*)malloc(n_rows * sizeof(double));
    if (!lut->lambda_nm || !lut->aw_m_inv || !lut->bw_m_inv) {
        rt_water_iop_lut_free(lut);
        fclose(f);
        return -4;
    }

    rewind(f);
    in_data = 0;
    int k = 0;
    while (fgets(line, sizeof line, f)) {
        if (!in_data) {
            if (strncmp(line, "/end_header", 11) == 0) in_data = 1;
            continue;
        }
        double a, b, c;
        if (sscanf(line, "%lf %lf %lf", &a, &b, &c) == 3) {
            lut->lambda_nm[k] = a;
            lut->aw_m_inv[k]  = b;
            lut->bw_m_inv[k]  = c;
            k++;
        }
    }
    fclose(f);

    /* Verify monotonic ascending in wavelength */
    for (int i = 1; i < lut->n; ++i) {
        if (lut->lambda_nm[i] <= lut->lambda_nm[i-1]) {
            fprintf(stderr, "rt_water_iop_lut_load: non-monotonic λ at i=%d "
                    "(%.4f <= %.4f)\n", i, lut->lambda_nm[i], lut->lambda_nm[i-1]);
            rt_water_iop_lut_free(lut);
            return -5;
        }
    }

    if (!rt_spectral_grid_covers(lut->lambda_nm[0],
                                 lut->lambda_nm[lut->n - 1])) {
        fprintf(stderr,
                "rt_water_iop_lut_load: table '%s' does not cover %.0f-%.0f nm\n",
                filepath, RT_SPECTRAL_MIN_NM, RT_SPECTRAL_MAX_NM);
        rt_water_iop_lut_free(lut);
        return -6;
    }
    fprintf(stderr, "rt_water_iop_lut: loaded %d rows from '%s' (λ %.2f..%.2f nm)\n",
            lut->n, filepath, lut->lambda_nm[0], lut->lambda_nm[lut->n - 1]);
    return 0;
}

void rt_water_iop_lut_free(rt_water_iop_lut_t* lut) {
    if (!lut) return;
    free(lut->lambda_nm); lut->lambda_nm = NULL;
    free(lut->aw_m_inv);  lut->aw_m_inv  = NULL;
    free(lut->bw_m_inv);  lut->bw_m_inv  = NULL;
    lut->n = 0;
}

/* ===== psi_T LUT loader (2-column file: lambda, psi_T) ===== */

int rt_water_iop_psi_T_load(const char* filepath, rt_water_iop_psi_T_lut_t* lut) {
    if (!filepath || !lut) return -1;
    lut->n = 0; lut->lambda_nm = NULL; lut->psi_T = NULL; lut->T_ref = 20.0;

    FILE* f = fopen(filepath, "r");
    if (!f) {
        fprintf(stderr, "rt_water_iop_psi_T_load: cannot open '%s'\n", filepath);
        return -2;
    }

    char line[1024];
    int in_data = 0, n_rows = 0;
    while (fgets(line, sizeof line, f)) {
        if (!in_data) {
            if (strncmp(line, "/end_header", 11) == 0) in_data = 1;
            continue;
        }
        double a, b;
        if (sscanf(line, "%lf %lf", &a, &b) == 2) n_rows++;
    }
    if (n_rows < 2) {
        fclose(f);
        fprintf(stderr, "rt_water_iop_psi_T_load: too few data rows (%d)\n", n_rows);
        return -3;
    }

    lut->n         = n_rows;
    lut->lambda_nm = (double*)malloc(n_rows * sizeof(double));
    lut->psi_T     = (double*)malloc(n_rows * sizeof(double));
    if (!lut->lambda_nm || !lut->psi_T) {
        rt_water_iop_psi_T_free(lut); fclose(f); return -4;
    }

    rewind(f);
    in_data = 0;
    int k = 0;
    while (fgets(line, sizeof line, f)) {
        if (!in_data) {
            if (strncmp(line, "/end_header", 11) == 0) in_data = 1;
            continue;
        }
        double a, b;
        if (sscanf(line, "%lf %lf", &a, &b) == 2) {
            lut->lambda_nm[k] = a;
            lut->psi_T[k]     = b;
            k++;
        }
    }
    fclose(f);

    for (int i = 1; i < lut->n; ++i) {
        if (lut->lambda_nm[i] <= lut->lambda_nm[i-1]) {
            fprintf(stderr, "rt_water_iop_psi_T_load: non-monotonic λ at i=%d\n", i);
            rt_water_iop_psi_T_free(lut);
            return -5;
        }
    }
    if (!rt_spectral_grid_covers(lut->lambda_nm[0],
                                 lut->lambda_nm[lut->n - 1])) {
        fprintf(stderr,
                "rt_water_iop_psi_T_load: table '%s' does not cover %.0f-%.0f nm\n",
                filepath, RT_SPECTRAL_MIN_NM, RT_SPECTRAL_MAX_NM);
        rt_water_iop_psi_T_free(lut);
        return -6;
    }
    fprintf(stderr, "rt_water_iop_psi_T: loaded %d rows from '%s' (λ %.1f..%.1f nm)\n",
            lut->n, filepath, lut->lambda_nm[0], lut->lambda_nm[lut->n - 1]);
    return 0;
}

void rt_water_iop_psi_T_free(rt_water_iop_psi_T_lut_t* lut) {
    if (!lut) return;
    free(lut->lambda_nm); lut->lambda_nm = NULL;
    free(lut->psi_T);     lut->psi_T     = NULL;
    lut->n = 0;
}

/* ===== Wavelength interpolation ===== */

/* Internal: in-range linear interpolation.  Endpoint values are valid data;
 * a query outside the table returns NAN and is marked extrapolated. */
static double interp_lin_nearest_ext(const double* x, const double* y, int n,
                                     double xq, int* extrapolated) {
    if (extrapolated) *extrapolated = 0;
    if (n <= 0 || !x || !y || !isfinite(xq)) return NAN;
    if (xq < x[0])     { if (extrapolated) *extrapolated = 1; return NAN; }
    if (xq > x[n - 1]) { if (extrapolated) *extrapolated = 1; return NAN; }
    if (xq == x[0]) return y[0];
    if (xq == x[n - 1]) return y[n - 1];

    /* Binary search */
    int lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (x[mid] <= xq) lo = mid; else hi = mid;
    }
    double t = (xq - x[lo]) / (x[hi] - x[lo]);
    return (1.0 - t) * y[lo] + t * y[hi];
}

double rt_water_iop_aw(const rt_water_iop_lut_t* lut, double lambda_nm,
                       int* extrapolated) {
    return interp_lin_nearest_ext(lut->lambda_nm, lut->aw_m_inv, lut->n,
                                  lambda_nm, extrapolated);
}

double rt_water_iop_bw(const rt_water_iop_lut_t* lut, double lambda_nm,
                       int* extrapolated) {
    return interp_lin_nearest_ext(lut->lambda_nm, lut->bw_m_inv, lut->n,
                                  lambda_nm, extrapolated);
}

/* Temperature-corrected absorption.
 * Linear interpolation in λ for both a_w(20) LUT and ψ_T LUT.
 * Tables must explicitly cover the OCRT support interval; no endpoint
 * extrapolation is allowed. */
double rt_water_iop_aw_T(const rt_water_iop_lut_t*       aw_lut,
                         const rt_water_iop_psi_T_lut_t* psi_T_lut,
                         double lambda_nm,
                         double T_celsius,
                         int*   extrapolated) {
    int ext_a = 0, ext_p = 0;
    double aw20 = interp_lin_nearest_ext(aw_lut->lambda_nm, aw_lut->aw_m_inv,
                                         aw_lut->n, lambda_nm, &ext_a);
    if (extrapolated) *extrapolated = ext_a;

    if (!psi_T_lut || psi_T_lut->n == 0) {
        /* No ψ_T data → no T correction */
        return aw20;
    }

    double T_ref = psi_T_lut->T_ref;
    double psi   = interp_lin_nearest_ext(psi_T_lut->lambda_nm, psi_T_lut->psi_T,
                                          psi_T_lut->n, lambda_nm, &ext_p);
    if (extrapolated && ext_p) *extrapolated = 1;

    return aw20 + psi * (T_celsius - T_ref);
}

/* ===== Refractive index ===== */

/* Real refractive index of seawater n_w(lambda, T, S).
 *
 * Quan & Fry (1995), Appl. Opt. 34(18), 3477-3480, empirical fit (their
 * Eq. (3) family; full T,S-dependent form evaluated here at FIXED
 * T = 20 degC, S = 38.4 g/kg):
 *
 *   n_w = n0 + (n1 + n2*T + n3*T^2)*S + n4*T^2
 *       + (n5 + n6*S + n7*T)/lambda
 *       + n8/lambda^2 + n9/lambda^3
 *
 * Symbols / units:
 *   lambda = wavelength in nm (vacuum), fit range ~200-1100 nm
 *   T      = temperature [degC], fit range 0-30 degC
 *   S      = salinity [g/kg = PSU], fit range 0-35 (S = 38.4 is a mild
 *            extrapolation of the fit; value chosen to match the frozen
 *            NASA water-IOP file conditions — provenance of that pairing
 *            is the upstream data file header, not re-derived here)
 *   n0..n9 = Quan & Fry global fit coefficients (dimensionless / nm powers)
 *
 * NOTE: T and S are compile-time constants here (Phase A pure-water scope);
 * a T,S-parameterized API can expose the same polynomial unchanged. */
double rt_water_iop_n_real(double lambda_nm) {
    /* Quan & Fry coefficients (full T,S form): */
    const double n0 = 1.31405;
    const double n1 = 1.779e-4;
    const double n2 = -1.05e-6;
    const double n3 = 1.6e-8;
    const double n4 = -2.02e-6;
    const double n5 = 15.868;
    const double n6 = 0.01155;
    const double n7 = -0.00423;
    const double n8 = -4382.0;
    const double n9 = 1.1455e6;

    const double T = 20.0;
    const double S = 38.4;

    double lambda = lambda_nm;  /* nm */
    double n = n0
             + (n1 + n2*T + n3*T*T) * S
             + n4*T*T
             + (n5 + n6*S + n7*T) / lambda
             + n8 / (lambda*lambda)
             + n9 / (lambda*lambda*lambda);
    return n;
}

/* Depolarization ratio of pure seawater delta_w(lambda)  [dimensionless].
 *
 * Constant approximation delta_w = 0.039 (Zhang & Hu / Zhang et al. 2009
 * family: weak lambda-dependence over 400-700 nm; Farinato & Rowell 1976
 * measured 0.051 @200nm, 0.039 @500nm, 0.038 @700nm — i.e. <3% variation
 * across the OCRT production bands 412-865 nm).
 *
 * CONSISTENCY: this same delta = 0.039 is hard-wired as DELTA_RAY in the
 * OCRT<->OSOAA consistency harness and patched into OSOAA (CTE_MDF_SEA =
 * 0.039); changing it here without changing both breaks the parity setup.
 * It enters the water Rayleigh phase via:
 *   P(theta) proportional to 1 + c2*cos^2(theta), c2 = (1-delta)/(1+delta)
 * (see rayleigh_p11 in the harness and the Mueller builder below). */
double rt_water_iop_depol(double lambda_nm) {
    (void)lambda_nm;  /* unused — constant approximation */
    return 0.039;
}

/* ===== Mueller scattering matrix (pure water, Rayleigh-like) =====
 *
 * Depolarized-Rayleigh phase matrix (Mishchenko et al. 2002 convention),
 * 3x3 {I,Q,U} block, angle theta = scattering angle.
 * Normalization: (1/4pi) * Integral P11(theta) dOmega = 1.
 * (Check: <cos^2> over the sphere = 1/3 =>
 *  norm*[(1-d)*(4/3) + 2d] = (4+2d)/(2(2+d)) = 1.  Exact.)
 *
 * With norm = 3 / (2*(2+delta)):
 *
 *   P11(theta) =  norm * [ (1-delta)*(1 + cos^2 theta) + 2*delta ]
 *   P12(theta) = -norm * (1-delta) * sin^2 theta          (I<->Q coupling)
 *   P33(theta) =  2*norm * (1-delta) * cos theta          (U<->U rotation)
 *   P22(theta) =  P11(theta)                              ** APPROXIMATION **
 *
 * P22 approximation note (documented deviation, do not "fix" silently):
 * in the exact depolarized-Rayleigh decomposition the isotropic
 * depolarization term (+2*delta inside the bracket) contributes to a1=P11
 * ONLY; the exact a2 is
 *   P22_exact(theta) = norm * (1-delta) * (1 + cos^2 theta)
 *                    = P11(theta) - 2*norm*delta,
 * a constant offset of magnitude 2*norm*delta ~ 0.057 (delta = 0.039).
 * Setting P22 = P11 therefore overestimates the Q->Q channel by that
 * additive constant.  Impact is limited to second-order polarized water
 * scattering (Q entering as source of Q); intensity and single-scatter
 * polarization are unaffected.  Full 4x4 (V, P44) is out of scope.
 *
 * Output layout: row-major 3x3, Z_out[r*3+c], symmetric with zeros in the
 * U row/column except P33 (Rayleigh has no I<->U / Q<->U coupling). */
void rt_water_iop_mueller(double theta_deg, double lambda_nm, double Z_out[9]) {
    double delta = rt_water_iop_depol(lambda_nm);
    double th    = theta_deg * M_PI / 180.0;
    double mu    = cos(th);
    double mu2   = mu * mu;
    double sin2  = 1.0 - mu2;

    double norm  = 3.0 / (2.0 * (2.0 + delta));

    double P11 = norm * ((1.0 - delta) * (1.0 + mu2) + 2.0 * delta);
    double P12 = -norm * (1.0 - delta) * sin2;
    double P33 = 2.0 * norm * (1.0 - delta) * mu;
    double P22 = P11;  /* 3-vector approximation; full 4x4 differs */

    /* Row-major 3x3: */
    Z_out[0] = P11; Z_out[1] = P12; Z_out[2] = 0.0;
    Z_out[3] = P12; Z_out[4] = P22; Z_out[5] = 0.0;
    Z_out[6] = 0.0; Z_out[7] = 0.0; Z_out[8] = P33;
}

/* ============================================================================
 * Component-based IOP evaluators (B.5+, 2026-05-23)
 * ----------------------------------------------------------------------------
 * 각 component 가 rt_iop_t (a, b, bb) 를 동일한 약속으로 채운다. RT 솔버는
 * rt_iop_total() 결과만 사용. CDOM 처럼 산란이 없는 component 는 b=bb=0.
 * ========================================================================= */

int rt_iop_pure_water_eval(const rt_water_iop_lut_t*       aw_lut,
                           const rt_water_iop_psi_T_lut_t* psi_T_lut,
                           double lambda_nm,
                           double T_celsius,
                           rt_iop_t*                       out,
                           int*                            extrapolated)
{
    if (!out || !aw_lut) return -1;
    int ex_local = 0;
    double a_w  = rt_water_iop_aw_T(aw_lut, psi_T_lut, lambda_nm, T_celsius, &ex_local);
    double b_w  = rt_water_iop_bw(aw_lut, lambda_nm, &ex_local);
    if (!(a_w > 0.0) || !(b_w > 0.0)) return -2;
    out->a  = a_w;
    out->b  = b_w;
    /* bb/b = 0.5 EXACTLY for any phase function even in cos(theta):
     * bb/b = (1/2) int_{90..180} P sin dtheta and the depolarized-Rayleigh
     * P(theta) above is symmetric about 90 deg, so forward and backward
     * halves integrate equally.  Not an approximation. */
    out->bb = 0.5 * b_w;
    if (extrapolated) *extrapolated = ex_local;
    return 0;
}

/* CDOM (colored dissolved organic matter) absorption  [m^-1].
 *
 *   a_cdom(lambda) = a_ref * exp( -S * (lambda - lambda_ref) )
 *
 * Symbols / units:
 *   a_ref      = params->a440, absorption at the reference wavelength [m^-1]
 *   lambda_ref = params->lambda_ref_nm [nm] (conventionally 440 nm)
 *   S          = params->S_nm_inv, spectral slope [nm^-1]
 *                (option default 0.014 nm^-1; typical ocean range 0.01-0.02,
 *                 e.g. Bricaud, Morel & Prieur 1981)
 * CDOM is purely absorbing: b = bb = 0 by construction.
 * params == NULL or a440 <= 0 => component off (all-zero IOP, rc 0). */
int rt_iop_cdom_eval(const rt_iop_cdom_params_t* params,
                     double lambda_nm,
                     rt_iop_t*                   out)
{
    if (!out) return -1;
    out->a = 0.0;
    out->b = 0.0;
    out->bb = 0.0;
    if (!params) return 0;                    /* params NULL → CDOM off */
    if (!(params->a440 > 0.0)) return 0;      /* explicit off */
    if (!(params->lambda_ref_nm > 0.0)) return -1;
    double a = params->a440 *
               exp(-params->S_nm_inv * (lambda_nm - params->lambda_ref_nm));
    if (a < 0.0) a = 0.0;
    out->a = a;
    return 0;
}


/* ---------------------------------------------------------------------------
 * CCRR constituent IOP adapter: chlorophyll absorption table
 *
 * The bundled canonical table is not a Bricaud (1998) A(lambda), E(lambda)
 * data set.  It is the user-supplied Morel (1988) normalized A_chl(lambda)
 * spectrum converted to the Morel-Maritorena / HydroLight classic Case-1
 * closure
 *
 *   a_p(lambda) = 0.06 A_chl(lambda) Chl^0.65
 *
 * and represented through the legacy five-column OCRT interface as
 * Aphi(lambda)=0.06 A_chl(lambda), Ephi(lambda)=0.65.  The historical
 * aph_bricaud_1998.txt filename is retained only as a compatibility fallback.
 * --------------------------------------------------------------------------- */
typedef struct { double wl, Aphi, Ephi; } ccrr_chl_abs_row_t;
static ccrr_chl_abs_row_t *g_ccrr_chl_abs = NULL;
static int g_ccrr_chl_abs_n = 0;

static int ccrr_load_chl_absorption(void)
{
    if (g_ccrr_chl_abs_n > 0) return 0;
    const char *paths[] = {
        "inputs/water_iop/aph_ccrr_morel1988_mm01.txt",
        "../inputs/water_iop/aph_ccrr_morel1988_mm01.txt",
        /* Legacy filename used by OCRT <= v1.18. */
        "inputs/water_iop/aph_bricaud_1998.txt",
        "../inputs/water_iop/aph_bricaud_1998.txt",
        NULL
    };
    FILE *f = NULL;
    const char *loaded_path = NULL;
    for (int ip = 0; paths[ip]; ++ip) {
        f = fopen(paths[ip], "r");
        if (f) { loaded_path = paths[ip]; break; }
    }
    if (!f) {
        fprintf(stderr,
                "rt_iop_ccrr: cannot open aph_ccrr_morel1988_mm01.txt "
                "(or legacy aph_bricaud_1998.txt)\n");
        return -1;
    }
    int cap = 256, n = 0;
    ccrr_chl_abs_row_t *rows = (ccrr_chl_abs_row_t*)calloc((size_t)cap, sizeof(*rows));
    if (!rows) { fclose(f); return -2; }
    char line[512];
    while (fgets(line, sizeof line, f)) {
        char *q = line;
        while (*q == ' ' || *q == '\t') ++q;
        if (*q == '\0' || *q == '\n' || *q == '\r' || *q == '!' || *q == '/') continue;
        double wl, Ap, Ep, Aphi, Ephi;
        if (sscanf(q, " %lf , %lf , %lf , %lf , %lf", &wl, &Ap, &Ep, &Aphi, &Ephi) == 5 ||
            sscanf(q, " %lf %lf %lf %lf %lf", &wl, &Ap, &Ep, &Aphi, &Ephi) == 5) {
            (void)Ap;
            (void)Ep;
            if (!(wl > 0.0) || !(Aphi >= 0.0) || !isfinite(Ephi)) {
                fprintf(stderr, "rt_iop_ccrr: invalid Chl absorption row in '%s'\n", loaded_path);
                free(rows); fclose(f); return -3;
            }
            if (n > 0 && !(wl > rows[n-1].wl)) {
                fprintf(stderr, "rt_iop_ccrr: non-monotonic wavelength in '%s'\n", loaded_path);
                free(rows); fclose(f); return -3;
            }
            if (n >= cap) {
                cap *= 2;
                ccrr_chl_abs_row_t *tmp = (ccrr_chl_abs_row_t*)realloc(rows, (size_t)cap * sizeof(*rows));
                if (!tmp) { free(rows); fclose(f); return -2; }
                rows = tmp;
            }
            rows[n].wl = wl; rows[n].Aphi = Aphi; rows[n].Ephi = Ephi; ++n;
        }
    }
    fclose(f);
    if (n < 2) { free(rows); return -3; }
    g_ccrr_chl_abs = rows;
    g_ccrr_chl_abs_n = n;
    return 0;
}

static double ccrr_interp_chl_absorption(double lambda_nm, int which)
{
    if (ccrr_load_chl_absorption() != 0) return NAN;
    if (lambda_nm <= g_ccrr_chl_abs[0].wl)
        return which == 0 ? g_ccrr_chl_abs[0].Aphi : g_ccrr_chl_abs[0].Ephi;
    if (lambda_nm >= g_ccrr_chl_abs[g_ccrr_chl_abs_n-1].wl)
        return which == 0 ? g_ccrr_chl_abs[g_ccrr_chl_abs_n-1].Aphi : g_ccrr_chl_abs[g_ccrr_chl_abs_n-1].Ephi;
    int lo = 0, hi = g_ccrr_chl_abs_n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (g_ccrr_chl_abs[mid].wl <= lambda_nm) lo = mid; else hi = mid;
    }
    double t = (lambda_nm - g_ccrr_chl_abs[lo].wl) /
               (g_ccrr_chl_abs[hi].wl - g_ccrr_chl_abs[lo].wl);
    double y0 = which == 0 ? g_ccrr_chl_abs[lo].Aphi : g_ccrr_chl_abs[lo].Ephi;
    double y1 = which == 0 ? g_ccrr_chl_abs[hi].Aphi : g_ccrr_chl_abs[hi].Ephi;
    return (1.0 - t) * y0 + t * y1;
}

int rt_iop_ccrr_pigment_eval(double lambda_nm, double chl_mg_m3, rt_iop_t* out)
{
    if (!out) return -1;
    out->a = out->b = out->bb = 0.0;
    const double c = chl_mg_m3;
    if (c <= 0.0) return 0;
    double Aphi = ccrr_interp_chl_absorption(lambda_nm, 0);
    double Ephi = ccrr_interp_chl_absorption(lambda_nm, 1);
    if (!isfinite(Aphi) || !isfinite(Ephi)) return -2;
    double a = Aphi * pow(c, Ephi);
    /* v1.052: scattering/backscatter switched to the smooth Morel & Maritorena
     * 2001 (MM01) Case-1 parameterization, replacing the previous
     * b = c_pig - a_pig construction (which imprinted the inverse of the
     * structured Bricaud absorption onto b -> spurious green-band peaks).
     * MM01 b_p and bb_p are spectrally smooth:
     *   b_p(550)     = 0.416 * Chl^0.766                  (MM01 eq.12; Loisel&Morel 1998)
     *   b_p(lambda)  = b_p(550) * (550/lambda)            (lambda^-1 dependence)
     *   bb_p(lambda) = [0.002 + 0.01*(0.5 - 0.25*log10 Chl)*(lambda/550)^v] * b_p(550)  (MM01 eq.13)
     *   v = 0.5*(log10 Chl - 0.3) for 0.02<Chl<2; v=0 for Chl>2; v=-1 for Chl<0.02 (MM01 eq.14)
     * Note: a_pig uses the configured CCRR Chl absorption table. bb/b is no longer
     * constant 0.006 but the
     * MM01 backscattering efficiency (wavelength- and Chl-dependent, ~0.3-1.1%). */
    double bp550 = 0.416 * pow(c, 0.766);
    double v_exp;
    if (c < 0.02)       v_exp = -1.0;
    else if (c <= 2.0)  v_exp = 0.5 * (log10(c) - 0.3);
    else                v_exp = 0.0;
    double b = bp550 * (550.0 / lambda_nm);
    double bb = (0.002 + 0.01 * (0.5 - 0.25 * log10(c)) * pow(lambda_nm / 550.0, v_exp)) * bp550;
    if (bb < 0.0) bb = 0.0;
    out->a = a;
    out->b = b;
    out->bb = bb;
    return 0;
}

/* CCRR/IOCCG-21 mineral (TSM) constituent IOPs  [all m^-1].
 *
 * Parameterization (CCRR reference package; constants below are the CCRR
 * values, not re-derived):
 *
 *   a_min(lambda)  = a*(443) * TSM * exp( -S_min * (lambda - 443) )
 *       a*(443) = 0.041 m^2 g^-1,  S_min = 0.0123 nm^-1
 *   c_min(555)     = a_min(555) + b*(555) * TSM,   b*(555) = 0.51 m^2 g^-1
 *   c_min(lambda)  = c_min(555) * (lambda/555)^(-gamma),  gamma = 0.3749
 *   b_min(lambda)  = max(0, c_min(lambda) - a_min(lambda))
 *   bb_min         = 0.026 * b_min        (fixed backscatter ratio 2.6%)
 *
 * Symbols / units:
 *   TSM = min_g_m3, total suspended mineral concentration [g m^-3]
 *   a*, b* = mass-specific absorption/scattering [m^2 g^-1]
 * Structure: absorption is an exponential in lambda; ATTENUATION c is the
 * power law; scattering is their difference (hence the max(0,.) guard for
 * far-red where the exponential tail could exceed the power law).
 * NOTE: a "+48% high-omega 412 nm delta-M defect" once recorded against
 * the consuming water-RT path was RETRACTED (no measurement record;
 * S-004, 2026-07-05).  These IOP formulas were never implicated either
 * way.  Validation matrix item #5 (mineral 4-species sweep) remains the
 * gate for any native-path truncation changes. */
int rt_iop_ccrr_mineral_eval(double lambda_nm, double min_g_m3, rt_iop_t* out)
{
    if (!out) return -1;
    out->a = out->b = out->bb = 0.0;
    const double m = min_g_m3;
    if (m <= 0.0) return 0;
    const double amin443_per_min = 0.041;
    const double mineral_abs_slope = 0.0123;
    const double bmin555_per_min = 0.51;
    const double cmin_power = 0.3749;
    double a443 = amin443_per_min * m;
    double a = a443 * exp(-mineral_abs_slope * (lambda_nm - 443.0));
    double a555 = a443 * exp(-mineral_abs_slope * (555.0 - 443.0));
    double c555 = a555 + bmin555_per_min * m;
    double c_raw = c555 * pow(lambda_nm / 555.0, -cmin_power);
    double b_raw = c_raw - a;
    double b = b_raw > 0.0 ? b_raw : 0.0;
    out->a = a;
    out->b = b;
    out->bb = 0.026 * b;
    return 0;
}

/* Backscatter fraction of the Henyey-Greenstein phase function  [0..0.5].
 *
 *   P_HG(theta; g) = (1 - g^2) / (1 + g^2 - 2 g cos theta)^(3/2)
 *   (normalized (1/2) int_{-1}^{1} P_HG dcos = 1)
 *
 * Closed-form backscatter fraction (integral of P_HG/2 over 90..180 deg,
 * i.e. cos theta in [-1, 0]; standard antiderivative of (A - 2g u)^(-3/2)):
 *
 *   bb/b (g) = (1 - g^2) / (2 g) * [ (1 + g^2)^(-1/2) - (1 + g)^(-1) ]
 *
 * Limits: g -> 0  => 1/2 (isotropic);  g -> 1 => 0 (fully forward).
 * The g -> 0 case is returned explicitly to avoid the 0/0 form. */
static double hg_backscatter_fraction(double g)
{
    if (fabs(g) < 1.0e-12) return 0.5;
    double a = 1.0 + g*g;
    return 0.5 * (1.0 - g*g) / g * (1.0 / sqrt(a) - 1.0 / (1.0 + g));
}

/* Inverse map: HG asymmetry g such that bb/b(g) = target.
 *
 * bb/b(g) is strictly decreasing on g in (0,1) from 0.5 to 0, so the
 * target must lie in (0, 0.5); out-of-range returns g = 0 (isotropic).
 * Solved by 100 bisection steps (interval width 1e-30 — far below double
 * precision noise; effectively exact).  Consumers: rt_water_rt.c uses this
 * to synthesize an HG proxy phase with a prescribed backscatter ratio
 * (delta2bb / hg-lut modes). */
double rt_iop_hg_g_for_backscatter_fraction(double target_bb_fraction)
{
    double t = target_bb_fraction;
    if (!(t > 0.0 && t < 0.5)) return 0.0;
    double lo = 0.0, hi = 0.999999;
    for (int it = 0; it < 100; ++it) {
        double mid = 0.5 * (lo + hi);
        if (hg_backscatter_fraction(mid) > t) lo = mid;
        else hi = mid;
    }
    return 0.5 * (lo + hi);
}


/* CCRR/IOCCG21 Chapter-4 scalar particle phase coefficients.
 *
 * These arrays are generated from ccrr_ocean_ref_package_v0.5.0:
 *   ccrr_ocean_ref.phase.particle_moment_array(..., lmax=24)
 * followed by Betal[l] = (2l+1) * chi_l, matching rt_kernel_phase_fourier().
 *
 * The port is intentionally pre-solver only: it injects scalar particle P11
 * coefficients into OCRT's aerosol/scalar slot while pure-water Rayleigh
 * polarization remains in the Rayleigh slot. Full particle Mueller matrices are
 * deferred until vector CCRR particle optics are available. */

/* Optional CCRR/IOCCG21 phase-moment CSV loader.
 * CSV schema expected from CCRR regenerated references:
 *   phase_name,l,chi_l[,bb_over_b_target,source_method,...]
 * Optional aliases accepted:
 *   phase/kind/name, ell/order, chi/moment, betal_l/betal/beta_l.
 *
 * Pure-water phase rows are ignored by design: OCRT CCRR-mode keeps the
 * OCRT vector Rayleigh-like seawater phase to preserve polarization support.
 */
#define CCRR_PHASE_LOADED_LMAX 24
static int    g_ccrr_phase_loaded = 0;
static char   g_ccrr_phase_loaded_path[1024] = {0};
static int    g_ccrr_phase_have[2] = {0, 0};
static int    g_ccrr_phase_lmax[2] = {-1, -1};
static double g_ccrr_phase_betal[2][CCRR_PHASE_LOADED_LMAX + 1];

static void ccrr_trim(char *s)
{
    if (!s) return;
    char *p = s;
    while (*p && isspace((unsigned char)*p)) ++p;
    if (p != s) memmove(s, p, strlen(p) + 1);
    size_t n = strlen(s);
    while (n > 0 && isspace((unsigned char)s[n-1])) s[--n] = '\0';
    if (n >= 2 && ((s[0] == '"' && s[n-1] == '"') || (s[0] == '\'' && s[n-1] == '\''))) {
        memmove(s, s + 1, n - 2);
        s[n - 2] = '\0';
    }
}

static void ccrr_lower(char *s)
{
    if (!s) return;
    for (; *s; ++s) *s = (char)tolower((unsigned char)*s);
}

static int ccrr_split_csv_line(char *line, char **fields, int max_fields)
{
    int n = 0;
    char *p = line;
    while (p && *p && n < max_fields) {
        fields[n++] = p;
        char *comma = strchr(p, ',');
        if (!comma) break;
        *comma = '\0';
        p = comma + 1;
    }
    for (int i = 0; i < n; ++i) ccrr_trim(fields[i]);
    return n;
}

static int ccrr_header_index(char **fields, int n, const char *name)
{
    for (int i = 0; i < n; ++i) {
        char tmp[128];
        snprintf(tmp, sizeof tmp, "%s", fields[i] ? fields[i] : "");
        ccrr_trim(tmp); ccrr_lower(tmp);
        if (!strcmp(tmp, name)) return i;
    }
    return -1;
}

static int ccrr_header_index_any(char **fields, int n, const char **names)
{
    for (int k = 0; names[k]; ++k) {
        int idx = ccrr_header_index(fields, n, names[k]);
        if (idx >= 0) return idx;
    }
    return -1;
}

static int ccrr_phase_kind_from_name(const char *phase_name)
{
    char tmp[256];
    snprintf(tmp, sizeof tmp, "%s", phase_name ? phase_name : "");
    ccrr_lower(tmp);
    if (strstr(tmp, "water") || strstr(tmp, "rayleigh")) return -3; /* explicitly ignored */
    if (strstr(tmp, "min") || strstr(tmp, "petzold") || strstr(tmp, "pet")) return 1;
    if (strstr(tmp, "pig") || strstr(tmp, "phyto") || strstr(tmp, "ff")) return 0;
    return -1;
}

int rt_iop_ccrr_phase_moments_load(const char* csv_path)
{
    if (!csv_path || !csv_path[0]) return -1;
    if (g_ccrr_phase_loaded && !strcmp(g_ccrr_phase_loaded_path, csv_path)) return 0;

    FILE *f = fopen(csv_path, "r");
    if (!f) {
        fprintf(stderr, "rt_iop_ccrr_phase_moments_load: cannot open '%s'\n", csv_path);
        return -2;
    }

    for (int k = 0; k < 2; ++k) {
        g_ccrr_phase_have[k] = 0;
        g_ccrr_phase_lmax[k] = -1;
        for (int l = 0; l <= CCRR_PHASE_LOADED_LMAX; ++l) g_ccrr_phase_betal[k][l] = 0.0;
    }

    char line[4096];
    if (!fgets(line, sizeof line, f)) { fclose(f); return -3; }
    char *fields[128];
    int nf = ccrr_split_csv_line(line, fields, 128);
    const char *phase_names[] = {"phase_name", "phase", "kind", "name", "component", NULL};
    const char *l_names[]     = {"l", "ell", "order", "legendre_order", NULL};
    const char *chi_names[]   = {"chi_l", "chi", "moment", "moment_l", "value", NULL};
    const char *betal_names[] = {"betal_l", "betal", "beta_l", "bet_l", NULL};
    int iphase = ccrr_header_index_any(fields, nf, phase_names);
    int il     = ccrr_header_index_any(fields, nf, l_names);
    int ichi   = ccrr_header_index_any(fields, nf, chi_names);
    int ibetal = ccrr_header_index_any(fields, nf, betal_names);
    if (iphase < 0 || il < 0 || (ichi < 0 && ibetal < 0)) {
        fclose(f);
        fprintf(stderr,
            "rt_iop_ccrr_phase_moments_load: '%s' missing required columns; "
            "need phase_name,l,chi_l or betal_l\n", csv_path);
        return -4;
    }

    int n_rows = 0;
    while (fgets(line, sizeof line, f)) {
        char *q = line;
        while (*q && isspace((unsigned char)*q)) ++q;
        if (!*q || *q == '#' || *q == '\n' || *q == '\r') continue;
        char *row[128];
        int nr = ccrr_split_csv_line(q, row, 128);
        if (nr <= iphase || nr <= il) continue;
        int kind = ccrr_phase_kind_from_name(row[iphase]);
        if (kind == -3) continue; /* pure water intentionally ignored */
        if (kind < 0 || kind > 1) continue;
        int l = atoi(row[il]);
        if (l < 0 || l > CCRR_PHASE_LOADED_LMAX) continue;
        double betal = 0.0;
        if (ibetal >= 0 && ibetal < nr && row[ibetal] && row[ibetal][0]) {
            betal = strtod(row[ibetal], NULL);
        } else if (ichi >= 0 && ichi < nr && row[ichi] && row[ichi][0]) {
            double chi = strtod(row[ichi], NULL);
            betal = (2.0 * (double)l + 1.0) * chi;
        } else {
            continue;
        }
        g_ccrr_phase_betal[kind][l] = betal;
        if (l > g_ccrr_phase_lmax[kind]) g_ccrr_phase_lmax[kind] = l;
        g_ccrr_phase_have[kind] = 1;
        n_rows++;
    }
    fclose(f);

    for (int kind = 0; kind < 2; ++kind) {
        if (g_ccrr_phase_have[kind]) {
            double b0 = g_ccrr_phase_betal[kind][0];
            if (isfinite(b0) && fabs(b0) > 0.0 && fabs(b0 - 1.0) > 1.0e-10) {
                for (int l = 0; l <= g_ccrr_phase_lmax[kind]; ++l)
                    g_ccrr_phase_betal[kind][l] /= b0;
                fprintf(stderr,
                    "warning: normalized CCRR phase kind=%d by Betal[0]=%.17g\n", kind, b0);
            }
        }
    }

    snprintf(g_ccrr_phase_loaded_path, sizeof g_ccrr_phase_loaded_path, "%s", csv_path);
    g_ccrr_phase_loaded = 1;
    if (n_rows <= 0 || (!g_ccrr_phase_have[0] && !g_ccrr_phase_have[1])) {
        fprintf(stderr,
            "rt_iop_ccrr_phase_moments_load: '%s' contained no pigment/mineral moments\n", csv_path);
        return -5;
    }
    return 0;
}

int rt_iop_ccrr_particle_betal(int kind, int n_out, double* betal_out)
{
    static const double pig_ff006_betal[] = {
        1.00000000000000000e+00,
        2.90750012926079471e+00,
        4.70793278658028491e+00,
        6.38661562443336894e+00,
        7.90585306645675345e+00,
        9.28250788200892529e+00,
        1.04976506848516369e+01,
        1.15163354688417492e+01,
        1.23361370937665793e+01,
        1.29713541408653725e+01,
        1.34201938303219670e+01,
        1.36621414097765381e+01,
        1.37182693490500824e+01,
        1.36053179123755719e+01,
        1.33344280737040375e+01,
        1.29055554231316112e+01,
        1.23457325580187671e+01,
        1.16719431567840815e+01,
        1.09030028348568955e+01,
        1.00532035146098426e+01,
        9.14440022379050887e+00,
        8.19172367638061338e+00,
        7.21530608731857193e+00,
        6.23454866949466879e+00,
        5.26207714102114288e+00
    };
    static const double min_petzold026_betal[] = {
        1.00000000000000000e+00,
        2.72878292758440821e+00,
        4.20985674554430922e+00,
        5.51211013032962605e+00,
        6.63506938685149183e+00,
        7.66777415933062034e+00,
        8.55358486974643384e+00,
        9.23020281473353243e+00,
        9.73762762328856368e+00,
        1.01222630384345127e+01,
        1.03662098779425058e+01,
        1.04063387781741792e+01,
        1.03076542197237853e+01,
        1.00965937357994626e+01,
        9.77653473445472443e+00,
        9.31767800966559889e+00,
        8.77310508412777779e+00,
        8.15433579489090654e+00,
        7.47932288395693057e+00,
        6.75184847407614086e+00,
        5.99773998498848027e+00,
        5.21948197553579973e+00,
        4.43983463120529898e+00,
        3.68256389550886531e+00,
        2.94895009618249748e+00
    };
    if (!betal_out || n_out <= 0) return -1;
    const double *src = NULL;
    int n_src = 0;
    if (kind == 0) {
        if (g_ccrr_phase_have[0]) {
            src = g_ccrr_phase_betal[0];
            n_src = g_ccrr_phase_lmax[0] + 1;
        } else {
            src = pig_ff006_betal;
            n_src = (int)(sizeof(pig_ff006_betal) / sizeof(pig_ff006_betal[0]));
        }
    } else if (kind == 1) {
        if (g_ccrr_phase_have[1]) {
            src = g_ccrr_phase_betal[1];
            n_src = g_ccrr_phase_lmax[1] + 1;
        } else {
            src = min_petzold026_betal;
            n_src = (int)(sizeof(min_petzold026_betal) / sizeof(min_petzold026_betal[0]));
        }
    } else {
        return -2;
    }
    for (int i = 0; i < n_out; ++i) {
        betal_out[i] = (i < n_src) ? src[i] : 0.0;
    }
    return n_src - 1;
}

/* Total bulk IOPs: constituents are INDEPENDENTLY ADDITIVE in a, b, bb
 * (standard IOP linearity; valid because absorption and scattering
 * coefficients are per-volume quantities and multiple-constituent
 * interference is neglected by construction of the constituent models).
 * The solver consumes only this total (plus a total phase function built
 * elsewhere as the b-weighted blend of constituent phases). */
rt_iop_t rt_iop_total(const rt_iop_components_t* comp)
{
    rt_iop_t t = { 0.0, 0.0, 0.0 };
    if (!comp) return t;
    t.a  = comp->pure_water.a  + comp->cdom.a  + comp->pigment.a +
           comp->eap_phyto.a   + comp->detritus.a + comp->mineral.a;
    t.b  = comp->pure_water.b  + comp->cdom.b  + comp->pigment.b +
           comp->eap_phyto.b   + comp->detritus.b + comp->mineral.b;
    t.bb = comp->pure_water.bb + comp->cdom.bb + comp->pigment.bb +
           comp->eap_phyto.bb  + comp->detritus.bb + comp->mineral.bb;
    return t;
}
