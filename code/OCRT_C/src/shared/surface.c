/* ============================================================================
 * shared/surface.c — Polarized Fresnel rough-ocean surface operators
 *
 * The implementation is organized from the physical definitions: Snell and
 * Fresnel amplitudes at a planar dielectric interface, an isotropic Gaussian
 * facet-slope distribution, microfacet reflection/transmission geometry, and
 * azimuthal Fourier projection of the resulting Mueller operators.  Scalar
 * and vector transmission paths share one geometry evaluator so their intensity
 * kernels cannot drift through duplicated formulas.
 * ============================================================================ */
#include "surface.h"
#include "surface_persistent_cache.h"
#include <stdlib.h>
#include "mat3.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* 자오면 회전각의 극 한계를 발동시키는 사인 문턱 (2026-07-24, 정합테스트 v4).
 *
 * Mishchenko 회전각 cos σ = (μ_상대 - μ_자기·cosΘ)/(sinΘ·sin θ_자기) 는
 * sin θ_자기 -> 0 에서 분자도 같은 차수로 0 이 된다.  0/0 이다.
 * 분모만 막으면 분자가 정확히 0 이 되는 극에서 회전이 방위각과 무관해지고,
 * 그 결과 천저에서 Q 에 있으면 안 되는 m=0 성분이 생기며 m=2 진폭이 약 2%
 * 낮아진다(2026-07-24 실측, 관측천정각 0 도 18조건).
 *
 * 전개하면 극한은 다음과 같다.  φ 는 두 방향의 방위각 차이다.
 *     cos σ -> -sign(μ_자기) · cos φ,      sin σ -> sin φ
 *
 * 문턱은 sqrt(DBL_EPSILON) ≈ 1.5e-8 근방으로 잡는다.  이보다 크면 원식의
 * 자릿수 상쇄 오차(≈ eps/sin θ)가 작고, 이보다 작으면 극한식의 절단 오차
 * (≈ sin θ)가 작다.  두 오차가 교차하는 지점이다.  실제 격자에 나타나는
 * sin θ 는 정확히 0 이거나 1e-2 이상이므로 문턱값에 민감하지 않다.
 *
 * [양쪽 방향이 동시에 극인 경우] 두 연직 방향 사이의 방위각 차는 정의되지
 * 않으므로 극한이 접근 경로에 의존한다.  위 식은 그 중 한 경로를 고른 것이다.
 * 이 모호성은 결과에 영향이 없음을 2026-07-24 에 실측으로 확인했다.
 *   - 세기 행·열은 이 자리에서도 영향을 받지 않는다.  II 가 자리 하나 바뀌지
 *     않고 IQ·IU·QI·UI 는 1e-13 이하다.  Q·U 블록의 기준면 회전만 달라진다.
 *   - 이 배치는 μ=1 이 입사와 출사 양쪽 지표로 동시에 나타날 때만 생기는데,
 *     μ=1 은 삽입 절점(관측 방향·특수 슬롯)으로만 링에 들어오고 삽입 절점의
 *     구적 가중치는 0 이다(rt_solver.c 의 atm->gb[+jt] = 0.0, 수중은
 *     rt_water_rt.h 의 "0 for zero-weight specials").  경계 조건 합에서
 *     이 쌍은 가중치 0 과 곱해진다.
 *   - 검증: 이 자리의 회전을 항등으로 바꾼 판(QQ 의 부호가 뒤집히고 크기가
 *     두 배가 되며 QU·UQ 결합이 사라지는 큰 교란)과 정본이 224 조건에서
 *     전부 비트 동일했다.  태양천정각 0~60, 관측천정각 0~40, 상대방위각
 *     0~180, 풍속 0~10, 443·555·865 nm, black_fresnel_ocean·coxmunk·ocean,
 *     생산급 수중 설정을 덮는다.  천저 수출광 Rrs0plus_Q 가 8.70e-04 로
 *     0 이 아닌 상태에서 확인했으므로 시험이 헛돈 것이 아니다. */
#define OCRT_SURF_POLE_EPS 1.0e-8

double surface_slope_variance(double wind_speed_ms, int sigma_type) {
    if (sigma_type == 1) {
        return 0.003 + 0.00512 * fmax(0.01, wind_speed_ms);
    }
    double sigma = 0.0731 * sqrt(fmax(0.0, wind_speed_ms));
    double sig2 = sigma * sigma;
    return fmax(sig2, 1e-10);
}

typedef struct {
    double sin_i;
    double sin_o;
    double incidence_cosine;
    double exit_dot_half;
    double facet_cosine;
    double slope_density;
    double half_denom_sq;
} surface_aw_microfacet_geometry_t;

/* Evaluate the air-to-water transmission half-vector and Beckmann/Cox-Munk
 * density.  Directions point away from the interface: the incident air ray is
 * upward and the transmitted water ray is downward.  Fresnel amplitudes and
 * Stokes rotations are deliberately left to the scalar/vector callers. */
static int surface_aw_microfacet_geometry_(double mu_i, double mu_o,
                                           double cphi, double sphi,
                                           double n_water,
                                           double slope_variance,
                                           surface_aw_microfacet_geometry_t *g)
{
    if (!g || slope_variance <= 1e-10 || mu_i <= 1e-9 || mu_o <= 1e-9)
        return -1;

    const double si = sqrt(fmax(0.0, 1.0 - mu_i * mu_i));
    const double so = sqrt(fmax(0.0, 1.0 - mu_o * mu_o));
    const double ix = si,        iy = 0.0,       iz =  mu_i;
    const double ox = so*cphi,   oy = so*sphi,   oz = -mu_o;

    double hx = -(ix + n_water * ox);
    double hy = -(iy + n_water * oy);
    double hz = -(iz + n_water * oz);
    const double hnorm = sqrt(hx*hx + hy*hy + hz*hz);
    if (hnorm < 1e-12) return -1;
    hx /= hnorm; hy /= hnorm; hz /= hnorm;
    if (hz < 0.0) { hx = -hx; hy = -hy; hz = -hz; }

    const double ih = ix*hx + iy*hy + iz*hz;
    const double oh = ox*hx + oy*hy + oz*hz;
    if (ih <= 0.0 || hz <= 1e-9) return -1;

    const double tan2 = (1.0 - hz*hz) / (hz*hz);
    const double D = exp(-tan2 / slope_variance) /
                     (M_PI * slope_variance * hz*hz*hz*hz);
    double denom = ih + n_water * oh;
    denom *= denom;

    g->sin_i = si;
    g->sin_o = so;
    g->incidence_cosine = ih;
    g->exit_dot_half = oh;
    g->facet_cosine = hz;
    g->slope_density = D;
    g->half_denom_sq = denom;
    return 0;
}

/* ----------------------------------------------------------------------------
 * surface_T_aw_coxmunk_direct — direct-beam air->water IRRADIANCE transmittance
 *   for a Cox-Munk rough surface (slope-integrated facet Fresnel transmission).
 *
 *   The slope integral accounts for facets whose local incidence differs from
 *   the global solar zenith, which is essential at grazing illumination.
 *
 *   T_aw(theta_s) = (1/mu_s) * INT p(zx,zy) * (mu_s + zx*sin_s) * T_F(theta_loc)
 *                              * H(mu_s + zx*sin_s > 0)  dzx dzy
 *     p(zx,zy)        = (1/(pi*sigma2)) exp(-(zx^2+zy^2)/sigma2)  Cox-Munk NDF
 *     cos(theta_loc)  = (mu_s + zx*sin_s)/sqrt(1+zx^2+zy^2)       local incidence
 *     (mu_s+zx*sin_s) = cos(theta_loc)*sqrt(1+..)  = facet projection / cos(beta)
 *     T_F             = 1 - R_unpol(theta_loc)  air->water Fresnel irrad. trans.
 *   Shadowing is omitted in this direct transmission integral.  As sigma2->0
 *   the expression tends to the planar value 1-R(theta_s).
 * -------------------------------------------------------------------------- */
double surface_T_aw_coxmunk_direct(double mu_sun_air, double n_water,
                                   double wind_speed, int sigma_type) {
    if (mu_sun_air <= 1e-9) return 0.0;
    double sig2 = surface_slope_variance(wind_speed, sigma_type);
    double sin_s = sqrt(fmax(0.0, 1.0 - mu_sun_air * mu_sun_air));
    if (sig2 <= 1e-10) {
        /* flat limit: 1 - R_unpol(theta_sun) */
        double cti0 = mu_sun_air;
        double stt0 = sin_s / n_water;
        double ctt0 = sqrt(fmax(0.0, 1.0 - stt0 * stt0));
        double rs0  = (cti0 - n_water * ctt0) / (cti0 + n_water * ctt0);
        double rp0  = (n_water * cti0 - ctt0) / (n_water * cti0 + ctt0);
        return 1.0 - 0.5 * (rs0 * rs0 + rp0 * rp0);
    }
    double sig = sqrt(sig2);
    double L   = 6.0 * sig;
    const int N = 101;
    double dz  = 2.0 * L / (double)(N - 1);
    double acc = 0.0;
    for (int i = 0; i < N; ++i) {
        double zx  = -L + i * dz;
        double wxi = (i == 0 || i == N - 1) ? 0.5 : 1.0;   /* trapezoid weight */
        double proj = mu_sun_air + zx * sin_s;             /* facet projection toward sun */
        if (proj <= 0.0) continue;                         /* facet not illuminated */
        for (int j = 0; j < N; ++j) {
            double zy   = -L + j * dz;
            double wyj  = (j == 0 || j == N - 1) ? 0.5 : 1.0;
            double r2   = zx * zx + zy * zy;
            double cti  = proj / sqrt(1.0 + r2);           /* cos(theta_local) */
            double sti2 = fmax(0.0, 1.0 - cti * cti);
            double stt  = sqrt(sti2) / n_water;
            double ctt  = sqrt(fmax(0.0, 1.0 - stt * stt));
            double rs   = (cti - n_water * ctt) / (cti + n_water * ctt);
            double rp   = (n_water * cti - ctt) / (n_water * cti + ctt);
            double TF   = 1.0 - 0.5 * (rs * rs + rp * rp);
            double p    = (1.0 / (M_PI * sig2)) * exp(-r2 / sig2);
            acc += wxi * wyj * p * proj * TF;
        }
    }
    acc *= dz * dz;
    double T_aw = acc / mu_sun_air;
    if (T_aw < 0.0) T_aw = 0.0;
    if (T_aw > 1.0) T_aw = 1.0;
    return T_aw;
}

/* ----------------------------------------------------------------------------
 * surface_T_aw_coxmunk_btdf_scalar — air->water microfacet transmission BTDF
 *   (Walter et al. 2007, "Microfacet Models for Refraction", eq. 21), scalar
 *   intensity, Cox-Munk (Beckmann) slope distribution, no shadowing (G=1).
 *
 *   Geometry (directions point AWAY from the surface point):
 *     i = (sin θ_i, 0, μ_i)            incident air dir, toward sky (up), φ_i=0
 *     o = (sin θ_o cos Δφ, sin θ_o sin Δφ, -μ_o)   transmitted water dir (down)
 *     h = -(η_i i + η_o o)/|·|   transmission half-vector (η_i=1, η_o=n), up
 *   f_t = (|i·h||o·h|)/(μ_i μ_o) · η_o²(1-F(i·h)) D(h) / (η_i(i·h)+η_o(o·h))²
 *     D(h) = exp(-tan²θ_h/σ²)/(π σ² cos⁴θ_h)   Beckmann (α²=σ²=Cox-Munk MSS)
 *   The η_o²=n² factor is the in-water radiance enhancement.  Validated against
 *   surface_T_aw_coxmunk_direct: ∫ f_t cos(θ_o) dΩ_o == T_aw_direct(μ_i)
 *   (ratio 1.00 for θ_i>=30°, agreement limited only by output-grid resolution
 *   near θ_i->0 where the transmission lobe collapses to the Snell direction).
 * -------------------------------------------------------------------------- */
double surface_T_aw_coxmunk_btdf_scalar(double mu_i, double mu_o, double dphi,
                                        double n_water, double wind_speed,
                                        int sigma_type) {
    double s2 = surface_slope_variance(wind_speed, sigma_type);
    surface_aw_microfacet_geometry_t g;
    if (surface_aw_microfacet_geometry_(mu_i, mu_o, cos(dphi), sin(dphi),
                                        n_water, s2, &g) != 0) return 0.0;
    double n = n_water;
    /* Fresnel reflectance at local incidence (air->water). */
    double cti = g.incidence_cosine;
    double sti = sqrt(fmax(0.0, 1.0 - cti * cti));
    double stt = sti / n, F;
    if (stt >= 1.0) {
        F = 1.0;
    } else {
        double ctt = sqrt(fmax(0.0, 1.0 - stt * stt));
        double rs = (cti - n * ctt) / (cti + n * ctt);
        double rp = (n * cti - ctt) / (n * cti + ctt);
        F = 0.5 * (rs * rs + rp * rp);
    }
    double denom = g.half_denom_sq;
    if (denom < 1e-12) return 0.0;
    double ft = (fabs(g.incidence_cosine) * fabs(g.exit_dot_half)) /
                (mu_i * mu_o)
              * (n * n * (1.0 - F) * g.slope_density) / denom;
    if (ft < 0.0) ft = 0.0;
    return ft;
}

static inline double surface_water_to_air_unpolarized_reflectance_(double cti,
                                                                      double n_water) {
    double sti2 = fmax(0.0, 1.0 - cti * cti);
    double sair = n_water * sqrt(sti2);
    if (sair >= 1.0) return 1.0;
    double cair = sqrt(fmax(0.0, 1.0 - sair * sair));
    double rs = (n_water * cti - cair) / (n_water * cti + cair);
    double rp = (cti - n_water * cair) / (cti + n_water * cair);
    return 0.5 * (rs * rs + rp * rp);
}

/* ----------------------------------------------------------------------------
 * surface_R_ww_coxmunk_direct — water->water (internal) IRRADIANCE reflectance
 *   for an upwelling beam at μ_up_water on a Cox-Munk rough surface (slope-
 *   integrated facet water->air Fresnel reflectance, including total internal
 *   reflection beyond the critical angle).  No shadowing (matches the air->water
 *   direct transmittance convention).  σ²->0 gives the flat internal reflectance
 *   R_wa(μ_up) (=1 for μ_up < μ_crit = sqrt(1-1/n²) = TIR).  By energy
 *   conservation surface_R_ww_coxmunk_direct + (water->air transmittance) = 1.
 *
 *   R_ww(μ_up) = (1/μ_up) ∫ p(zx,zy) (μ_up + zx sin_up) R_wa(θ_loc)
 *                          H(μ_up + zx sin_up > 0) dzx dzy
 *   cos θ_loc = (μ_up + zx sin_up)/sqrt(1+zx²+zy²)
 *   R_wa: incidence from water (n) to air (1); sin θ_air = n sin θ_loc;
 *         sin θ_air >= 1 -> R=1 (TIR). -------------------------------------- */
#ifdef OCRT_FAST_KERNELS
static double surface_R_ww_coxmunk_direct_eval_(double, double, double, int);
#endif
double surface_R_ww_coxmunk_direct(double mu_up_water, double n_water,
                                   double wind_speed, int sigma_type) {
#ifdef OCRT_FAST_KERNELS
    /* v1.09-opt S2 (bit-safe memo): within a solve this is called ~1e3x per
     * node with mu drawn from the fixed quadrature set and (n_water, wind,
     * sigma_type) constant.  Caching the previously computed double returns
     * the EXACT same bits the direct evaluation produced.  Thread-local so
     * OMP batch threads never share entries. */
    /* open-addressing, linear probe, NEVER evict: the key set per process is
     * small (node mu x wind classes); a direct-mapped slot ping-pong measured
     * a 32% miss rate, so colliding keys must coexist. Probe cap then fall
     * through to a plain eval (correct, uncached). */
    enum { OCRT_RWWD_N = 512, OCRT_RWWD_PROBE = 16 };
    static _Thread_local struct {
        double mu, nw, ws, val; int sig, used;
    } ocrt_rwwd_cache[OCRT_RWWD_N];
    {
        unsigned long long h;
        memcpy(&h, &mu_up_water, sizeof(h));
        h ^= h >> 33; h *= 0xff51afd7ed558ccdULL; h ^= h >> 33;
        unsigned idx = (unsigned)(h % OCRT_RWWD_N);
        int free_slot = -1;
        for (int pr = 0; pr < OCRT_RWWD_PROBE; pr++) {
            unsigned i2 = (idx + (unsigned)pr) % OCRT_RWWD_N;
            if (!ocrt_rwwd_cache[i2].used) { if (free_slot < 0) free_slot = (int)i2; break; }
            if (ocrt_rwwd_cache[i2].mu == mu_up_water &&
                ocrt_rwwd_cache[i2].nw == n_water &&
                ocrt_rwwd_cache[i2].ws == wind_speed &&
                ocrt_rwwd_cache[i2].sig == sigma_type) {
                return ocrt_rwwd_cache[i2].val;
            }
        }
        double v = surface_R_ww_coxmunk_direct_eval_(mu_up_water, n_water,
                                                     wind_speed, sigma_type);
        if (free_slot >= 0) {
            ocrt_rwwd_cache[free_slot].mu = mu_up_water;
            ocrt_rwwd_cache[free_slot].nw = n_water;
            ocrt_rwwd_cache[free_slot].ws = wind_speed;
            ocrt_rwwd_cache[free_slot].sig = sigma_type;
            ocrt_rwwd_cache[free_slot].val = v;
            ocrt_rwwd_cache[free_slot].used = 1;
        }
        return v;
    }
}
static double surface_R_ww_coxmunk_direct_eval_(double mu_up_water, double n_water,
                                                double wind_speed, int sigma_type) {
#endif
    if (mu_up_water <= 1e-9) return 1.0;          /* grazing up-beam: fully reflected */
    double s2  = surface_slope_variance(wind_speed, sigma_type);
    double sin_u = sqrt(fmax(0.0, 1.0 - mu_up_water * mu_up_water));
    if (s2 <= 1e-10) {
        return surface_water_to_air_unpolarized_reflectance_(mu_up_water, n_water);
    }
    /* v1.10 OPTION-2 CLOSURE (2026-07-10): this function's ledger is now the
     * PHYSICAL water-side reflectance R_eff = 1 - T_up, where T_up is the
     * multi-encounter upward escape from the cascade engine (re-entrant
     * transmitted light that refracts back INTO the water is booked as
     * reflection, exactly matching its role in the Ed(0-) bookkeeping at the
     * single call site).  Uniform at all angles; below the critical cone the
     * correction vanishes smoothly (cascade returns T_up = T_h there). */
    /* SPRINT GUARD (2026-07-10): the closed-ledger cascade is OPT-IN ONLY
     * during the SPEED sprint (OCRT_MB_CLOSURE=1).  Default off restores the
     * pre-surgery first-encounter integral bit-exactly AND avoids a ~400k-ray
     * cascade per call on the STRICT build (which has no S2 memo).  getenv
     * snapshotted once per thread (commit #21: repeated getenv = lock storm).
     * Accuracy work parked - see BUILD_SPEEDUP carry-over list. */
    {
        static _Thread_local int mb_on = -1;
        if (mb_on < 0) mb_on = (getenv("OCRT_MB_CLOSURE") != NULL);
        if (mb_on) {
            double Th_c = 0.0, Tup_c = 0.0, Wret_c = 0.0;
            if (surface_wa_multibounce_cascade(n_water, wind_speed, mu_up_water,
                                               400000L, &Th_c, &Tup_c, &Wret_c,
                                               NULL, NULL, NULL, 0) == 0)
                return 1.0 - Tup_c;
            /* engine failure — fall through to first-encounter form */
        }
    }
    double sig = sqrt(s2);
    double L   = 6.0 * sig;
    const int N = 101;
    double dz  = 2.0 * L / (double)(N - 1);
    double acc = 0.0;
    for (int i = 0; i < N; ++i) {
        double zx   = -L + i * dz;
        double wxi  = (i == 0 || i == N - 1) ? 0.5 : 1.0;
        double proj = mu_up_water + zx * sin_u;
        if (proj <= 0.0) continue;
        for (int j = 0; j < N; ++j) {
            double zy  = -L + j * dz;
            double wyj = (j == 0 || j == N - 1) ? 0.5 : 1.0;
            double r2  = zx * zx + zy * zy;
            double cti = proj / sqrt(1.0 + r2);
            double Rwa = surface_water_to_air_unpolarized_reflectance_(cti, n_water);
            double p   = (1.0 / (M_PI * s2)) * exp(-r2 / s2);
            acc += wxi * wyj * p * proj * Rwa;
        }
    }
    acc *= dz * dz;
    double Rww = acc / mu_up_water;
    if (Rww < 0.0) Rww = 0.0;
    if (Rww > 1.0) Rww = 1.0;
    return Rww;
}

/* ----------------------------------------------------------------------------
 * Fresnel core helper (Phase B.2)
 *   Compute Snell-refracted angle and Fresnel s/p reflection amplitudes for a
 *   flat interface between media of indices n1 and n2. Used by both the
 *   reflection and transmission Mueller routines and (after refactor) by the
 *   Cox-Munk microfacet Fresnel call.
 *
 *   Convention: amplitudes follow the standard form
 *     rs = (n1*mu_i - n2*mu_t) / (n1*mu_i + n2*mu_t)
 *     rp = (n2*mu_i - n1*mu_t) / (n2*mu_i + n1*mu_t)
 *   With n1=1, n2=n_water this reduces exactly to the pre-Phase-B.2 inline
 *   formulas in surface_flat_fresnel_matrix_raw (line 48-49 of pre-refactor
 *   source) — refactor is therefore bit-exact for the air-side reflection.
 * ---------------------------------------------------------------------------- */
typedef struct {
    double mu_i;     /* incidence cosine in medium 1 (>= 0) */
    double mu_t;     /* transmission cosine in medium 2 (>= 0, 0 if TIR) */
    double rs_amp;   /* s-pol reflection amplitude */
    double rp_amp;   /* p-pol reflection amplitude */
    int    is_TIR;   /* 1 if total internal reflection (n1>n2 case only) */
} fresnel_core_t;

static void fresnel_compute_core(double mu_i, double n1, double n2,
                                  fresnel_core_t *fc) {
    fc->mu_i = mu_i;
    double sin_i_sq = fmax(0.0, 1.0 - mu_i * mu_i);
    /* Snell: n1 sin_i = n2 sin_t  =>  sin_t = (n1/n2) sin_i
     * sin_t_sq = (n1/n2)^2 * sin_i_sq */
    double sin_t_sq = (n1 / n2) * (n1 / n2) * sin_i_sq;
    if (sin_t_sq >= 1.0) {
        fc->is_TIR  = 1;
        fc->mu_t    = 0.0;
        fc->rs_amp  = 0.0;  /* |r| = 1 in TIR but amplitudes are complex; */
        fc->rp_amp  = 0.0;  /* caller (reflection branch) returns identity. */
        return;
    }
    fc->is_TIR = 0;
    fc->mu_t   = sqrt(1.0 - sin_t_sq);
    /* General Fresnel amplitudes */
    fc->rs_amp = (n1 * mu_i - n2 * fc->mu_t) / (n1 * mu_i + n2 * fc->mu_t);
    fc->rp_amp = (n2 * mu_i - n1 * fc->mu_t) / (n2 * mu_i + n1 * fc->mu_t);
}

/* Reflection Mueller matrix 3x3 in meridian basis, built from Fresnel core. */
static void fresnel_R_mueller_from_core(const fresnel_core_t *fc,
                                         int q_convention, double *MR) {
    mat3_zero(MR);
    if (fc->is_TIR) {
        MR[0] = MR[4] = MR[8] = 1.0;  /* identity (TIR) */
        return;
    }
    double Rs = fc->rs_amp * fc->rs_amp;
    double Rp = fc->rp_amp * fc->rp_amp;
    double Q_kernel = (q_convention == 1) ? 0.5 * (Rp - Rs) : 0.5 * (Rs - Rp);
    MR[0*3+0] = 0.5 * (Rs + Rp);
    MR[0*3+1] = Q_kernel;
    MR[1*3+0] = Q_kernel;
    MR[1*3+1] = 0.5 * (Rs + Rp);
    MR[2*3+2] = fc->rs_amp * fc->rp_amp;
}

/* Radiance transmission Mueller matrix 3x3, built from Fresnel core.
 *   Fresnel transmission amplitudes derived from reflection amplitudes via
 *   identities:  ts_amp = 1 + rs_amp,  tp_amp = (n1/n2)(1 + rp_amp).
 *   Flux transmittance: T_X_flux = (n2*mu_t)/(n1*mu_i) * t_X^2
 *   Radiance n^2 law:   T_X_rad  = (n2/n1)^2 * T_X_flux
 *   Energy conservation: R_X + T_X_flux = 1 per polarization (verified).
 */
static void fresnel_T_mueller_radiance_from_core(const fresnel_core_t *fc,
                                                  double n1, double n2,
                                                  int q_convention,
                                                  double *MT) {
    mat3_zero(MT);
    if (fc->is_TIR) {
        return;  /* zero matrix - no transmission */
    }
    double ts_amp = 1.0 + fc->rs_amp;
    double tp_amp = (n1 / n2) * (1.0 + fc->rp_amp);

    double flux_factor = (n2 * fc->mu_t) / (n1 * fc->mu_i);
    double n_ratio_sq  = (n2 * n2) / (n1 * n1);
    double radiance_factor = n_ratio_sq * flux_factor;

    double Ts = radiance_factor * ts_amp * ts_amp;
    double Tp = radiance_factor * tp_amp * tp_amp;

    double Q_kernel = (q_convention == 1) ? 0.5 * (Tp - Ts) : 0.5 * (Ts - Tp);
    MT[0*3+0] = 0.5 * (Ts + Tp);
    MT[0*3+1] = Q_kernel;
    MT[1*3+0] = Q_kernel;
    MT[1*3+1] = 0.5 * (Ts + Tp);
    MT[2*3+2] = radiance_factor * ts_amp * tp_amp;
}

/* ----------------------------------------------------------------------------
 * Flat-water Fresnel Mueller (3x3, meridian basis) — refactored
 *   Backward-compatible wrapper. Calls Fresnel core with (n1, n2) = (1, n_water)
 *   to reproduce the pre-Phase-B.2 air-side reflection Mueller bit-exactly.
 * ---------------------------------------------------------------------------- */
void surface_flat_fresnel_matrix_raw(double mu, double n_water,
                                      int q_convention, double *MF) {
    fresnel_core_t fc;
    fresnel_compute_core(mu, 1.0, n_water, &fc);
    fresnel_R_mueller_from_core(&fc, q_convention, MF);
}

/* ----------------------------------------------------------------------------
 * Generalized R / T Mueller matrices (Phase B.2 public API)
 * ---------------------------------------------------------------------------- */
void surface_flat_fresnel_R_matrix_general(double mu_i, double n1, double n2,
                                            int q_convention, double *MR) {
    fresnel_core_t fc;
    fresnel_compute_core(mu_i, n1, n2, &fc);
    fresnel_R_mueller_from_core(&fc, q_convention, MR);
}

void surface_flat_fresnel_T_matrix(double mu_i, double n1, double n2,
                                    int q_convention, double *MT) {
    fresnel_core_t fc;
    fresnel_compute_core(mu_i, n1, n2, &fc);
    fresnel_T_mueller_radiance_from_core(&fc, n1, n2, q_convention, MT);
}

double surface_flat_fresnel_check_energy(double mu_i, double n1, double n2) {
    fresnel_core_t fc;
    fresnel_compute_core(mu_i, n1, n2, &fc);
    if (fc.is_TIR) return 0.0;
    double Rs = fc.rs_amp * fc.rs_amp;
    double Rp = fc.rp_amp * fc.rp_amp;
    double ts_amp = 1.0 + fc.rs_amp;
    double tp_amp = (n1 / n2) * (1.0 + fc.rp_amp);
    double flux_factor = (n2 * fc.mu_t) / (n1 * fc.mu_i);
    double Ts_flux = flux_factor * ts_amp * ts_amp;
    double Tp_flux = flux_factor * tp_amp * tp_amp;
    double res_s = fabs(1.0 - (Rs + Ts_flux));
    double res_p = fabs(1.0 - (Rp + Tp_flux));
    return res_s + res_p;
}

/* ----------------------------------------------------------------------------
 * Cox-Munk + Fresnel BRDF kernel (3x3 Mueller R, full per-direction)
 *   V1 R_sfc_coxmunk_trig (vrt_solver.c L2983)
 *
 * Returns R such that I_up = R · I_dn at one (μ_o, μ_i, φ) direction triple.
 * Includes:
 *   - isotropic Gaussian wave-slope distribution (isotropic, selectable slope-variance law (sigma_type 0/1))
 *   - Sancer bistatic shadowing factor S_bi
 *   - Fresnel Mueller at microfacet incidence angle (omega = half-angle)
 *   - Hovenier-style basis rotation (incoming meridian → microfacet
 *     scatter plane → outgoing meridian) via L1, L2 matrices
 * ---------------------------------------------------------------------------- */
void surface_R_coxmunk_trig(double mu_out, double mu_in_dn,
                             double cos_phi, double sin_phi,
                             double ws, int sigma_type_int, double n_water,
                             int q_convention, double *R) {
    mat3_zero(R);

    double sigma_sq = surface_slope_variance(ws, sigma_type_int);

    double mu_in_actual = -mu_in_dn;
    double s1 = sqrt(fmax(0.0, 1.0 - mu_out * mu_out));
    double s2 = sqrt(fmax(0.0, 1.0 - mu_in_actual * mu_in_actual));

    double cT = mu_out * mu_in_actual + s1 * s2 * cos_phi;
    if      (cT > +1.0) cT = +1.0;
    else if (cT < -1.0) cT = -1.0;
    double cos_omega = sqrt(fmax(0.0, 0.5 * (1.0 - cT)));
    if (cos_omega <= 1e-6) return;
    double cos_beta = (mu_out + mu_in_dn) / (2.0 * cos_omega);
    if (cos_beta <= 1e-6) return;

    double sT = sqrt(fmax(0.0, 1.0 - cT * cT));
    double denom = fmax(1e-12, sT);

    /* Stokes rotation angles (Refactoring Tasks Task 2 fix, 2026-05-10):
     *   표준 (Mishchenko 2002 Appendix B):
     *     cos σ_1 = (μ_out - μ_in·cos Θ) / (sin Θ · sin θ_in)
     *     cos σ_2 = (μ_in - μ_out·cos Θ) / (sin Θ · sin θ_out)
     *   기존 V3 식 가 분모 에 sin θ 누락 — typical case s ≈ 1 에서 small,
     *   asymmetric geometry (sza ≠ vza) 에서 큰 deviation 발생 가능.
     *   Fix: ci1 분모 에 s2 (= sin θ_in), ci2 분모 에 s1 (= sin θ_out) 추가. */
    double s2_safe = fmax(1e-12, s2);
    double s1_safe = fmax(1e-12, s1);

    /* 극 한계 (2026-07-24, 정합테스트 v4).  OCRT_SURF_POLE_EPS 정의부의 주석을
     * 함께 볼 것.  극이 아닌 각도의 식은 이전 판과 글자 하나 다르지 않다. */
    double ci1, si1, ci2, si2;
    if (s2 <= OCRT_SURF_POLE_EPS) {
        ci1 = -copysign(1.0, mu_in_actual) * cos_phi;
        si1 = sin_phi;
    } else {
        ci1 = (mu_out - mu_in_actual * cT) / (denom * s2_safe);
        si1 = s1 * sin_phi / denom;
        double n1 = sqrt(ci1*ci1 + si1*si1);
        if (n1 < 1e-12) { ci1 = 1.0; si1 = 0.0; }
        else { ci1 /= n1; si1 /= n1; }
    }

    if (s1 <= OCRT_SURF_POLE_EPS) {
        ci2 = -copysign(1.0, mu_out) * cos_phi;
        si2 = sin_phi;
    } else {
        ci2 = (mu_in_actual - mu_out * cT) / (denom * s1_safe);
        si2 = s2 * sin_phi / denom;
        double n2 = sqrt(ci2*ci2 + si2*si2);
        if (n2 < 1e-12) { ci2 = 1.0; si2 = 0.0; }
        else { ci2 /= n2; si2 /= n2; }
    }

    /* Fresnel Mueller at microfacet incidence angle omega
     * (Phase B.2 refactor: use shared fresnel core helper. Bit-exact to the
     *  pre-refactor inline computation for n1=1, n2=n_water.) */
    double MF[9];
    {
        fresnel_core_t fc;
        fresnel_compute_core(cos_omega, 1.0, n_water, &fc);
        fresnel_R_mueller_from_core(&fc, q_convention, MF);
    }

    /* Cox-Munk isotropic Gaussian slope */
    double cos_beta_sq = cos_beta * cos_beta;
    double tan_beta_sq = (1.0 - cos_beta_sq) / fmax(cos_beta_sq, 1e-8);
    double P_slope = (1.0 / (M_PI * sigma_sq)) * exp(-tan_beta_sq / sigma_sq);

    /* Sancer (1969) bistatic shadowing for an isotropic Gaussian slope
     * distribution.  OCRT uses the additive two-direction combination
     * S = 1/(1 + lambda_i + lambda_v). */
    double sigma_len = sqrt(sigma_sq);
    double lam_i = 0.0, lam_v = 0.0;

    if (mu_in_dn < 1.0 - 1e-10) {
        double sin_i = sqrt(fmax(1e-20, 1.0 - mu_in_dn * mu_in_dn));
        double nu_i = (mu_in_dn / sin_i) / sigma_len;
        if (nu_i <= 20.0) {
            lam_i = 0.5 * (exp(-nu_i * nu_i) / (nu_i * sqrt(M_PI)) - erfc(nu_i));
        }
    }
    if (mu_out < 1.0 - 1e-10) {
        double sin_v = sqrt(fmax(1e-20, 1.0 - mu_out * mu_out));
        double nu_v = (mu_out / sin_v) / sigma_len;
        if (nu_v <= 20.0) {
            lam_v = 0.5 * (exp(-nu_v * nu_v) / (nu_v * sqrt(M_PI)) - erfc(nu_v));
        }
    }
    double S_bi = 1.0 / (1.0 + lam_i + lam_v);

    /* Cox-Munk/Wang-Gordon radiance BRDF. */
    double brdf = (P_slope * S_bi) / (4.0 * mu_in_dn * mu_out * cos_beta_sq * cos_beta_sq);

    /* Assemble: R = L2 · MF · L1 · brdf */
    double L1[9], L2[9], MFL1[9], L2MFL1[9];
    build_rotation_L(ci1, si1, L1);
    build_rotation_L(ci2, si2, L2);
    mat3_mul(MF, L1, MFL1);
    mat3_mul(L2, MFL1, L2MFL1);
    for (int i = 0; i < 9; ++i) R[i] = L2MFL1[i] * brdf;
}

/* ----------------------------------------------------------------------------
 * Cox-Munk + Fresnel BTDF kernel for water -> air TRANSMISSION (3x3 Mueller T)
 *
 * Transmission analog of surface_R_coxmunk_trig. Returns T such that
 *   I_air(mu_out_air) = T . I_water(mu_in_water)   at one (mu2, mu1, phi) triple.
 *
 * Geometry (refraction facet, not the reflection half-vector):
 *   p1 = in-water upwelling propagation dir (cosine mu1 > 0)
 *   p2 = air outgoing  propagation dir (cosine mu2 > 0)
 *   Snell tangential-momentum continuity => facet normal along
 *       N = n_w * p1 - n_a * p2      (n_a = 1)
 *   |N|^2 = n^2 + 1 - 2 n cosPsi,   cosPsi = p1.p2 = mu1 mu2 + s1 s2 cos_phi
 *   cos_beta = (n mu1 - mu2)/|N|              (facet tilt from vertical)
 *   cos_omega_w = (n - cosPsi)/|N|            (in-water incidence on facet)
 *   Snell at facet: sin_omega_a = n sin_omega_w ; TIR (>crit) => zero transmit.
 *
 * BTDF (Walter et al. 2007, radiance form), same NDF/shadowing framework as the
 * reflection kernel (D = P_slope / cos^4 beta, Sancer S_bi, Hovenier rotations):
 *   f_t,X = T_flux,X(omega_w) * (mu_wf mu_af)/(mu1 mu2)
 *               * eta_o^2 / (eta_i mu_wf + eta_o mu_af)^2 * D * G
 *         = t_X^2 * mu_af^2 / (n mu1 mu2 (n mu_wf + mu_af)^2) * (P/cos^4beta) * G
 *   with eta_i=n_w, eta_o=1, mu_wf=cos_omega_w, mu_af=cos_omega_a,
 *   t_s = 1+rs, t_p = n_w(1+rp)  (Fresnel transmission amplitudes from core).
 * Reuses fresnel_compute_core for the per-facet Fresnel; reuses
 * surface_slope_variance, the Sancer shadowing, and build_rotation_L exactly as
 * the reflection branch, so R and T share one consistent surface model.
 * ---------------------------------------------------------------------------- */
void surface_T_coxmunk_trig(double mu_out_air, double mu_in_water,
                            double cos_phi, double sin_phi,
                            double ws, int sigma_type_int, double n_water,
                            int q_convention, double *T) {
    mat3_zero(T);

    double mu1 = mu_in_water;   /* in-water upwelling cosine (>0) */
    double mu2 = mu_out_air;    /* air outgoing cosine        (>0) */
    if (mu1 <= 1e-9 || mu2 <= 1e-9) return;

    double sigma_sq = surface_slope_variance(ws, sigma_type_int);
    double n = n_water;

    double s1 = sqrt(fmax(0.0, 1.0 - mu1 * mu1));
    double s2 = sqrt(fmax(0.0, 1.0 - mu2 * mu2));

    double cosPsi = mu1 * mu2 + s1 * s2 * cos_phi;
    if      (cosPsi > +1.0) cosPsi = +1.0;
    else if (cosPsi < -1.0) cosPsi = -1.0;

    double Nsq = n * n + 1.0 - 2.0 * n * cosPsi;
    if (Nsq <= 1e-12) return;
    double Nmag = sqrt(Nsq);

    double cos_beta    = (n * mu1 - mu2) / Nmag;     /* facet tilt cosine */
    double cos_omega_w = (n - cosPsi)    / Nmag;     /* in-water facet incidence */
    if (cos_beta <= 1e-6 || cos_omega_w <= 1e-6) return;

    /* Per-facet Fresnel (reuse core): water(n) -> air(1) at incidence omega_w */
    fresnel_core_t fc;
    fresnel_compute_core(cos_omega_w, n, 1.0, &fc);
    if (fc.is_TIR) return;                 /* beyond critical angle: no transmit */
    double mu_af = fc.mu_t;                /* = cos_omega_a (Snell) */

    double ts = 1.0 + fc.rs_amp;           /* s transmission amplitude */
    double tp = n * (1.0 + fc.rp_amp);     /* p transmission amplitude (n1/n2=n) */

    /* Cox-Munk isotropic Gaussian slope NDF: D = P_slope / cos^4 beta */
    double cos_beta_sq  = cos_beta * cos_beta;
    double tan_beta_sq  = (1.0 - cos_beta_sq) / fmax(cos_beta_sq, 1e-8);
    double P_slope = (1.0 / (M_PI * sigma_sq)) * exp(-tan_beta_sq / sigma_sq);

    /* Sancer (1969) bistatic shadowing, identical to the reflection branch. */
    double sigma_len = sqrt(sigma_sq);
    double lam_i = 0.0, lam_v = 0.0;
    if (mu1 < 1.0 - 1e-10) {
        double sin_i = sqrt(fmax(1e-20, 1.0 - mu1 * mu1));
        double nu_i = (mu1 / sin_i) / sigma_len;
        if (nu_i <= 20.0)
            lam_i = 0.5 * (exp(-nu_i * nu_i) / (nu_i * sqrt(M_PI)) - erfc(nu_i));
    }
    if (mu2 < 1.0 - 1e-10) {
        double sin_v = sqrt(fmax(1e-20, 1.0 - mu2 * mu2));
        double nu_v = (mu2 / sin_v) / sigma_len;
        if (nu_v <= 20.0)
            lam_v = 0.5 * (exp(-nu_v * nu_v) / (nu_v * sqrt(M_PI)) - erfc(nu_v));
    }
    double S_bi = 1.0 / (1.0 + lam_i + lam_v);

    /* Walter radiance BTDF common factor (per-polarization scaled by t_X^2).
     * Denominator (eta_i(i.h)+eta_o(o.h))^2 reduces, with signed dot products,
     * to |N|^2 = Nsq exactly:  eta_i(i.h)+eta_o(o.h) = -(n^2+1-2n cosPsi)/|N| = -|N|.
     * (Using +n*mu_wf+mu_af here would be a sign error, ~30x too small.) */
    double common = (mu_af * mu_af) /
                    (n * mu1 * mu2 * Nsq) *
                    (P_slope / (cos_beta_sq * cos_beta_sq)) * S_bi;

    double ft_s  = common * ts * ts;
    double ft_p  = common * tp * tp;
    double ft_sp = common * ts * tp;

    double MT[9];
    mat3_zero(MT);
    double Q_kernel = (q_convention == 1) ? 0.5 * (ft_p - ft_s)
                                          : 0.5 * (ft_s - ft_p);
    MT[0*3+0] = 0.5 * (ft_s + ft_p);
    MT[0*3+1] = Q_kernel;
    MT[1*3+0] = Q_kernel;
    MT[1*3+1] = 0.5 * (ft_s + ft_p);
    MT[2*3+2] = ft_sp;

    /* Stokes basis rotations (same Mishchenko form as reflection branch),
     * scattering cosine cT = cosPsi between p1 (in) and p2 (out). */
    double sT = sqrt(fmax(0.0, 1.0 - cosPsi * cosPsi));
    double denom = fmax(1e-12, sT);
    double s1_safe = fmax(1e-12, s1);
    double s2_safe = fmax(1e-12, s2);

    /* 극 한계 (2026-07-24, 정합테스트 v4).  OCRT_SURF_POLE_EPS 정의부의 주석을
     * 함께 볼 것.  극이 아닌 각도의 식은 이전 판과 글자 하나 다르지 않다. */
    double ci1, si1, ci2, si2;
    if (s1 <= OCRT_SURF_POLE_EPS) {
        ci1 = -copysign(1.0, mu1) * cos_phi;
        si1 = sin_phi;
    } else {
        ci1 = (mu2 - mu1 * cosPsi) / (denom * s1_safe);
        si1 = s2 * sin_phi / denom;
        double r1 = sqrt(ci1*ci1 + si1*si1);
        if (r1 < 1e-12) { ci1 = 1.0; si1 = 0.0; } else { ci1 /= r1; si1 /= r1; }
    }

    if (s2 <= OCRT_SURF_POLE_EPS) {
        ci2 = -copysign(1.0, mu2) * cos_phi;
        si2 = sin_phi;
    } else {
        ci2 = (mu1 - mu2 * cosPsi) / (denom * s2_safe);
        si2 = s1 * sin_phi / denom;
        double r2 = sqrt(ci2*ci2 + si2*si2);
        if (r2 < 1e-12) { ci2 = 1.0; si2 = 0.0; } else { ci2 /= r2; si2 /= r2; }
    }

    double L1[9], L2[9], MTL1[9], L2MTL1[9];
    build_rotation_L(ci1, si1, L1);
    build_rotation_L(ci2, si2, L2);
    mat3_mul(MT, L1, MTL1);
    mat3_mul(L2, MTL1, L2MTL1);
    for (int i = 0; i < 9; ++i) T[i] = L2MTL1[i];
}

/* ----------------------------------------------------------------------------
 * Direct sunglint single-bounce TOA reflectance contribution
 *   ρ_glint = R(mu_v, mu_0, Δφ) · F_solar · exp(-tau/μ_0) · exp(-tau/μ_v) / μ_0
 *
 * In the V3 ρ convention (F_sun = π normalized), the contribution to TOA ρ
 * from direct surface reflection of the solar beam is:
 *   ρ_glint[k] = R[k][0] · μ_0 · exp(-tau/μ_0) · exp(-tau/μ_v)
 *
 * (R has units sr^-1 absorbed into BRDF; μ_0 cancellation per V1 vrt_solver
 *  reflectance_from_field.)
 *
 * Used post-hoc when decouple_sunglint=1.
 * ---------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------------
 * Phase B: Cox-Munk Fourier-mode surface kernel
 *
 * Compute m-mode coefficients R^m_kl(μ_o, μ_i) by numerical φ-quadrature.
 * For Stokes basis I,Q ~ cos(mφ), U ~ sin(mφ):
 *
 *   - cos integration for entries: R_II, R_IQ, R_QI, R_QQ, R_UU
 *   - sin integration for entries: R_IU, R_QU, R_UI, R_UQ
 *
 * Normalization (consistent with (2-δ_{m,0}) reconstruction):
 *   m = 0:  R^0_kl = (1/2π) ∫ R_kl(φ) dφ
 *   m > 0:  R^m_kl = (1/π)  ∫ R_kl(φ) · cos(mφ) dφ      (cos entries)
 *           R^m_kl = (1/π)  ∫ R_kl(φ) · sin(mφ) dφ      (sin entries)
 *
 * Uniform φ-quadrature on [0, 2π): φ_p = (p + 0.5) · 2π/N, w_p = 2π/N
 * (midpoint rule — exactly integrates trig polynomials up to order N-1).
 * ---------------------------------------------------------------------------- */
#ifdef OCRT_FAST_KERNELS
/* v1.11-speed S1 (2026-08-25): FKC type + prototype hoisted above the AIR-side
 * kernel entry point so it can serve from the all-m cache too (definition and
 * design notes below at "v1.09-opt S3").  g_fkc_build_mmax_hint caps the
 * batched build at the caller's true m_max (default: full OCRT_FKC_MMAX,
 * i.e. the pre-hint behavior). */
#define OCRT_FKC_MMAX 64
typedef void (*ocrt_fkc_trig_fn)(double, double, double, double,
                                 double, int, double, int, double *);
typedef struct {
    int valid, n_o, n_i, n_phi, sigma_type, q_conv;
    int built_mmax;               /* modes 0..built_mmax present in coef */
    double ws, nw;
    double *mu_o, *mu_i;
    double *coef;                 /* (built_mmax+1) slabs of n_o*n_i*9, m-major */
} ocrt_fkc_t;
static int ocrt_fkc_serve(ocrt_fkc_t *C, uint32_t pcache_operator,
                          ocrt_fkc_trig_fn trig,
                          const double *mu_o, int n_o,
                          const double *mu_i, int n_i,
                          int m, int n_phi_quad,
                          double ws, int sigma_type,
                          double n_water, int q_convention,
                          double *out_m);
static _Thread_local int g_fkc_build_mmax_hint = OCRT_FKC_MMAX;
#endif

/* Public: declare the true mode-loop upper bound before a Fourier-kernel mode
 * loop, so the batched all-m cache build stops there instead of projecting all
 * OCRT_FKC_MMAX+1 modes.  Values outside [0, OCRT_FKC_MMAX] restore the
 * uncapped default.  Thread-local; call before EACH kernel-building m loop
 * (atm boundary, R_ww, T_wa) with that loop's m_max.  Bit-transparent for all
 * modes actually built. */
void surface_fkc_set_build_mmax(int m_max) {
#ifdef OCRT_FAST_KERNELS
    g_fkc_build_mmax_hint = (m_max >= 0 && m_max <= OCRT_FKC_MMAX)
                                ? m_max : OCRT_FKC_MMAX;
#else
    (void)m_max;
#endif
}

static int surface_coxmunk_fourier_kernel_uncached_(const double *mu_o, int n_o,
                                                      const double *mu_i, int n_i,
                                                      int m, int n_phi_quad,
                                                      double ws, int sigma_type,
                                                      double n_water, int q_convention,
                                                      double *R_m) {
    if (!mu_o || !mu_i || !R_m) return -1;
    if (n_o < 1 || n_i < 1 || m < 0 || n_phi_quad < 4) return -1;

    /* Fourier normalization: m=0 → 1/(2π); m>0 → 1/π */
    const double norm = (m == 0) ? (1.0 / (2.0 * M_PI)) : (1.0 / M_PI);

    /* Which Mueller entries use cos vs sin Fourier?
     *   cos entries: indices (k*3 + l) where (k,l) ∈
     *     (II=0, IQ=1, QI=3, QQ=4, UU=8)
     *   sin entries: (IU=2, QU=5, UI=6, UQ=7)
     *
     * For m = 0:  sin(0)=0 → all sin-entry m=0 coeffs are zero.
     * cos integration with cos(0)=1 → just average.  */
    static const int cos_entry[9] = {1,1,0, 1,1,0, 0,0,1};

    /* Uniform midpoint φ-quadrature on [0, 2π) */
    const double dphi = 2.0 * M_PI / (double)n_phi_quad;

    /* Initialize output to zero */
    memset(R_m, 0, (size_t)n_o * (size_t)n_i * 9 * sizeof(double));

    for (int j_o = 0; j_o < n_o; j_o++) {
        const double mo = mu_o[j_o];
        for (int j_i = 0; j_i < n_i; j_i++) {
            const double mi = mu_i[j_i];
            double *R_pair = R_m + ((size_t)j_o * (size_t)n_i + (size_t)j_i) * 9;

            /* Numerical φ-integration */
            for (int p = 0; p < n_phi_quad; p++) {
                const double phi = (p + 0.5) * dphi;
                const double cphi = cos(phi);
                const double sphi = sin(phi);
                const double cmphi = cos(m * phi);
                const double smphi = sin(m * phi);

                double R_local[9];
                surface_R_coxmunk_trig(mo, mi, cphi, sphi,
                                        ws, sigma_type, n_water,
                                        q_convention, R_local);

                /* Accumulate cos- or sin-weighted entries */
                for (int kl = 0; kl < 9; kl++) {
                    double weight = cos_entry[kl] ? cmphi : smphi;
                    R_pair[kl] += R_local[kl] * weight;
                }
            }
            /* Apply quadrature weight (dphi) and Fourier normalization */
            for (int kl = 0; kl < 9; kl++) {
                R_pair[kl] *= dphi * norm;
            }

            /* OSOAA surface RAA files use the opposite azimuthal sine basis
             * for the incoming-U column entries of the air-side Cox-Munk
             * reflection Fourier matrix.  Keep the reflected-U row entries
             * (UI,UQ) unchanged: direct RAA matrix parity tests show only
             * R13/R23 (IU/QU; indices 2/5) require this sign conversion for
             * m>0. */
            if (m > 0) {
                R_pair[2] = -R_pair[2];
                R_pair[5] = -R_pair[5];
            }
        }
    }
    return 0;
}

/* Exact, bounded per-thread cache for the air-side rough-Fresnel Fourier
 * operator.  A coupled ocean solve requests the same grid-grid and grid-solar
 * slabs in atmosphere pass 1 and pass 2.  Re-integrating the Cox-Munk Mueller
 * kernel over azimuth dominates clean/high-turbidity runs even though the
 * result depends only on the immutable angular grid and surface options.
 *
 * The key is verified by byte comparison of every scalar and both input arrays;
 * the cache therefore cannot return a false hit from a hash collision.  Slots
 * own their copies and are bounded/replaced LRU-style.  Thread-local storage
 * avoids locks in the OpenMP case scheduler.  OCRT_AIR_SURFACE_CACHE_OFF=1 is
 * retained as an authoritative no-cache regression path. */
#define OCRT_AIR_SURFACE_CACHE_SLOTS 128

typedef struct {
    int valid;
    int n_o, n_i, m, n_phi_quad, sigma_type, q_convention;
    double ws, n_water;
    double *mu_o;
    double *mu_i;
    double *kernel;
    size_t kernel_count;
    unsigned long long stamp;
} ocrt_air_surface_cache_entry_t;

static _Thread_local ocrt_air_surface_cache_entry_t
    g_air_surface_cache[OCRT_AIR_SURFACE_CACHE_SLOTS];
static _Thread_local unsigned long long g_air_surface_cache_tick;
static _Thread_local unsigned long long g_air_surface_cache_hits;
static _Thread_local unsigned long long g_air_surface_cache_builds;

static int ocrt_air_surface_cache_enabled_(void) {
    /* Startup option, cached per OpenMP worker to avoid repeated getenv locks. */
    static _Thread_local int enabled = -1;
    if (enabled < 0) {
        const char *off = getenv("OCRT_AIR_SURFACE_CACHE_OFF");
        enabled = !(off && off[0] && strcmp(off, "0") != 0 &&
                    strcmp(off, "off") != 0 && strcmp(off, "false") != 0);
    }
    return enabled;
}

static int ocrt_double_bits_equal_(double a, double b) {
    return memcmp(&a, &b, sizeof(double)) == 0;
}

static int ocrt_air_surface_cache_match_(
    const ocrt_air_surface_cache_entry_t *entry,
    const double *mu_o, int n_o, const double *mu_i, int n_i,
    int m, int n_phi_quad, double ws, int sigma_type,
    double n_water, int q_convention)
{
    if (!entry || !entry->valid || !entry->mu_o || !entry->mu_i ||
        !entry->kernel) return 0;
    if (entry->n_o != n_o || entry->n_i != n_i || entry->m != m ||
        entry->n_phi_quad != n_phi_quad || entry->sigma_type != sigma_type ||
        entry->q_convention != q_convention ||
        !ocrt_double_bits_equal_(entry->ws, ws) ||
        !ocrt_double_bits_equal_(entry->n_water, n_water)) return 0;
    return memcmp(entry->mu_o, mu_o, (size_t)n_o * sizeof(double)) == 0 &&
           memcmp(entry->mu_i, mu_i, (size_t)n_i * sizeof(double)) == 0;
}

static void ocrt_air_surface_cache_release_(
    ocrt_air_surface_cache_entry_t *entry)
{
    if (!entry) return;
    free(entry->mu_o);
    free(entry->mu_i);
    free(entry->kernel);
    memset(entry, 0, sizeof(*entry));
}

static void ocrt_air_surface_cache_diag_(void) {
    const char *diag = getenv("OCRT_AIR_SURFACE_CACHE_DIAG");
    if (diag && diag[0] && strcmp(diag, "0") != 0) {
        fprintf(stderr,
                "[air-surface-cache] hits=%llu builds=%llu slots=%d\n",
                g_air_surface_cache_hits, g_air_surface_cache_builds,
                OCRT_AIR_SURFACE_CACHE_SLOTS);
    }
}

int surface_coxmunk_fourier_kernel(const double *mu_o, int n_o,
                                    const double *mu_i, int n_i,
                                    int m, int n_phi_quad,
                                    double ws, int sigma_type,
                                    double n_water, int q_convention,
                                    double *R_m) {
    if (!mu_o || !mu_i || !R_m) return -1;
    if (n_o < 1 || n_i < 1 || m < 0 || n_phi_quad < 4) return -1;

    const size_t kernel_count = (size_t)n_o * (size_t)n_i * 9u;
    const size_t kernel_bytes = kernel_count * sizeof(double);

#ifdef OCRT_FAST_KERNELS
    /* v1.11-speed S1 (2026-08-25): serve the AIR-side rough-reflection kernel
     * from the all-m FKC instead of re-integrating the m-independent Cox-Munk
     * Mueller kernel once per Fourier mode (measured: 43% of a fair-profile
     * coupled fullgrid run, 162M trig calls at n_mu=48/nphi=1024/m_max=32).
     * Bit-transparent: the FKC per-m accumulation arithmetic is exactly the
     * uncached_ loop's, and the identical m>0 IU/QU sign fold is applied
     * below.  The atm mode loop alternates grid-grid and grid-solar keys, so
     * a small LRU slot array (not a single slot) prevents rebuild thrash.
     * Kill switch for regression comparison: OCRT_AIR_FKC_OFF=1. */
    if (!getenv("OCRT_AIR_FKC_OFF")) {
        enum { OCRT_FKC_AIR_SLOTS = 8 };
        static _Thread_local ocrt_fkc_t g_fkc_air[OCRT_FKC_AIR_SLOTS];
        static _Thread_local unsigned long long air_stamp[OCRT_FKC_AIR_SLOTS];
        static _Thread_local unsigned long long air_tick;
        ++air_tick;
        int slot = -1;
        for (int s = 0; s < OCRT_FKC_AIR_SLOTS; ++s) {
            const ocrt_fkc_t *C = &g_fkc_air[s];
            if (C->valid && C->n_o == n_o && C->n_i == n_i &&
                C->n_phi == n_phi_quad && C->sigma_type == sigma_type &&
                C->q_conv == q_convention && C->ws == ws && C->nw == n_water &&
                memcmp(C->mu_o, mu_o, (size_t)n_o * sizeof(double)) == 0 &&
                memcmp(C->mu_i, mu_i, (size_t)n_i * sizeof(double)) == 0) {
                slot = s;
                break;
            }
        }
        if (slot < 0) {
            unsigned long long oldest = ~0ULL;
            slot = 0;
            for (int s = 0; s < OCRT_FKC_AIR_SLOTS; ++s) {
                if (!g_fkc_air[s].valid) { slot = s; break; }
                if (air_stamp[s] < oldest) { oldest = air_stamp[s]; slot = s; }
            }
        }
        const int fkc_rc = ocrt_fkc_serve(
                           &g_fkc_air[slot], OCRT_SURFACE_PCACHE_OP_R_AA_RAW_ALLM,
                           surface_R_coxmunk_trig,
                           mu_o, n_o, mu_i, n_i, m, n_phi_quad,
                           ws, sigma_type, n_water, q_convention, R_m);
        if (fkc_rc < 0) return -3;
        if (fkc_rc > 0) {
            air_stamp[slot] = air_tick;
            if (m > 0) {
                /* identical to the uncached_ tail: OSOAA-basis IU/QU (indices
                 * 2/5) sign conversion for m>0. */
                const size_t np = (size_t)n_o * (size_t)n_i;
                for (size_t q = 0; q < np; ++q) {
                    R_m[q * 9 + 2] = -R_m[q * 9 + 2];
                    R_m[q * 9 + 5] = -R_m[q * 9 + 5];
                }
            }
            return 0;
        }
        /* serve declined (mode above build cap, or alloc failure):
         * fall through to the legacy per-m cache + uncached_ path. */
    }
#endif

    if (!ocrt_air_surface_cache_enabled_()) {
        return surface_coxmunk_fourier_kernel_uncached_(
            mu_o, n_o, mu_i, n_i, m, n_phi_quad, ws, sigma_type,
            n_water, q_convention, R_m);
    }

    ++g_air_surface_cache_tick;
    for (int slot = 0; slot < OCRT_AIR_SURFACE_CACHE_SLOTS; ++slot) {
        ocrt_air_surface_cache_entry_t *entry = &g_air_surface_cache[slot];
        if (ocrt_air_surface_cache_match_(entry, mu_o, n_o, mu_i, n_i,
                                          m, n_phi_quad, ws, sigma_type,
                                          n_water, q_convention)) {
            memcpy(R_m, entry->kernel, kernel_bytes);
            entry->stamp = g_air_surface_cache_tick;
            ++g_air_surface_cache_hits;
            ocrt_air_surface_cache_diag_();
            return 0;
        }
    }

    const int status = surface_coxmunk_fourier_kernel_uncached_(
        mu_o, n_o, mu_i, n_i, m, n_phi_quad, ws, sigma_type,
        n_water, q_convention, R_m);
    if (status != 0) return status;

    int victim = -1;
    unsigned long long oldest = ~0ULL;
    for (int slot = 0; slot < OCRT_AIR_SURFACE_CACHE_SLOTS; ++slot) {
        const ocrt_air_surface_cache_entry_t *entry = &g_air_surface_cache[slot];
        if (!entry->valid) { victim = slot; break; }
        if (entry->stamp < oldest) { oldest = entry->stamp; victim = slot; }
    }
    if (victim >= 0) {
        ocrt_air_surface_cache_entry_t *entry = &g_air_surface_cache[victim];
        ocrt_air_surface_cache_release_(entry);
        entry->mu_o = (double *)malloc((size_t)n_o * sizeof(double));
        entry->mu_i = (double *)malloc((size_t)n_i * sizeof(double));
        entry->kernel = (double *)malloc(kernel_bytes);
        if (entry->mu_o && entry->mu_i && entry->kernel) {
            memcpy(entry->mu_o, mu_o, (size_t)n_o * sizeof(double));
            memcpy(entry->mu_i, mu_i, (size_t)n_i * sizeof(double));
            memcpy(entry->kernel, R_m, kernel_bytes);
            entry->n_o = n_o;
            entry->n_i = n_i;
            entry->m = m;
            entry->n_phi_quad = n_phi_quad;
            entry->sigma_type = sigma_type;
            entry->q_convention = q_convention;
            entry->ws = ws;
            entry->n_water = n_water;
            entry->kernel_count = kernel_count;
            entry->stamp = g_air_surface_cache_tick;
            entry->valid = 1;
            ++g_air_surface_cache_builds;
        } else {
            ocrt_air_surface_cache_release_(entry);
        }
    }
    ocrt_air_surface_cache_diag_();
    return 0;
}

/* ----------------------------------------------------------------------------
 * surface_R_ww_coxmunk_trig — WATER-SIDE internal reflection Mueller (3x3) at
 * one (mu_out, mu_in_dn, phi) triple. z-flip analogue of surface_R_coxmunk_trig
 * (the external/air-side kernel): identical Cox-Munk slope geometry, Hovenier
 * rotations, Sancer shadowing, and Wang-Gordon BRDF, with ONLY the facet
 * Fresnel changed to the water->air (internal, TIR) reflection used by A2.
 * Semantics here: mu_out = |mu| of the reflected DOWNWARD ray (into water),
 * mu_in_dn = |mu| of the incident UPWARD ray (the in-water upwelling).
 * ---------------------------------------------------------------------------- */
void surface_R_ww_coxmunk_trig(double mu_out, double mu_in_dn,
                             double cos_phi, double sin_phi,
                             double ws, int sigma_type_int, double n_water,
                             int q_convention, double *R) {
    mat3_zero(R);

    double sigma_sq = surface_slope_variance(ws, sigma_type_int);

    double mu_in_actual = -mu_in_dn;
    double s1 = sqrt(fmax(0.0, 1.0 - mu_out * mu_out));
    double s2 = sqrt(fmax(0.0, 1.0 - mu_in_actual * mu_in_actual));

    double cT = mu_out * mu_in_actual + s1 * s2 * cos_phi;
    if      (cT > +1.0) cT = +1.0;
    else if (cT < -1.0) cT = -1.0;
    double cos_omega = sqrt(fmax(0.0, 0.5 * (1.0 - cT)));
    if (cos_omega <= 1e-6) return;
    double cos_beta = (mu_out + mu_in_dn) / (2.0 * cos_omega);
    if (cos_beta <= 1e-6) return;

    double sT = sqrt(fmax(0.0, 1.0 - cT * cT));
    double denom = fmax(1e-12, sT);

    /* Stokes rotation angles (Refactoring Tasks Task 2 fix, 2026-05-10):
     *   표준 (Mishchenko 2002 Appendix B):
     *     cos σ_1 = (μ_out - μ_in·cos Θ) / (sin Θ · sin θ_in)
     *     cos σ_2 = (μ_in - μ_out·cos Θ) / (sin Θ · sin θ_out)
     *   기존 V3 식 가 분모 에 sin θ 누락 — typical case s ≈ 1 에서 small,
     *   asymmetric geometry (sza ≠ vza) 에서 큰 deviation 발생 가능.
     *   Fix: ci1 분모 에 s2 (= sin θ_in), ci2 분모 에 s1 (= sin θ_out) 추가. */
    double s2_safe = fmax(1e-12, s2);
    double s1_safe = fmax(1e-12, s1);

    /* 극 한계 (2026-07-24, 정합테스트 v4).  OCRT_SURF_POLE_EPS 정의부의 주석을
     * 함께 볼 것.  극이 아닌 각도의 식은 이전 판과 글자 하나 다르지 않다. */
    double ci1, si1, ci2, si2;
    if (s2 <= OCRT_SURF_POLE_EPS) {
        ci1 = -copysign(1.0, mu_in_actual) * cos_phi;
        si1 = sin_phi;
    } else {
        ci1 = (mu_out - mu_in_actual * cT) / (denom * s2_safe);
        si1 = s1 * sin_phi / denom;
        double n1 = sqrt(ci1*ci1 + si1*si1);
        if (n1 < 1e-12) { ci1 = 1.0; si1 = 0.0; }
        else { ci1 /= n1; si1 /= n1; }
    }

    if (s1 <= OCRT_SURF_POLE_EPS) {
        ci2 = -copysign(1.0, mu_out) * cos_phi;
        si2 = sin_phi;
    } else {
        ci2 = (mu_in_actual - mu_out * cT) / (denom * s1_safe);
        si2 = s2 * sin_phi / denom;
        double n2 = sqrt(ci2*ci2 + si2*si2);
        if (n2 < 1e-12) { ci2 = 1.0; si2 = 0.0; }
        else { ci2 /= n2; si2 /= n2; }
    }

    /* Facet Fresnel Mueller at microfacet incidence omega — WATER-SIDE
     * internal reflection: water(n_water) -> air(1.0), including TIR.
     * Reuses surface_flat_fresnel_R_matrix_general, the *same* routine the
     * flat M_Rww (rt_air_water_R_ww, A2) uses, so the slope-integrated rough
     * kernel reduces to the validated flat per-node limit as wind -> 0. */
    double MF[9];
    surface_flat_fresnel_R_matrix_general(cos_omega, n_water, 1.0, q_convention, MF);

    /* Cox-Munk isotropic Gaussian slope */
    double cos_beta_sq = cos_beta * cos_beta;
    double tan_beta_sq = (1.0 - cos_beta_sq) / fmax(cos_beta_sq, 1e-8);
    double P_slope = (1.0 / (M_PI * sigma_sq)) * exp(-tan_beta_sq / sigma_sq);

    /* Sancer (1969) bistatic shadowing for an isotropic Gaussian slope
     * distribution.  OCRT uses the additive two-direction combination
     * S = 1/(1 + lambda_i + lambda_v). */
    double sigma_len = sqrt(sigma_sq);
    double lam_i = 0.0, lam_v = 0.0;

    if (mu_in_dn < 1.0 - 1e-10) {
        double sin_i = sqrt(fmax(1e-20, 1.0 - mu_in_dn * mu_in_dn));
        double nu_i = (mu_in_dn / sin_i) / sigma_len;
        if (nu_i <= 20.0) {
            lam_i = 0.5 * (exp(-nu_i * nu_i) / (nu_i * sqrt(M_PI)) - erfc(nu_i));
        }
    }
    if (mu_out < 1.0 - 1e-10) {
        double sin_v = sqrt(fmax(1e-20, 1.0 - mu_out * mu_out));
        double nu_v = (mu_out / sin_v) / sigma_len;
        if (nu_v <= 20.0) {
            lam_v = 0.5 * (exp(-nu_v * nu_v) / (nu_v * sqrt(M_PI)) - erfc(nu_v));
        }
    }
    double S_bi = 1.0 / (1.0 + lam_i + lam_v);

    /* Cox-Munk/Wang-Gordon radiance BRDF. */
    double brdf = (P_slope * S_bi) / (4.0 * mu_in_dn * mu_out * cos_beta_sq * cos_beta_sq);

    /* Assemble: R = L2 · MF · L1 · brdf */
    double L1[9], L2[9], MFL1[9], L2MFL1[9];
    build_rotation_L(ci1, si1, L1);
    build_rotation_L(ci2, si2, L2);
    mat3_mul(MF, L1, MFL1);
    mat3_mul(L2, MFL1, L2MFL1);
    for (int i = 0; i < 9; ++i) R[i] = L2MFL1[i] * brdf;
}


/* ----------------------------------------------------------------------------
 * surface_R_ww_coxmunk_fourier_kernel — m-mode Fourier coefficients of the
 * WATER-SIDE internal-reflection Mueller (surface_R_ww_coxmunk_trig). Identical
 * Fourier convention to surface_coxmunk_fourier_kernel (norm 1/2pi for m=0,
 * 1/pi for m>0; cos/sin entry split). Output R_m[n_o*n_i*9], packed
 * R_m[(j_o*n_i + j_i)*9 + (k*3+l)]. mu_o = downward (reflected) |mu|,
 * mu_i = upward (incident) |mu|.
 * ---------------------------------------------------------------------------- */
#ifdef OCRT_FAST_KERNELS
/* ---- v1.09-opt S3: all-m Fourier-kernel cache (OSOAA SURF_MATR analog) ----
 * The per-azimuth Mueller M(mu_o,mu_i,phi) is m-INDEPENDENT; the original
 * kernels recompute it for every Fourier mode (34x per solve -> 19M trig
 * calls, the measured layer-independent floor).  First call for a geometry
 * does ONE phi pass projecting onto ALL modes 0..OCRT_FKC_MMAX; later m
 * calls memcpy their slab.  BIT-SAFETY: weights use the same cos(m*phi)/
 * sin(m*phi) libm expressions and the p-accumulation order per (m,pair,kl)
 * matches the original loop exactly, so served coefficients are the exact
 * bits the direct evaluation produces.  m > MMAX or alloc failure falls
 * through to the original body (always correct).  Thread-local.
 * (v1.11-speed S1: type/prototype hoisted above the air kernel entry; build
 * capped at surface_fkc_set_build_mmax; built_mmax records the built range.) */
static int ocrt_fkc_serve(ocrt_fkc_t *C, uint32_t pcache_operator,
                          ocrt_fkc_trig_fn trig,
                          const double *mu_o, int n_o,
                          const double *mu_i, int n_i,
                          int m, int n_phi_quad,
                          double ws, int sigma_type,
                          double n_water, int q_convention,
                          double *out_m)
{
    if (m > OCRT_FKC_MMAX) return 0;
    const size_t pair9 = (size_t)n_o * (size_t)n_i * 9u;
    int keymatch = C->valid && C->n_o == n_o && C->n_i == n_i &&
              C->n_phi == n_phi_quad && C->sigma_type == sigma_type &&
              C->q_conv == q_convention && C->ws == ws && C->nw == n_water &&
              memcmp(C->mu_o, mu_o, (size_t)n_o * sizeof(double)) == 0 &&
              memcmp(C->mu_i, mu_i, (size_t)n_i * sizeof(double)) == 0;
    if (keymatch && m <= C->built_mmax) {
        memcpy(out_m, C->coef + (size_t)m * pair9, pair9 * sizeof(double));
        return 1;
    }

    /* Build only the caller-declared mode range.  The persistent representation
     * stores the PRE-fold all-mode slab; public R_aa/T_wa wrappers apply their
     * existing incoming-U convention fold exactly once after serving. */
    const int bmax = g_fkc_build_mmax_hint;
    if (m > bmax) return 0;
    const size_t all_count = (size_t)(bmax + 1) * pair9;

    free(C->mu_o); free(C->mu_i); free(C->coef);
    C->valid = 0; C->built_mmax = -1;
    C->mu_o = (double*)malloc((size_t)n_o * sizeof(double));
    C->mu_i = (double*)malloc((size_t)n_i * sizeof(double));
    C->coef = (double*)malloc(all_count * sizeof(double));
    if (!C->mu_o || !C->mu_i || !C->coef) {
        free(C->mu_o); free(C->mu_i); free(C->coef);
        C->mu_o = C->mu_i = NULL; C->coef = NULL; C->valid = 0;
        return ocrt_surface_pcache_mode() == OCRT_SURFACE_PCACHE_BUILD ? -1 : 0;
    }
    memcpy(C->mu_o, mu_o, (size_t)n_o * sizeof(double));
    memcpy(C->mu_i, mu_i, (size_t)n_i * sizeof(double));

    ocrt_surface_pcache_key_t pkey;
    memset(&pkey, 0, sizeof pkey);
    pkey.operator_kind = pcache_operator;
    pkey.mode_first = 0;
    pkey.mode_count = bmax + 1;
    pkey.n_o = n_o; pkey.n_i = n_i; pkey.n_phi_quad = n_phi_quad;
    pkey.sigma_type = sigma_type; pkey.q_convention = q_convention;
    pkey.wind_speed = ws; pkey.n_water = n_water;
    pkey.mu_o = mu_o; pkey.mu_i = mu_i;

    const int phit = ocrt_surface_pcache_load(&pkey, C->coef, all_count);
    if (phit < 0) {
        free(C->mu_o); free(C->mu_i); free(C->coef);
        C->mu_o = C->mu_i = NULL; C->coef = NULL; C->valid = 0;
        return -1;
    }
    if (phit > 0) {
        C->n_o = n_o; C->n_i = n_i; C->n_phi = n_phi_quad;
        C->sigma_type = sigma_type; C->q_conv = q_convention;
        C->ws = ws; C->nw = n_water; C->built_mmax = bmax; C->valid = 1;
        memcpy(out_m, C->coef + (size_t)m * pair9, pair9 * sizeof(double));
        return 1;
    }

    double *W = (double*)malloc((size_t)2 * (size_t)(bmax + 1)
                                * (size_t)n_phi_quad * sizeof(double));
    if (!W) {
        free(C->mu_o); free(C->mu_i); free(C->coef);
        C->mu_o = C->mu_i = NULL; C->coef = NULL; C->valid = 0;
        return ocrt_surface_pcache_mode() == OCRT_SURFACE_PCACHE_BUILD ? -1 : 0;
    }
    const double dphi = 2.0 * M_PI / (double)n_phi_quad;
    for (int mm = 0; mm <= bmax; mm++)
        for (int p = 0; p < n_phi_quad; p++) {
            const double phi = (p + 0.5) * dphi;
            W[((size_t)2*mm    ) * (size_t)n_phi_quad + p] = cos(mm * phi);
            W[((size_t)2*mm + 1) * (size_t)n_phi_quad + p] = sin(mm * phi);
        }
    static const int cos_entry[9] = {1,1,0, 1,1,0, 0,0,1};
    for (int j_o = 0; j_o < n_o; j_o++) {
        const double mo = mu_o[j_o];
        for (int j_i = 0; j_i < n_i; j_i++) {
            const double mi = mu_i[j_i];
            const size_t pair = (size_t)j_o * (size_t)n_i + (size_t)j_i;
            double acc[OCRT_FKC_MMAX + 1][9];
            memset(acc, 0, sizeof(acc));
            for (int p = 0; p < n_phi_quad; p++) {
                const double phi = (p + 0.5) * dphi;
                const double cphi = cos(phi);
                const double sphi = sin(phi);
                double R_local[9];
                trig(mo, mi, cphi, sphi, ws, sigma_type,
                     n_water, q_convention, R_local);
                for (int mm = 0; mm <= bmax; mm++) {
                    const double cm = W[((size_t)2*mm    ) * (size_t)n_phi_quad + p];
                    const double sm = W[((size_t)2*mm + 1) * (size_t)n_phi_quad + p];
                    for (int kl = 0; kl < 9; kl++)
                        acc[mm][kl] += R_local[kl] * (cos_entry[kl] ? cm : sm);
                }
            }
            for (int mm = 0; mm <= bmax; mm++) {
                const double norm_mm = (mm == 0) ? (1.0/(2.0*M_PI)) : (1.0/M_PI);
                double *dst = C->coef + (size_t)mm * pair9 + pair * 9;
                for (int kl = 0; kl < 9; kl++)
                    dst[kl] = acc[mm][kl] * (dphi * norm_mm);
            }
        }
    }
    free(W);
    C->n_o = n_o; C->n_i = n_i; C->n_phi = n_phi_quad;
    C->sigma_type = sigma_type; C->q_conv = q_convention;
    C->ws = ws; C->nw = n_water; C->built_mmax = bmax; C->valid = 1;

    if (ocrt_surface_pcache_store(&pkey, C->coef, all_count) != 0)
        return -1;
    memcpy(out_m, C->coef + (size_t)m * pair9, pair9 * sizeof(double));
    return 1;
}
static _Thread_local ocrt_fkc_t g_fkc_rww;
static _Thread_local ocrt_fkc_t g_fkc_twa;
#endif

int surface_R_ww_coxmunk_fourier_kernel(const double *mu_o, int n_o,
                                    const double *mu_i, int n_i,
                                    int m, int n_phi_quad,
                                    double ws, int sigma_type,
                                    double n_water, int q_convention,
                                    double *R_m) {
    if (!mu_o || !mu_i || !R_m) return -1;
    if (n_o < 1 || n_i < 1 || m < 0 || n_phi_quad < 4) return -1;
#ifdef OCRT_FAST_KERNELS
    {
        const int fkc_rc = ocrt_fkc_serve(
                           &g_fkc_rww, OCRT_SURFACE_PCACHE_OP_R_WW_RAW_ALLM,
                           surface_R_ww_coxmunk_trig,
                           mu_o, n_o, mu_i, n_i, m, n_phi_quad,
                           ws, sigma_type, n_water, q_convention, R_m);
        if (fkc_rc < 0) return -3;
        if (fkc_rc > 0) return 0;
    }
#endif

    /* Fourier normalization: m=0 → 1/(2π); m>0 → 1/π */
    const double norm = (m == 0) ? (1.0 / (2.0 * M_PI)) : (1.0 / M_PI);

    /* Which Mueller entries use cos vs sin Fourier?
     *   cos entries: indices (k*3 + l) where (k,l) ∈
     *     (II=0, IQ=1, QI=3, QQ=4, UU=8)
     *   sin entries: (IU=2, QU=5, UI=6, UQ=7)
     *
     * For m = 0:  sin(0)=0 → all sin-entry m=0 coeffs are zero.
     * cos integration with cos(0)=1 → just average.  */
    static const int cos_entry[9] = {1,1,0, 1,1,0, 0,0,1};

    /* Uniform midpoint φ-quadrature on [0, 2π) */
    const double dphi = 2.0 * M_PI / (double)n_phi_quad;

    /* Initialize output to zero */
    memset(R_m, 0, (size_t)n_o * (size_t)n_i * 9 * sizeof(double));

    for (int j_o = 0; j_o < n_o; j_o++) {
        const double mo = mu_o[j_o];
        for (int j_i = 0; j_i < n_i; j_i++) {
            const double mi = mu_i[j_i];
            double *R_pair = R_m + ((size_t)j_o * (size_t)n_i + (size_t)j_i) * 9;

            /* Numerical φ-integration */
            for (int p = 0; p < n_phi_quad; p++) {
                const double phi = (p + 0.5) * dphi;
                const double cphi = cos(phi);
                const double sphi = sin(phi);
                const double cmphi = cos(m * phi);
                const double smphi = sin(m * phi);

                double R_local[9];
                surface_R_ww_coxmunk_trig(mo, mi, cphi, sphi,
                                        ws, sigma_type, n_water,
                                        q_convention, R_local);

                /* Accumulate cos- or sin-weighted entries */
                for (int kl = 0; kl < 9; kl++) {
                    double weight = cos_entry[kl] ? cmphi : smphi;
                    R_pair[kl] += R_local[kl] * weight;
                }
            }
            /* Apply quadrature weight (dphi) and Fourier normalization */
            for (int kl = 0; kl < 9; kl++) {
                R_pair[kl] *= dphi * norm;
            }
        }
    }
    return 0;
}

/* WATER->AIR TRANSMISSION Fourier-mode Mueller kernel (Milestone 2b, 2026-06-04).
 * Azimuth-mode (m) coefficients T^m_kl(mu_o_air, mu_i_water) of the Cox-Munk
 * water->air BTDF, obtained by phi-quadrature integration of the per-azimuth
 * transmission Mueller surface_T_coxmunk_trig. This is the exact transmission
 * analog of surface_R_ww_coxmunk_fourier_kernel: identical Fourier normalization
 * (m=0 -> 1/2pi, m>0 -> 1/pi) and identical cos/sin Mueller-entry layout; the
 * ONLY change is the per-azimuth kernel (reflection surface_R_ww_coxmunk_trig ->
 * transmission surface_T_coxmunk_trig, which already carries the n^2 radiance
 * factor and returns the zero matrix on TIR).
 *   mu_o[n_o] : AIR-side outgoing |mu| (>0)
 *   mu_i[n_i] : in-water upwelling |mu| (>0)
 *   T_m[n_o*n_i*9] : packed row-major  T_m[(j_o*n_i + j_i)*9 + (k*3+l)]
 * Used by the rough (wind>0) reverse coupling rt_air_water_couple_water_to_atm.
 * Returns 0 on success. */
/* Polarized air-to-water facet transmission, without shadowing.  The shared
 * microfacet geometry above supplies the same half vector and slope density as
 * the scalar BTDF.  Fresnel transmission is then split into s/p components and
 * rotated between the incident/outgoing meridian bases. */
static void surface_T_aw_coxmunk_trig(double mu_i_air, double mu_o_water,
                                      double cphi, double sphi,
                                      double ws, int sigma_type,
                                      double n_water, int q_convention,
                                      double *T) {
    mat3_zero(T);
    double s2v = surface_slope_variance(ws, sigma_type);
    double mu_i = mu_i_air, mu_o = mu_o_water;
    double n = n_water;
    surface_aw_microfacet_geometry_t g;
    if (surface_aw_microfacet_geometry_(mu_i, mu_o, cphi, sphi,
                                        n, s2v, &g) != 0) return;
    double si = g.sin_i, so = g.sin_o;
    /* air(1)->water(n) Fresnel at local incidence */
    double cti = g.incidence_cosine;
    double sti = sqrt(fmax(0.0, 1.0 - cti*cti));
    double stt = sti / n;
    double Ts, Tp;
    if (stt >= 1.0) { Ts = 0.0; Tp = 0.0; }
    else {
        double ctt = sqrt(fmax(0.0, 1.0 - stt*stt));
        double rs = (cti - n*ctt) / (cti + n*ctt);
        double rp = (n*cti - ctt) / (n*cti + ctt);
        Ts = 1.0 - rs*rs;
        Tp = 1.0 - rp*rp;
    }
    double denom_h = g.half_denom_sq;
    if (denom_h < 1e-12) return;
    double common = (fabs(g.incidence_cosine) * fabs(g.exit_dot_half)) /
                    (mu_i * mu_o)
                  * (n * n * g.slope_density) / denom_h;
    double ft_s  = common * Ts;
    double ft_p  = common * Tp;
    double ft_sp = common * sqrt(fmax(0.0, Ts * Tp));

    double MT[9];
    mat3_zero(MT);
    double Q_kernel = (q_convention == 1) ? 0.5 * (ft_p - ft_s)
                                          : 0.5 * (ft_s - ft_p);
    MT[0*3+0] = 0.5 * (ft_s + ft_p);
    MT[0*3+1] = Q_kernel;
    MT[1*3+0] = Q_kernel;
    MT[1*3+1] = 0.5 * (ft_s + ft_p);
    MT[2*3+2] = ft_sp;

    /* Mishchenko meridian rotations, roles: incoming = air, outgoing = water.
     * surface_aw_microfacet_geometry_ uses directions pointing away from the
     * interface.  Stokes rotations, however, must use photon propagation
     * directions: the incident air ray is the negative of that geometry
     * vector, while the transmitted water ray is already downward.  Hence
     * both signed direction cosines are negative and the relative azimuth is
     * shifted by pi. */
    const double cphi_r = -cphi;
    const double sphi_r = -sphi;
    const double mu_i_s = -mu_i;
    const double mu_o_s = -mu_o;
    double cosPsi = mu_i_s * mu_o_s + si * so * cphi_r;
    if      (cosPsi > +1.0) cosPsi = +1.0;
    else if (cosPsi < -1.0) cosPsi = -1.0;
    double sT = sqrt(fmax(0.0, 1.0 - cosPsi * cosPsi));
    double dnm = fmax(1e-12, sT);
    double s1_safe = fmax(1e-12, si);
    double s2_safe = fmax(1e-12, so);
    /* 극 한계 (2026-07-24, 정합테스트 v4).  OCRT_SURF_POLE_EPS 정의부의 주석을
     * 함께 볼 것.  극이 아닌 각도의 식은 이전 판과 글자 하나 다르지 않다. */
    double ci1, si1, ci2, si2;
    if (si <= OCRT_SURF_POLE_EPS) {
        ci1 = -copysign(1.0, mu_i_s) * cphi_r;
        si1 = sphi_r;
    } else {
        ci1 = (mu_o_s - mu_i_s * cosPsi) / (dnm * s1_safe);
        si1 = so * sphi_r / dnm;
        double r1 = sqrt(ci1*ci1 + si1*si1);
        if (r1 < 1e-12) { ci1 = 1.0; si1 = 0.0; } else { ci1 /= r1; si1 /= r1; }
    }
    if (so <= OCRT_SURF_POLE_EPS) {
        ci2 = -copysign(1.0, mu_o_s) * cphi_r;
        si2 = sphi_r;
    } else {
        ci2 = (mu_i_s - mu_o_s * cosPsi) / (dnm * s2_safe);
        si2 = si * sphi_r / dnm;
        double r2 = sqrt(ci2*ci2 + si2*si2);
        if (r2 < 1e-12) { ci2 = 1.0; si2 = 0.0; } else { ci2 /= r2; si2 /= r2; }
    }
    double L1[9], L2[9], MTL1[9], OUT[9];
    build_rotation_L(ci1, si1, L1);
    build_rotation_L(ci2, si2, L2);
    mat3_mul(MT, L1, MTL1);
    mat3_mul(L2, MTL1, OUT);
    for (int i = 0; i < 9; ++i) T[i] = OUT[i];
}

/* test export */
void surface_T_aw_coxmunk_trig_test(double a,double b,double c,double d,double e,int f,double g,int h,double*T){
    surface_T_aw_coxmunk_trig(a,b,c,d,e,f,g,h,T);
}
#ifdef OCRT_FAST_KERNELS
static unsigned long long taw_fnv64(const double *v, int n) {
    const unsigned char *p = (const unsigned char*)v;
    unsigned long long h = 1469598103934665603ULL;
    for (size_t i = 0; i < (size_t)n * sizeof(double); ++i) {
        h ^= p[i]; h *= 1099511628211ULL;
    }
    return h;
}
#endif

int surface_T_aw_coxmunk_fourier_kernel(const double *mu_o, int n_o,
                                        const double *mu_i, int n_i,
                                        int m, int n_phi_quad,
                                        double ws, int sigma_type,
                                        double n_water, int q_convention,
                                        double *T_m) {
    if (!mu_o || !mu_i || !T_m) return -1;
    if (n_o < 1 || n_i < 1 || m < 0 || n_phi_quad < 4) return -1;
    if (surface_slope_variance(ws, sigma_type) <= 1e-10) return -2; /* flat */
    const size_t need = (size_t)n_o * (size_t)n_i * 9u;
#ifdef OCRT_FAST_KERNELS
    /* v1.10 B-0c.2a: SELF-memoization - caches this function's OWN output,
     * hence bit-transparent by construction.  ocrt_fkc_serve was tried and
     * REJECTED here: its internal fourier machinery is specific to the
     * R_ww/T_wa kernels and produced a materially wrong T_aw table
     * (II element errors up to ~48 vs lobe peak ~12). */
    static _Thread_local struct {
        int    valid, n_o, n_i, m, nphi, st, qc;
        double ws, nw; unsigned long long ho, hi;
        double *tab; size_t cap;
    } taw_memo[8];
    {
        int sl = m % 8;
        if (taw_memo[sl].valid &&
            taw_memo[sl].n_o == n_o && taw_memo[sl].n_i == n_i &&
            taw_memo[sl].m == m && taw_memo[sl].nphi == n_phi_quad &&
            taw_memo[sl].st == sigma_type && taw_memo[sl].qc == q_convention &&
            taw_memo[sl].ws == ws && taw_memo[sl].nw == n_water &&
            taw_memo[sl].ho == taw_fnv64(mu_o, n_o) &&
            taw_memo[sl].hi == taw_fnv64(mu_i, n_i) &&
            taw_memo[sl].cap >= need) {
            memcpy(T_m, taw_memo[sl].tab, need * sizeof(double));
            return 0;
        }
    }
#endif

    /* Persistent T_aw stores the FINAL per-mode operator, including the
     * incoming-U column fold.  This path deliberately does not reuse the
     * generic all-m R_ww/T_wa builder, preserving the validated T_aw numerics. */
    ocrt_surface_pcache_key_t pkey;
    memset(&pkey, 0, sizeof pkey);
    pkey.operator_kind = OCRT_SURFACE_PCACHE_OP_T_AW_FINAL_M;
    pkey.mode_first = m; pkey.mode_count = 1;
    pkey.n_o = n_o; pkey.n_i = n_i; pkey.n_phi_quad = n_phi_quad;
    pkey.sigma_type = sigma_type; pkey.q_convention = q_convention;
    pkey.wind_speed = ws; pkey.n_water = n_water;
    pkey.mu_o = mu_o; pkey.mu_i = mu_i;
    const int phit = ocrt_surface_pcache_load(&pkey, T_m, need);
    if (phit < 0) return -3;
    if (phit > 0) {
#ifdef OCRT_FAST_KERNELS
        int sl = m % 8;
        if (taw_memo[sl].cap < need) {
            free(taw_memo[sl].tab);
            taw_memo[sl].tab = (double*)malloc(need * sizeof(double));
            taw_memo[sl].cap = taw_memo[sl].tab ? need : 0;
        }
        if (taw_memo[sl].tab) {
            memcpy(taw_memo[sl].tab, T_m, need * sizeof(double));
            taw_memo[sl].valid = 1;
            taw_memo[sl].n_o = n_o; taw_memo[sl].n_i = n_i;
            taw_memo[sl].m = m; taw_memo[sl].nphi = n_phi_quad;
            taw_memo[sl].st = sigma_type; taw_memo[sl].qc = q_convention;
            taw_memo[sl].ws = ws; taw_memo[sl].nw = n_water;
            taw_memo[sl].ho = taw_fnv64(mu_o, n_o);
            taw_memo[sl].hi = taw_fnv64(mu_i, n_i);
        }
#endif
        return 0;
    }

    const double norm = (m == 0) ? (1.0 / (2.0 * M_PI)) : (1.0 / M_PI);
    static const int cos_entry[9] = {1,1,0, 1,1,0, 0,0,1};
    const double dphi = 2.0 * M_PI / (double)n_phi_quad;
    memset(T_m, 0, need * sizeof(double));

    for (int j_o = 0; j_o < n_o; j_o++) {           /* in-water outgoing */
        const double mw = mu_o[j_o];
        for (int j_i = 0; j_i < n_i; j_i++) {       /* air incident */
            const double ma = mu_i[j_i];
            double *T_pair = T_m + ((size_t)j_o * (size_t)n_i + (size_t)j_i) * 9;
            for (int p = 0; p < n_phi_quad; p++) {
                const double phi  = (p + 0.5) * dphi;
                const double cm   = cos(m * phi);
                const double sm   = sin(m * phi);
                double Tl[9];
                surface_T_aw_coxmunk_trig(ma, mw, cos(phi), sin(phi),
                                          ws, sigma_type, n_water,
                                          q_convention, Tl);
                for (int kl = 0; kl < 9; kl++) {
                    double weight = cos_entry[kl] ? cm : sm;
                    T_pair[kl] += Tl[kl] * weight;
                }
            }
            for (int kl = 0; kl < 9; kl++) T_pair[kl] *= norm * dphi;
            /* The OCRT Fourier field stores U as a sine coefficient.  In the
             * relative-azimuth convolution, the incoming-U column acquires one
             * minus sign for m>0.  Fold it into the stored operator once so the
             * downstream 3x3 contraction remains an ordinary matrix-vector
             * product, matching OSOAA Step-7 TAW storage/contraction. */
            if (m > 0) {
                T_pair[0*3 + 2] = -T_pair[0*3 + 2];
                T_pair[1*3 + 2] = -T_pair[1*3 + 2];
            }
        }
    }
#ifdef OCRT_FAST_KERNELS
    {
        int sl = m % 8;
        if (taw_memo[sl].cap < need) {
            free(taw_memo[sl].tab);
            taw_memo[sl].tab = (double*)malloc(need * sizeof(double));
            taw_memo[sl].cap = taw_memo[sl].tab ? need : 0;
        }
        if (taw_memo[sl].tab) {
            memcpy(taw_memo[sl].tab, T_m, need * sizeof(double));
            taw_memo[sl].valid = 1;
            taw_memo[sl].n_o = n_o;  taw_memo[sl].n_i = n_i;
            taw_memo[sl].m = m;      taw_memo[sl].nphi = n_phi_quad;
            taw_memo[sl].st = sigma_type; taw_memo[sl].qc = q_convention;
            taw_memo[sl].ws = ws;    taw_memo[sl].nw = n_water;
            taw_memo[sl].ho = taw_fnv64(mu_o, n_o);
            taw_memo[sl].hi = taw_fnv64(mu_i, n_i);
        }
    }
#endif
    if (ocrt_surface_pcache_store(&pkey, T_m, need) != 0) return -3;
    return 0;
}
/* FIX 2026-08-25 (near-nadir water-leaving Q deficit): m>0 U-column sign fold
 * for the water->air Fourier transmission operator.  The identical fold exists
 * in surface_T_aw_coxmunk_fourier_kernel ("matching OSOAA Step-7 TAW
 * storage/contraction") but was missing here, so the m>=1 polarized modes of
 * the transmitted water-leaving field lost their spin-2 closure: Q was damped
 * (0.57x at VZA<=10, wind 3) while U stayed correct.  Unit gate: synthetic
 * m=2 field through rt_air_water_couple_water_to_atm must give
 * (T_I, T_Q2, T_U2) = rough/flat = (1,1,1) at all near-nadir angles.
 * See docs/validation cross_check_fr631_20260823 Q_TRUNCATION_FINDINGS §8-9. */
static void ocrt_twa_fold_ucol_m(double *T_m, int n_o, int n_i, int m) {
    if (m <= 0) return;
    const size_t np = (size_t)n_o * (size_t)n_i;
    for (size_t q = 0; q < np; ++q) {
        T_m[q * 9 + 2] = -T_m[q * 9 + 2];   /* IU */
        T_m[q * 9 + 5] = -T_m[q * 9 + 5];   /* QU */
    }
}

int surface_T_wa_coxmunk_fourier_kernel(const double *mu_o, int n_o,
                                        const double *mu_i, int n_i,
                                        int m, int n_phi_quad,
                                        double ws, int sigma_type,
                                        double n_water, int q_convention,
                                        double *T_m) {
    if (!mu_o || !mu_i || !T_m) return -1;
    if (n_o < 1 || n_i < 1 || m < 0 || n_phi_quad < 4) return -1;
#ifdef OCRT_FAST_KERNELS
    if (!getenv("OCRT_TWA_FKC_OFF")) {
        const int fkc_rc = ocrt_fkc_serve(
                           &g_fkc_twa, OCRT_SURFACE_PCACHE_OP_T_WA_RAW_ALLM,
                           surface_T_coxmunk_trig,
                           mu_o, n_o, mu_i, n_i, m, n_phi_quad,
                           ws, sigma_type, n_water, q_convention, T_m);
        if (fkc_rc < 0) return -3;
        if (fkc_rc > 0) {
            ocrt_twa_fold_ucol_m(T_m, n_o, n_i, m);
            return 0;
        }
    }
#endif

    /* Fourier normalization: m=0 -> 1/(2pi); m>0 -> 1/pi (same as R_ww kernel) */
    const double norm = (m == 0) ? (1.0 / (2.0 * M_PI)) : (1.0 / M_PI);

    /* cos entries: II,IQ,QI,QQ,UU ; sin entries: IU,QU,UI,UQ (same as R_ww) */
    static const int cos_entry[9] = {1,1,0, 1,1,0, 0,0,1};

    const double dphi = 2.0 * M_PI / (double)n_phi_quad;
    memset(T_m, 0, (size_t)n_o * (size_t)n_i * 9 * sizeof(double));

    for (int j_o = 0; j_o < n_o; j_o++) {
        const double mo = mu_o[j_o];   /* air-side outgoing */
        for (int j_i = 0; j_i < n_i; j_i++) {
            const double mi = mu_i[j_i];   /* in-water upwelling */
            double *T_pair = T_m + ((size_t)j_o * (size_t)n_i + (size_t)j_i) * 9;

            for (int p = 0; p < n_phi_quad; p++) {
                const double phi = (p + 0.5) * dphi;
                const double cphi = cos(phi);
                const double sphi = sin(phi);
                const double cmphi = cos(m * phi);
                const double smphi = sin(m * phi);

                double T_local[9];
                surface_T_coxmunk_trig(mo, mi, cphi, sphi,
                                       ws, sigma_type, n_water,
                                       q_convention, T_local);

                for (int kl = 0; kl < 9; kl++) {
                    double weight = cos_entry[kl] ? cmphi : smphi;
                    T_pair[kl] += T_local[kl] * weight;
                }
            }
            for (int kl = 0; kl < 9; kl++) {
                T_pair[kl] *= dphi * norm;
            }
        }
    }
    ocrt_twa_fold_ucol_m(T_m, n_o, n_i, m);
    return 0;
}

/* ----------------------------------------------------------------------------
 * Direct sunglint single-bounce TOA reflectance contribution.
 * 
 * Two regimes:
 *   ws = 0:  TRUE flat Fresnel — specular delta function at μ_v=μ_0, Δφ=π.
 *            Outside exact specular: zero. Otherwise diverges. Implementation
 *            checks proximity (|μ_v-μ_0| < tol AND |cos(Δφ)+1| < tol) and
 *            applies F_Fresnel(μ_0) · μ_0 · trans / (small geometric factor).
 *            Strictly speaking the integral over the delta gives:
 *              ρ_glint = R_F(μ_0) · trans
 *            where R_F is the Fresnel intensity reflectance at angle of
 *            incidence θ_0. Magnitude depends on detector field-of-view.
 *            AF1982 reference appears to use an effective specular kernel
 *            (large but finite at exact specular) — match by extrapolating
 *            limit ws→0 of Cox-Munk.
 *
 *   ws > 0:  Cox-Munk slope-weighted Fresnel (use existing R_coxmunk_trig).
 * ---------------------------------------------------------------------------- */
void surface_direct_sunglint_rho(double mu_v, double phi_v,
                                  double mu_0, double phi_0,
                                  double ws, int sigma_type, double n_water,
                                  double tau_total, int q_convention,
                                  double *rho_glint) {
    double dphi = phi_v - phi_0;
    double cos_dphi = cos(dphi), sin_dphi = sin(dphi);

    double T_dn = exp(-tau_total / mu_0);
    double T_up = exp(-tau_total / mu_v);
    double trans = T_dn * T_up;

    if (ws < 1e-6) {
        /* Refactoring Tasks Task 5 path C fix (2026-05-10):
         * Peak-preserving Gaussian smoothing for wind=0 specular delta.
         *
         * 이전 식: tol_mu/tol_phi step 함수 — grid 가 specular 에서 미세히 mismatch 시
         * energy 통째로 유실 (vza=30 → 30.1 사이 50× drop 확인됨).
         *
         * Path C: 2D Gaussian smoothing for *finite instrument FOV equivalent*
         *   ρ_glint(μ_v, φ_v) = ρ_total × exp(-0.5·((Δμ/σ_μ)² + (Δφ/σ_φ)²))
         *   ρ_total = 4π · F_F · trans / (μ_0 · μ_v · sin(2θ_0))
         *   Peak (exact specular): G=1, ρ=ρ_total (이전 식 보존)
         *   Off-specular: smooth Gaussian decay, *step discontinuity 제거*
         *
         * σ_FOV: instrument FOV equivalent (default 0.5° = 0.0087 rad).
         * env V3_GLINT_FOV_DEG 으로 override (e.g. GOCI-III ≈ 0.18°).
         *
         * Physical 의미: *Real instrument 의 finite FOV* 가 *delta function* 을
         * *Gaussian-shaped peak 으로 sample*. *peak 단위 보존* 로 이전 식 결과 의
         * *exact specular case* 그대로 reproduce, *near-specular 의 step
         * discontinuity* 사라짐.
         *
         * Energy 분석: *peak-preserving form* 은 *normalized Gaussian* 보다
         * *integrated total energy 가 다름* (Gaussian spread 만큼 추가 energy).
         * 그러나 *real instrument 의 *single direction sampling* 가 *peak value*
         * 측정 — *user expectation 정합* (single grid 의 ρ_glint = peak ρ at
         * that direction). */
        double sigma_fov_deg = 0.5;     /* default FOV: 0.5° */
        const char *fov_env = getenv("V3_GLINT_FOV_DEG");
        if (fov_env) {
            double tmp = atof(fov_env);
            if (tmp > 0.0) sigma_fov_deg = tmp;
        }
        const double sigma_fov_rad = sigma_fov_deg * M_PI / 180.0;
        /* σ in μ-space: differential of cos(θ) at specular = sin(θ_0)·dθ */
        const double sin_theta_0 = sqrt(fmax(0.0, 1.0 - mu_0*mu_0));
        const double sigma_mu = sigma_fov_rad * fmax(sin_theta_0, 1e-3);
        const double sigma_phi = sigma_fov_rad;

        /* dphi wrap to [-π, π] */
        double dphi_eff = dphi;
        while (dphi_eff >  M_PI) dphi_eff -= 2.0 * M_PI;
        while (dphi_eff < -M_PI) dphi_eff += 2.0 * M_PI;
        const double dmu = mu_v - mu_0;

        /* 2D peak-preserving Gaussian (peak = 1 at specular) */
        const double exponent = -0.5 * ((dmu*dmu) / (sigma_mu*sigma_mu) +
                                         (dphi_eff*dphi_eff) / (sigma_phi*sigma_phi));
        if (exponent < -50.0) {
            rho_glint[0] = 0.0;
            rho_glint[1] = 0.0;
            rho_glint[2] = 0.0;
            return;
        }
        const double gauss_weight = exp(exponent);

        double MF[9];
        surface_flat_fresnel_matrix_raw(mu_0, n_water, q_convention, MF);

        const double sin2_theta_0 = 2.0 * mu_0 * sin_theta_0;
        const double denom = fmax(sin2_theta_0 * mu_0 * mu_v, 1e-12);
        const double rho_total = 4.0 * M_PI * trans / denom;
        const double rho_local = rho_total * gauss_weight;

        rho_glint[0] = MF[0*3+0] * rho_local;
        rho_glint[1] = MF[1*3+0] * rho_local;
        rho_glint[2] = MF[2*3+0] * rho_local;
        return;
    }

    /* ws > 0: Cox-Munk slope-weighted */
    double R[9];
    surface_R_coxmunk_trig(mu_v, mu_0, cos_dphi, sin_dphi,
                            ws, sigma_type, n_water, q_convention, R);

    /* TOA reflectance from solar specular bounce:
     *   L_v = R(μ_v,μ_0,Δφ) · μ_0 · F_solar · trans          [radiance]
     *   ρ_TOA = π · L_v / (μ_0 · F_solar)
     *         = π · R · trans
     * (μ_0 cancels.) Note: R here has units of sr^-1 (BRDF). */
    rho_glint[0] = M_PI * R[0*3+0] * trans;
    rho_glint[1] = M_PI * R[1*3+0] * trans;
    rho_glint[2] = M_PI * R[2*3+0] * trans;
}
