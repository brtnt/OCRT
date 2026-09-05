/* rt_water_iop.h
 *
 * Pure water inherent optical properties (IOP) module for OCRT v1.02
 * underwater radiative transfer (Phase A).
 *
 * Sources:
 *   - Absorption a_w(λ):  Pope & Fry 1997 (380-700 nm) + Kou et al. 1993 (others).
 *                          Distributed via NASA OceanColor as water_coef_fine_z09.txt
 *                          (source: 200-2450 nm at 0.1 nm; distributed LUT: 200-2449 nm at 1 nm band averages).
 *   - Scattering b_w(λ):  Zhang, Hu & He 2009 (J. Geophys. Res. Oceans 114, C12025;
 *                          T=20°C, S=38.4 g/kg). bbw = 0.5*bw per Rayleigh-like
 *                          symmetric phase function.
 *   - Refractive index n_w(λ): Quan & Fry 1995 (T=20°C, S=38.4 g/kg) — real part only.
 *                              Zhang 2009 does NOT provide imaginary k_w; k_w can be
 *                              derived from a_w as k_w = a_w·λ / (4π) but is not
 *                              currently used (in-water RT uses a_w directly).
 *   - Depolarization δ_w(λ): Zhang 2009 ~ 0.039 (weakly wavelength-dependent).
 *
 * Mueller matrix model (pure-water scattering, vector mode):
 *   Z(θ) = (1/2) · [
 *      P11(θ),     P12(θ),     0       ;
 *      P12(θ),     P22(θ),     0       ;
 *      0,          0,          P33(θ)
 *   ]
 *   with Rayleigh-like depolarized scattering:
 *     P11(θ) = K · [(1 - δ_w) · (1 + cos²θ) + 2δ_w]
 *     P12(θ) = -K · (1 - δ_w) · sin²θ
 *     P22(θ) = K · [(1 - δ_w) · (1 + cos²θ) + 2δ_w]   (same as P11 here)
 *     P33(θ) = K · 2(1 - δ_w) · cosθ
 *     K      = 3 / (8π · (1 + 2δ_w / (1 - δ_w)))    [normalization to ∫P11 dΩ = 4π]
 *
 * Linear interpolation in λ over the loaded LUT. Nearest only in extrapolation
 * (outside LUT range; logged as warning).
 *
 * Phase A (current): pure water IOP only. Future Phase B will add seawater
 * constituents (Chl, TSM, CDOM, etc.) via the same Mueller-matrix infrastructure.
 *
 * Author: OCRT team
 * Date: 2026-05-22 KST
 */
#ifndef OCRT_RT_WATER_IOP_H
#define OCRT_RT_WATER_IOP_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================
 * Water IOP lookup table (loaded from NASA water_coef_fine_z09.txt)
 * ============================================================ */
typedef struct {
    int     n;              /* number of wavelength samples */
    double* lambda_nm;      /* [n] wavelengths in nm */
    double* aw_m_inv;       /* [n] absorption coefficient, m⁻¹ */
    double* bw_m_inv;       /* [n] total scattering coefficient, m⁻¹ */
    /* bbw = 0.5*bw computed on demand (per file note) */
} rt_water_iop_lut_t;

/* ============================================================
 * Temperature correction LUT (ψ_T(λ), units m⁻¹·°C⁻¹)
 * ============================================================
 * a_w(T, λ) = a_w(20, λ) + ψ_T(λ) · (T - 20)
 *
 * Reference T_ref = 20°C (matches NASA water_coef_fine_z09.txt).
 * Possible source datasets:
 *   - Röttgers et al. 2014 (Opt. Express 22, 25093): 400-2700 nm hyperspectral
 *   - Pegau, Gray & Zaneveld 1997 (Appl. Opt. 36, 6035): 15 wavelengths VIS+NIR
 *   - Sullivan et al. 2006 (Appl. Opt. 45, 5294): 400-750 nm hyperspectral
 *
 * Default behavior when no ψ_T LUT is provided: ψ_T = 0 (no T correction,
 * a_w identical to NASA file at all T).
 *
 * Expected text file format:
 *   /begin_header
 *   ! comments...
 *   /fields=wavelength,psi_T
 *   /units=nm,m^-1_per_degC
 *   /end_header
 *   <lambda_nm>  <psi_T_value>
 *   ...
 */
typedef struct {
    int     n;              /* number of wavelength samples (0 = unset) */
    double* lambda_nm;      /* [n] wavelengths in nm */
    double* psi_T;          /* [n] dψ_T/dT  (m⁻¹/°C) */
    double  T_ref;          /* reference temperature, °C (default 20) */
} rt_water_iop_psi_T_lut_t;

/* Allocate and parse ψ_T LUT file (same format as NASA water file).
 * Returns 0 on success, nonzero on error. */
int  rt_water_iop_psi_T_load(const char* filepath, rt_water_iop_psi_T_lut_t* lut);

void rt_water_iop_psi_T_free(rt_water_iop_psi_T_lut_t* lut);

/* Allocate and parse NASA water_coef_fine_z09.txt format file.
 * Returns 0 on success, nonzero on error.
 *
 * Expected format:
 *   /begin_header
 *   ! comments...
 *   /fields=wavelength,aw,bw
 *   /units=nm,m^-1,m^-1
 *   /end_header
 *   <lambda_nm> <aw_m_inv> <bw_m_inv>
 *   ...
 */
int  rt_water_iop_lut_load(const char* filepath, rt_water_iop_lut_t* lut);

void rt_water_iop_lut_free(rt_water_iop_lut_t* lut);

/* ============================================================
 * Wavelength-dependent pure water IOP (linear interp over LUT)
 * ============================================================ */

/* Returns interpolated a_w at lambda_nm (m⁻¹).
 * Outside LUT range: returns nearest-edge value and sets *extrapolated=1 (warning). */
double rt_water_iop_aw(const rt_water_iop_lut_t* lut, double lambda_nm,
                       int* extrapolated /* OUT, may be NULL */);

/* Temperature-corrected absorption at user temperature T (°C).
 * Formula: a_w(T, λ) = a_w(20, λ) + ψ_T(λ) · (T - T_ref)
 * psi_T_lut may be NULL → ψ_T = 0 (returns a_w(20) unchanged).
 * Linear interp in λ for both LUTs; nearest in extrapolation. */
double rt_water_iop_aw_T(const rt_water_iop_lut_t*       aw_lut,
                         const rt_water_iop_psi_T_lut_t* psi_T_lut /* may be NULL */,
                         double lambda_nm,
                         double T_celsius,
                         int*   extrapolated /* OUT, may be NULL */);

double rt_water_iop_bw(const rt_water_iop_lut_t* lut, double lambda_nm,
                       int* extrapolated /* OUT, may be NULL */);

/* bbw = 0.5 * bw */
static inline double rt_water_iop_bbw(const rt_water_iop_lut_t* lut, double lambda_nm,
                                       int* extrapolated) {
    return 0.5 * rt_water_iop_bw(lut, lambda_nm, extrapolated);
}

/* Real refractive index of pure seawater (Quan & Fry 1995, T=20°C, S=38.4 g/kg).
 * Formula (lambda in nm):
 *   n_w = n0 + (n1 + n2*T + n3*T²) * S + n4*T² + (n5 + n6*S + n7*T)/λ
 *         + n8/λ² + n9/λ³
 * For T=20, S=38.4 default coefficients reduce to a simple λ-only polynomial. */
double rt_water_iop_n_real(double lambda_nm);

/* Depolarization ratio δ_w(λ).  Zhang 2009 reports ~0.039 with weak λ-dependence;
 * we use a 2-term fit from their Table 1 (or constant 0.039 as fallback). */
double rt_water_iop_depol(double lambda_nm);

/* ============================================================
 * Pure-water Mueller scattering matrix Z(θ;λ)
 * ============================================================
 * Fills Z[9] (row-major 3x3): [I,Q,U] convention, V channel ignored.
 *
 *   Z[0] = P11    Z[1] = P12   Z[2] = 0
 *   Z[3] = P12    Z[4] = P22   Z[5] = 0
 *   Z[6] = 0      Z[7] = 0     Z[8] = P33
 *
 * Normalized such that (1/4π) ∫ P11 dΩ = 1.
 */
void rt_water_iop_mueller(double theta_deg, double lambda_nm,
                          double Z_out[9]);

/* ============================================================
 * Component-based IOP interface (B.5+, 2026-05-23)
 * ============================================================
 *
 * 각 in-water component (pure water, CDOM, phytoplankton, NAP …) 가
 * 동일한 (a, b, bb) 자료구조를 채우는 통일된 인터페이스. aggregator 가
 * 단순 합산하여 RT 솔버에 전달.
 *
 * 약속:
 *   - 모든 단위 m⁻¹.
 *   - "scattering 이 없는" component (예: CDOM 흡광) 는 b=bb=0 으로 둔다.
 *   - "absorption 이 없는" component (가상의 순수 산란체) 는 a=0.
 *   - 일관성: bb ≤ 0.5·b 가 권장되나 (back-scattering 분율 ≤ 1) RT 솔버는
 *     bb 를 m=0 차수에서만 별도 사용하지 않으므로 강제 검사하지 않는다.
 *
 * 새 component 추가 절차:
 *   1) rt_iop_components_t 에 필드 추가
 *   2) rt_iop_<component>_eval() 함수 작성 (a/b/bb 계산)
 *   3) rt_iop_components_eval_all() 안에 분기 추가
 * Mueller matrix 단계 통합은 B.6 이후 단계로 분리.
 */
typedef struct {
    double a;       /* absorption coefficient (m⁻¹) */
    double b;       /* total scattering coefficient (m⁻¹) */
    double bb;      /* back-scattering coefficient (m⁻¹) */
} rt_iop_t;

/* Public water-input branch selected before entering the common OCRT
 * underwater RT solver.  UNSET is intentionally a real state so the CLI can
 * distinguish an omitted branch selector from an explicitly selected model.
 *
 * OCRT/CCRR are constituent-to-optics adapters.  When their three required
 * constituent values (Chl, TSM, aDOM440) are all exactly zero, main.c lowers
 * the request to the native pure-water execution path while retaining this
 * selection as provenance.  IOP injects total bulk a/b/bb directly. */
typedef enum {
    RT_WATER_INPUT_UNSET = 0,
    RT_WATER_INPUT_OCRT  = 1,
    RT_WATER_INPUT_CCRR  = 2,
    RT_WATER_INPUT_IOP   = 3
} rt_water_input_mode_t;

/* Constituent-to-optics adapter selected before the common OCRT RT solver. */
typedef enum {
    RT_WATER_CONSTITUENT_CCRR = 0,
    RT_WATER_CONSTITUENT_OCRT = 1
} rt_water_constituent_model_t;

/* 모든 component 의 a/b/bb 집합. RT 솔버는 rt_iop_total() 결과만 사용. */
typedef struct {
    rt_iop_t pure_water;    /* a_w(λ,T) + b_w(λ); bb_w = 0.5·b_w (Rayleigh-like) */
    rt_iop_t cdom;          /* 2-param exponential: a only, b=bb=0 */
    rt_iop_t pigment;       /* CCRR pigment contribution */
    rt_iop_t eap_phyto;     /* OCRT EAP phytoplankton */
    rt_iop_t detritus;      /* OCRT Chl-covarying organic detritus */
    rt_iop_t mineral;       /* TSM/mineral: OCRT Ahn or CCRR empirical contribution */
} rt_iop_components_t;

/* CDOM 모델 파라미터 (Phase B.5):
 *   a_CDOM(λ) = a440 · exp(-S · (λ - λ_ref))
 * a440 ≤ 0 → CDOM 비활성 (a=0 반환). */
typedef struct {
    double a440;            /* CDOM absorption at λ_ref (m⁻¹) */
    double S_nm_inv;        /* spectral slope (nm⁻¹), e.g. 0.014 (Bricaud 1981) */
    double lambda_ref_nm;   /* reference wavelength (nm), default 440 */
} rt_iop_cdom_params_t;

/* ---- Per-component evaluators ------------------------------- */

/* Pure seawater (a_w via Pope & Fry + Röttgers ψ_T, b_w via Zhang 2009).
 * NULL psi_T_lut → no T correction. extrapolated 는 NULL 허용. */
int rt_iop_pure_water_eval(const rt_water_iop_lut_t*       aw_lut,
                           const rt_water_iop_psi_T_lut_t* psi_T_lut,
                           double lambda_nm,
                           double T_celsius,
                           rt_iop_t*                       out,
                           int*                            extrapolated /* may be NULL */);

/* CDOM 2-parameter exponential absorption. b=bb=0. */
int rt_iop_cdom_eval(const rt_iop_cdom_params_t* params,
                     double lambda_nm,
                     rt_iop_t*                   out);

/* CCRR pigment/mineral IOP evaluators. These port only the pre-RT
 * constituent-to-IOP conversion; the RT solver remains OCRT. Units:
 *   chl_mg_m3 : mg m^-3
 *   min_g_m3  : g m^-3
 * The canonical Chl table is inputs/water_iop/aph_ccrr_morel1988_mm01.txt.
 * It represents the Morel (1988) normalized spectrum through the classic
 * Case-1 closure a_p=0.06*A_chl(lambda)*Chl^0.65. The historical filename
 * aph_bricaud_1998.txt is accepted only as a compatibility fallback. */
int rt_iop_ccrr_pigment_eval(double lambda_nm, double chl_mg_m3, rt_iop_t* out);
int rt_iop_ccrr_mineral_eval(double lambda_nm, double min_g_m3, rt_iop_t* out);

/* HG scalar phase helper for CCRR particle placeholder phase in OCRT. Returns
 * Legendre coefficient B_l=(2l+1)g^l where g is chosen to match the requested
 * backscatter fraction for a Henyey-Greenstein phase. */
double rt_iop_hg_g_for_backscatter_fraction(double target_bb_fraction);

/* CCRR/IOCCG21 scalar particle phase moments.
 * Preferred path: call rt_iop_ccrr_phase_moments_load() with a CSV exported
 * by the CCRR pre-solver.  Required columns are:
 *   phase_name,l,chi_l
 * Optional aliases: chi, moment, betal_l, betal, beta_l.  If betal is given,
 * it is used directly; if chi_l is given, Betal[l]=(2l+1)chi_l.
 * Phase names containing pig/ff map to kind=0; names containing min/petzold
 * map to kind=1.  Pure-water moments are ignored: OCRT keeps its vector
 * Rayleigh-like seawater phase in CCRR mode.
 *
 * Fallback path: if no CSV is loaded, use embedded proxy arrays.
 * Returns Betal coefficients for P11(mu)=Σ_l Betal[l]P_l(mu). */
int rt_iop_ccrr_phase_moments_load(const char* csv_path);
int rt_iop_ccrr_particle_betal(int kind, int n_out, double* betal_out);

/* ---- Aggregator -------------------------------------------- */

/* 모든 component IOP 의 단순 합산.
 *   a_total  = Σ a_i
 *   b_total  = Σ b_i
 *   bb_total = Σ bb_i
 *
 * 미래에 phase function 통합이 필요해질 때까지는 이 합산만으로 충분
 * (CDOM 은 b=bb=0 이므로 phase function 가중치 0). */
rt_iop_t rt_iop_total(const rt_iop_components_t* comp);

#ifdef __cplusplus
}
#endif

#endif /* OCRT_RT_WATER_IOP_H */
