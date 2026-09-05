/* rt_iop_organic.h
 *
 * OCRT chlorophyll-linked organic-particle optical-property adapter.
 *
 * Current production policy (2026-08-22): EAP species scattering is disabled.
 * Chl contributes the common approved phytoplankton absorption spectrum, while
 * Stramski organic detritus supplies the Chl-linked particle phase/scattering.
 * EAP .mie files are archival/generator products and are not runtime inputs.
 *
 * DOC-REF: OCRT_EAP_RUNTIME_DISABLED_AND_WINDOWS_RUN_POLICY_2026-08-22.
 */
#ifndef RT_IOP_ORGANIC_H
#define RT_IOP_ORGANIC_H

#include "rt_water_iop.h"
#include "rt_spectral_contract.h"

#ifdef __cplusplus
extern "C" {
#endif

/* 앞의 세 항목(pico/nano/micro)은 번호와 파일명을 그대로 둔다.  기존 배치
 * 격자와 회귀 결과가 이 번호에 묶여 있기 때문이다.  2026-07-26 에 추가한
 * EAP 17종은 3번부터 이어 붙였다.  새 종은 층상구형 산란 생성기로 만든
 * 위상 파일을 쓴다(크기분포 표본 간격 0.05 um). */
typedef enum {
    ORGANIC_PHYTO_PICO  = 0,
    ORGANIC_PHYTO_NANO  = 1,
    ORGANIC_PHYTO_MICRO = 2,
    /* --- 신규 EAP 17종 (2026-07-26) --- */
    ORGANIC_PHYTO_EAP_DIATOMS_PENNATE = 3,
    ORGANIC_PHYTO_EAP_CHLOROPHYTES,
    ORGANIC_PHYTO_EAP_DIATOMS_CENTRIC,
    ORGANIC_PHYTO_EAP_CRYPTOPHYTES,
    ORGANIC_PHYTO_EAP_CYANO_BLUE,
    ORGANIC_PHYTO_EAP_CYANO_RED,
    ORGANIC_PHYTO_EAP_DINOFLAGELLATES,
    ORGANIC_PHYTO_EAP_EUSTIGMATOPHYTES,
    ORGANIC_PHYTO_EAP_HAPTO_PAVLOVACEAE,
    ORGANIC_PHYTO_EAP_PELAGOPHYTES,
    ORGANIC_PHYTO_EAP_PRASINOPHYTES,
    ORGANIC_PHYTO_EAP_PROCHLOROCOCCUS,
    ORGANIC_PHYTO_EAP_HAPTO_PRYMNESIACEAE,
    ORGANIC_PHYTO_EAP_RAPHIDOPHYTES,
    ORGANIC_PHYTO_EAP_RHODOPHYTES,
    ORGANIC_PHYTO_EAP_SYNECHOCOCCUS,
    ORGANIC_PHYTO_EAP_MICROCYSTIS,
    ORGANIC_PHYTO_GROUP_COUNT
} organic_phyto_group_t;

#define ORGANIC_DETRITUS_SLOPE_DEFAULT 0.0109
#define ORGANIC_DETRITUS_SLOPE_MIN     0.0024
#define ORGANIC_DETRITUS_SLOPE_MAX     0.0170
#define ORGANIC_WAVELENGTH_MIN_NM RT_SPECTRAL_MIN_NM
#define ORGANIC_WAVELENGTH_MAX_NM RT_SPECTRAL_MAX_NM
/* The distributed PLOPS-derived table is explicit through 1100 nm.  The
 * project taper reaches exact zero at 800 nm and all longer wavelengths are
 * represented by actual zero rows in the data file, not a runtime branch. */
#define ORGANIC_PHYTO_ABS_ZERO_START_NM 800.0
#define ORGANIC_PHYTO_ABS_TABLE_MAX_NM RT_SPECTRAL_MAX_NM

/* Load the common phytoplankton absorption table and Stramski detritus.
 * No EAP .mie file is opened. */
int  rt_iop_organic_init(const char *data_dir);
/* ABI-compatible initializer. selected_group/load_all_phyto are ignored while
 * EAP runtime scattering is disabled. */
int  rt_iop_organic_init_selected(const char *data_dir,
                                  organic_phyto_group_t selected_group,
                                  int load_all_phyto);
void rt_iop_organic_free(void);
int  rt_iop_organic_ready(void);

int rt_iop_organic_group_parse(const char *name, organic_phyto_group_t *out);
organic_phyto_group_t rt_iop_organic_group_from_name(const char *name);
const char *rt_iop_organic_group_name(organic_phyto_group_t group);

/* 대형 식물플랑크톤 종인지 알려준다.  1 이면 대형이라 지금은 쓸 수 없다.
 *
 * 이유는 수치 표현 한계다.  수중 솔버는 위상함수를 L=200 차 르장드르 전개로
 * 받는데(FIXED_BULK_BETAL_LMAX), 유효직경이 큰 종은 전방 첨두가 가팔라 이
 * 차수로 담기지 않는다.  그 결과 전개를 되살린 위상이 중간 산란각에서 음수가
 * 되며, 이는 물리적으로 불가능한 값이다.  구성모델 경로는 델타 절단을 걸 수
 * 없으므로(혼합 유기 위상 미지원) 이를 피할 방법이 현재 없다.
 *
 * 판정 근거(2026-07-26 측정): 412·443·555 nm 에서 L=200 전개를 되살려 30~150도
 * 구간의 음수 개수를 세었다.  유효직경 1.2 um 이하 세 종만 음수가 0개였고,
 * 3 um 이상은 모두 42~121개의 음수가 나왔다. */
int rt_iop_organic_group_is_large(organic_phyto_group_t group);
const char *rt_iop_organic_phyto_phase_path(organic_phyto_group_t group);
const char *rt_iop_organic_detritus_phase_path(void);
/* The active production detritus phase retains its validated native range.
 * The 330--1100 recipe-regenerated candidate is distributed separately under
 * inputs/water_iop/candidates after failing the L=200/SOS acceptance gate. */
int rt_iop_organic_detritus_wavelength_supported(double lambda_nm);
double rt_iop_organic_detritus_wavelength_min_nm(void);
double rt_iop_organic_detritus_wavelength_max_nm(void);
int rt_iop_organic_wavelength_supported(double lambda_nm);

/* Huot et al. (2008) total particulate backscattering [m^-1]. */
double rt_iop_huot_bbp(double lambda_nm, double chl_mg_m3);

/* Phytoplankton fraction of total particulate bbp. */
double rt_iop_fph_fraction(double chl_mg_m3);

/* Disabled EAP scattering API. Returns -8 (EAP_RUNTIME_DISABLED). */
int rt_iop_eap_phyto_eval(double lambda_nm, double chl_mg_m3,
                          organic_phyto_group_t group, rt_iop_t *out);

/* Disabled EAP scattering API. Returns -8 (EAP_RUNTIME_DISABLED). */
int rt_iop_eap_phyto_eval_with_phase_ratios(
    double lambda_nm, double chl_mg_m3, organic_phyto_group_t group,
    double bb_ratio_lambda, double bb_ratio_550, rt_iop_t *out);
/* Absorption-only production path used while species-specific EAP scattering is
 * disabled.  This deliberately does not evaluate phase moments or bb/b. */
int rt_iop_eap_phyto_absorption_eval(
    double lambda_nm, double chl_mg_m3, organic_phyto_group_t group,
    double *a_out_m_inv);

/* Organic detritus IOP at one wavelength.
 * Scattering follows the Huot residual allocation.  Absorption uses the
 * caller-supplied a_d(440) and the Bricaud-Stramski exponential slope.
 * a_d440=0 is valid and means that only the Chl-covarying detrital scattering
 * is included. */
int rt_iop_detritus_eval(double lambda_nm, double chl_mg_m3,
                         double a_d440_m_inv, double slope_nm_inv,
                         rt_iop_t *out);

int rt_iop_detritus_eval_with_phase_ratios(
    double lambda_nm, double chl_mg_m3,
    double a_d440_m_inv, double slope_nm_inv,
    double bb_ratio_lambda, double bb_ratio_550, rt_iop_t *out);

double rt_iop_detritus_abs_shape(double lambda_nm, double slope_nm_inv);

#ifdef __cplusplus
}
#endif
#endif /* RT_IOP_ORGANIC_H */
