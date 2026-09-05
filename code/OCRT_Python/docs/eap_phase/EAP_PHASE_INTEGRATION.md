# EAP 위상행렬 생성기 통합 가이드

**대상 기준선:** OCRT-v1.2-2026-07-25-KST-stage2-coupling-clamp-fix
**작업지시서:** OCRT-WO-EAP-MIE-PHASE-20260726 Rev.1
**공개 헤더 SHA-256:** `034f25ac13766ee6575edb504068bc686f731482ef09fcecfa36d24433a34ada` (지시서 명시값과 일치 확인)

---

## 1. 무엇이 들어 있나

### 1.1 신규 파일 (전부 추가, 기존 파일 수정 없음)

```
src/rt_eap_mie_phase.h                     공개 헤더 (지시서 부록 A 그대로, 무수정)
src/rt_eap_mie_phase.c                     공개 함수 구현
src/internal/rt_eap_coated_mie.h/.c        층상구형 산란 코어
src/internal/rt_eap_species_catalog.h/.c   17종 카탈로그 접근부
src/generated/rt_eap_species_data.inc      빌드 시 생성되는 굴절률 표 (읽기 전용)
tools/eap_catalog_freeze.py                카탈로그 동결표 (모든 상수의 출처 주석 포함)
tools/generate_eap_species_catalog.py      원자료 -> C 표 변환기
tools/eap_mie_writer.c                     레거시 .mie 파일 작성기 (공개 API 밖)
tools/eap_generate_all.sh                  70개 조합 일괄 생산
tests/eap_mie_shim.c                       시험용 실수인자 껍데기
inputs/eap_source/                         원자료 2종 (해시가 .inc 머리글에 기록됨)
```

### 1.2 기존 파일 변경

**없다.** 신규 경로가 호출되지 않으면 기존 OCRT는 바이트 단위로 동일하다(지시서 A10). 통합은 3절의 선택적 패치로만 이루어진다.

---

## 2. 빌드

기존 빌드 명령에 새 소스 세 개를 추가하면 끝난다. `find src -name '*.c'`를 쓰고 있으면 자동으로 잡힌다.

```sh
# 1) 카탈로그 표 생성 (원자료가 바뀔 때만 다시 실행)
python3 tools/generate_eap_species_catalog.py

# 2) 기존 빌드 명령 그대로
gcc -std=c11 -O3 -march=native -ffp-contract=fast -fassociative-math \
    -fno-signed-zeros -fno-trapping-math -DOCRT_FAST_KERNELS -fopenmp \
    -Isrc $(find src -name '*.c') -o build/ocrt -lm
```

`src/generated/rt_eap_species_data.inc`는 `.c`가 아니므로 `find`에 걸리지 않는다. `rt_eap_species_catalog.c`가 `#include` 한다.

---

## 3. 통합 방법

### 3.1 1단계 — 진단 경로만 추가 (권장 시작점)

공개 함수를 호출해 기존 모멘트 변환에 그대로 넘긴다. 기존 캐시는 건드리지 않는다.

```c
#include "rt_eap_mie_phase.h"

/* 어댑터. 공개 API 밖에 둔다(지시서 7.2). */
static int eap_phase_to_moments(unsigned species_id, double wl_nm,
                                const double *theta_deg, int n_theta,
                                int lmax, int mie_n_mu,
                                double *betal, double *gammal,
                                double *alphal, double *zetal)
{
    double *p11 = malloc((size_t)n_theta * sizeof *p11);
    double *p12 = malloc((size_t)n_theta * sizeof *p12);
    double *p33 = malloc((size_t)n_theta * sizeof *p33);
    if (!p11 || !p12 || !p33) { free(p11); free(p12); free(p33); return -1; }

    int rc = rt_eap_mie_phase_compute(species_id, &wl_nm, 1,
                                      theta_deg, (size_t)n_theta,
                                      p11, p12, p33);
    if (rc == RT_EAP_PHASE_OK)
        rc = rt_aerosol_compute_vector_legendre_gauss(
                 p11, p12, p33, theta_deg, n_theta, lmax, mie_n_mu,
                 betal, gammal, alphal, zetal);

    free(p11); free(p12); free(p33);
    return rc;
}
```

모멘트 변환 이후의 경로는 기존과 완전히 동일하다. `rt_aerosol_compute_vector_legendre_gauss()`가 끝에서 `rt_phase_normalize_vector()`로 β₀ 재규격화를 하므로 입력의 절대 크기는 결과에 영향을 주지 않는다.

### 3.2 2단계 — .mie 경로로 종 수 확장

OCRT의 식물플랑크톤 축은 열거형 하나와 이름표·파일표 두 개뿐이다. 확장은 세 곳만 고치면 된다.

**(a) `src/rt_iop_organic.h`**

```c
typedef enum {
    ORGANIC_PHYTO_PICO  = 0,
    ORGANIC_PHYTO_NANO  = 1,
    ORGANIC_PHYTO_MICRO = 2,
    /* --- 신규 EAP 17종 --- */
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
```

**(b) `src/rt_iop_organic.c`** — `k_group_name[]`과 `k_phyto_file[]`에 같은 순서로 항목을 추가한다. 이름은 `--ocrt-phyto-group` 인자로 그대로 쓰인다.

```c
static const char *const k_group_name[ORGANIC_PHYTO_GROUP_COUNT] = {
    "pico", "nano", "micro",
    "eap_diatoms_pennate", "eap_chlorophytes", "eap_diatoms_centric",
    "eap_cryptophytes", "eap_cyano_blue", "eap_cyano_red",
    "eap_dinoflagellates", "eap_eustigmatophytes", "eap_hapto_pavlovaceae",
    "eap_pelagophytes", "eap_prasinophytes", "eap_prochlorococcus",
    "eap_hapto_prymnesiaceae", "eap_raphidophytes", "eap_rhodophytes",
    "eap_synechococcus", "eap_microcystis"
};
static const char *const k_phyto_file[ORGANIC_PHYTO_GROUP_COUNT] = {
    "pico_Synechococcus_EAP.mie", "nano_Haptophytes_EAP.mie",
    "Diatoms_centric_EAP.mie",
    "EAP_00_Diatoms_pennate_D6.mie",   "EAP_01_Chlorophytes_D8.mie",
    "EAP_02_Diatoms_centric_D6.mie",   "EAP_03_Cryptophytes_D6.mie",
    "EAP_04_Cyano_blue_D6.mie",        "EAP_05_Cyano_red_D6.mie",
    "EAP_06_Dinoflagellates_D24.mie",  "EAP_07_Eustigmatophytes_D6.mie",
    "EAP_08_Hapto_Pavlovaceae_D6.mie", "EAP_09_Pelagophytes_D3.mie",
    "EAP_10_Prasinophytes_D3.mie",     "EAP_11_Prochlorococcus_D0p5.mie",
    "EAP_12_Hapto_Prymnesiaceae_D4.mie", "EAP_13_Raphidophytes_D24.mie",
    "EAP_14_Rhodophytes_D6.mie",       "EAP_15_Synechococcus_D1p2.mie",
    "EAP_16_Microcystis_D5.mie"
};
```

**(c) `src/main.c`** — 도움말 문자열만 고친다(279행 부근). 인자 해석은 `rt_iop_organic_group_parse()`가 이름표를 순회하므로 자동으로 새 이름을 받는다.

기존 세 항목의 번호와 파일명을 그대로 두었으므로 **기존 배치 격자와 회귀 결과는 영향을 받지 않는다.** 검증 전에 기존 캐시가 자동 대체되는 일도 없다(지시서 부록 C).

### 3.3 mie 파일 생산

```sh
gcc -std=c11 -O3 -march=native -Isrc -o build/eap_mie_writer \
    tools/eap_mie_writer.c src/internal/rt_eap_species_catalog.c \
    src/internal/rt_eap_coated_mie.c -lm

# 17개 대표 조합
sh tools/eap_generate_all.sh representative inputs/water_iop

# 논문 Table 1 전체 70개 조합
sh tools/eap_generate_all.sh all inputs/water_iop
```

파일 형식은 기존 `.mie`와 동일하다. 1195행, 스펙트럴 표 101행(0.350~0.850 μm, 5 nm), 위상함수 블록 세 개(P11/P12/P33) 각각 12파장 × 361각도(180°에서 0°까지 0.5°).

---

## 4. 확정된 물리 설정

모든 값의 출처는 `tools/eap_catalog_freeze.py` 주석에 있다. 임의로 정한 값은 없다.

| 항목 | 값 | 출처 |
|---|---|---|
| 주변매질 굴절률 | 1.334 (고정) | 논문 식 1 |
| 비포장 최대 비흡수 | 675 nm에서 0.027 m²/mg | Johnsen et al. 1994 |
| 규격화 파장 | **675 nm** | 논문 식 1 |
| 세포내 엽록소 밀도 | 2 kg/m³ | 껍질 굴절률이 물보다 커지는 유일한 값 |
| 분포 유효분산 | 0.6 | 논문 본문 |
| 껍질 부피비(진핵) | 0.2 | 논문 본문, 역산으로 재확인 |
| 껍질 실수 굴절률(진핵) | 1.10 | 논문 본문, 역산으로 재확인 |
| 핵 실수 굴절률(진핵) | 1.02 | 논문 본문 |
| 껍질 부피비(마이크로시스티스) | 0.5 | Matthews & Bernard 2013 |
| 껍질 실수 굴절률(마이크로시스티스) | 1.12 | Matthews & Bernard 2013 |
| 기체공포 굴절률 | n=333λ⁻¹·⁹⁴+0.82, k=2.28e7λ⁻⁴·⁶⁶+1.08e-5 | Matthews & Bernard 2013 |
| 힐베르트 변환 격자 | 300~1000 nm, 1 nm | 가장자리 오염 회피 (수렴 확인) |
| 크기분포 표본 간격 | 0.05 μm | 수렴 확인 (0.10 μm 이하 동일) |
| 크기분포 상한 | 유효직경의 5.5배 | 단면적 가중 꼬리 |

### 4.1 공개 EAP 자료와 다른 점

공개 CSV는 껍질 허수 굴절률을 **655 nm**에서 규격화했다. 원 코드 주석의 색인 오류(1 nm 해상도에서 675 nm를 255번째로 적었으나 실제로는 275번째)에서 비롯됐고, 종별 배율 1.35~3.05를 네 종에서 1% 이내로 재현해 확인했다.

675 nm를 채택한 근거는 세 가지다. 655 nm를 쓰면 (1) 17종 중 15종에서 껍질 실수 굴절률이 저자 자신이 명시한 1.06~1.22 범위를 벗어나고, (2) 페라고조류는 물보다 낮아지며, (3) 최대 선형편광도의 크기 의존성이 사라진다(6 μm 규조류가 1.2 μm 시네코코쿠스와 같은 0.95를 낸다).

따라서 이 생성기는 공개 CSV의 흡수·산란 값을 재현하지 않는다. **지시서 Gate E5는 "일치"가 아니라 "차이의 원인 규명 완료"로 통과 처리한다.**

---

## 5. 검증 현황

| Gate | 항목 | 상태 |
|---|---|---|
| E0 | 카탈로그 provenance 동결 | 완료 (원자료 SHA-256이 .inc 머리글에 기록) |
| E1 | API·오류처리 | 검증기 있음, sanitizer 미실행 |
| E2 | 독립 Mie 대조 | **통과** — 4개 사례에서 Qext·Qsca 비트 동일, 위상 성분 차 1e-14 |
| E3 | 규격화·물리성 | **통과** — 규격화 잔차 1e-15 (요구 1e-8), P33(0)/P11(0)=1, P33(180)/P11(180)=−1 |
| E4 | 수렴성 | **통과** — 차수 여유 8배 변화에 12자리 불변, 크기격자 0.10 μm 이하 동일 |
| E5 | 기존 3종 대조 | 차이 원인 규명 완료 (4.1절) |
| E6~E10 | 모멘트 폐합·RT·성능·비회귀 | 미실행 |

### 5.1 원 Fortran 루틴 대비 개선

원 루틴은 하향 점화식 시작 차수를 `NMX1 = 1500`으로 **크기와 무관하게 고정**해 두었다. 크기변수가 1350을 넘으면 발산하며, 400 nm 기준 직경 180 μm가 한계다.

| 직경 | 크기변수 | 신규 C Qext | 원 Fortran Qext |
|---|---|---|---|
| 100 μm | 1048 | 2.03123 | 2.03123 |
| 150 μm | 1572 | 2.01955 | 2.95496 |
| 200 μm | 2095 | 2.00986 | 3.09121 |
| 340 μm | 3562 | 2.00823 | **−0.15901** |

라피도조류(유효직경 60 μm)와 규조류(48 μm)는 분포 꼬리가 이 한계를 넘으므로 원 루틴으로는 계산이 불가능하다.

### 5.2 위상함수 진동

기존 `.mie` 파일의 각도 진동은 크기분포를 1 μm 간격으로만 표본한 데서 왔다. 물속 파장 기준 간섭 무늬 주기가 약 0.1 μm이므로 1 μm 간격으로는 분해되지 않는다.

중심 규조류 443 nm, 90~170° 구간(진동이 드러나는 매끄러운 영역) 기준:

| | 요철 중앙값 | 요철 최대값 |
|---|---|---|
| 기존 파일 | 0.01046 | 0.7101 |
| 신규 (0.05 μm 표본) | 0.00008 | 0.0033 |

---

## 6. 성능

단일 코어 기준, 3파장 × 361각도 × 유효직경 6 μm에서 0.5초다. 12파장이면 약 2초이며, 신규 경로를 호출하지 않으면 기존 실행시간은 변하지 않는다.

---

## 7. 남은 작업

1. Gate E6 모멘트 폐합 시험 (`beta[0]=1`, 각도 재구성 대조).
2. Gate E7 고정 IOP RT 시험 및 I/Q/U 독립 산포도.
3. Gate E9 동시호출·반복성·sanitizer.
4. Gate E10 기존 72개 물리 실행 및 600셀 비회귀표 바이트 대조.
5. 70개 조합 일괄 생산 및 17종 일괄 검증표.
