> **현재 상태 (OCRT v1.18, 2026-07-18)**
> 아래 문서는 첨부 유기모듈의 원래 통합 가이드를 보존한다. 실제 OCRT 트리에는
> EAP phytoplankton, organic detritus, Ahn TSM의 P11/P12/P33 산란계수 가중 혼합과
> phase cache가 이미 통합되어 있다. 현재 public CLI는 `--water-model ocrt`와
> `--ocrt-chl`, `--ocrt-tsm`, `--ocrt-adom440`을 사용한다. 정확한 실행 계약은
> `WATER_INPUT_INTERFACE_v1.18_2026-07-18.md`가 우선한다.

# OCRT 유기 성분 통합 가이드 (EAP 식물플랑크톤 + 유기 detritus)

이 패키지는 OCRT v1.11에 수중 유기 입자(식물플랑크톤 + 유기 detritus) IOP를
추가한다. 무기물(TSM, `rt_iop_ahn_mineral`)은 별도 모듈이며 이미 전달됐다.
본 모듈(`rt_iop_organic`)은 TSM과 독립이다.

식물플랑크톤과 detritus는 Chl과 covary하므로 한 모듈로 묶었다(f_ph·Huot 공유).
무기물은 Chl 독립이라 별도 모듈이다. 이 분리가 Case-1/Case-2 구조에 부합한다.

물리 모델·레퍼런스는 `docs/OCRT_ORGANIC_IOP_SPEC.md` 및 `src/rt_iop_organic.h`
헤더 주석 참조. 본 문서는 통합 절차만 기술한다.

---

## 1. 패키지 구성

```
src/
  rt_iop_organic.c    유기 성분 IOP 평가 (EAP phyto + detritus)
  rt_iop_organic.h    인터페이스 + 문서 수준 주석
inputs/water_iop/
  pico_Synechococcus_EAP.mie    식물플랑크톤 pico (Brewin 소형)
  nano_Haptophytes_EAP.mie      식물플랑크톤 nano (Brewin 중형)
  Diatoms_centric_EAP.mie       식물플랑크톤 micro (Brewin 대형, 연안 우점)
  Detritus_Stramski2001.mie     유기 detritus
docs/
  OCRT_ORGANIC_IOP_SPEC.md      물리 모델 명세
```

모든 `.mie`는 350-850nm, 361각, P11/P12/P33 포함. OCRT 시작밴드 350nm와 정합.

주의: detritus `.mie`의 bulk Ext/Sca는 placeholder(=1.0)이다. detritus 절대
스케일은 Huot 2008(산란)과 Bricaud & Stramski 1990(흡광)에서 오며, `.mie`는
위상함수(P11/P12/P33)와 g만 사용한다. 컬럼헤더에 명시돼 있다.

---

## 2. 빌드

`rt_iop_organic.c`를 OCRT `src/`에 넣고 기존 빌드에 포함한다. 추가 의존성 없다
(`rt_water_iop.h`의 `rt_iop_t`만 사용).

```
gcc -std=c11 -O3 -march=native -ffp-contract=fast -fassociative-math \
    -fno-signed-zeros -fno-trapping-math -DOCRT_FAST_KERNELS -fopenmp \
    -Isrc $(find src -name '*.c') -o build/<bin> -lm
```

`find src -name '*.c'`가 `rt_iop_organic.c`를 자동 포함하므로 별도 조치 불필요.
STRICT 빌드(`-Wall -Wextra -Werror`) 통과 확인됨.

---

## 3. OCRT 본체 연결 (3단계)

`rt_water_iop.h`에 문서화된 성분 추가 절차를 따른다. TSM과 동일한 패턴이다.

### 3.1 성분 필드 추가

`rt_iop_components_t`(rt_water_iop.h)에 두 필드를 추가한다.

```c
typedef struct {
    rt_iop_t pure_water;
    rt_iop_t cdom;
    rt_iop_t pigment;      /* 기존 CCRR (선택적 대체) */
    rt_iop_t mineral;      /* 기존 CCRR (TSM 모듈로 대체) */
    rt_iop_t eap_phyto;    /* 신규: EAP 식물플랑크톤 */
    rt_iop_t detritus;     /* 신규: 유기 detritus */
} rt_iop_components_t;
```

기존 CCRR `pigment`를 EAP로 대체할지, 병행할지는 선택이다. EAP가 편광 위상함수를
제공하므로 대기보정 논문 목적에는 EAP 사용을 권장한다.

### 3.2 평가 분기 추가

`rt_iop_components_eval_all()`(또는 해당 집계 함수)에 두 성분 호출을 추가한다.

```c
#include "rt_iop_organic.h"

/* EAP 식물플랑크톤: .mie 경로는 Brewin 계급 옵션(pico/nano/micro) */
rt_iop_eap_phyto_eval(lambda_nm, chl_mg_m3, eap_mie_path, &comp->eap_phyto);

/* 유기 detritus: a_d440은 Chl covary(지역보정), S_d는 옵션(≤0→기본 0.0109) */
rt_iop_detritus_eval(lambda_nm, chl_mg_m3, a_d440_chl, S_d_nm_inv,
                     detritus_mie_path, &comp->detritus);
```

`rt_iop_total()`의 단순 합산(a_total=Σa, b_total=Σb, bb_total=Σbb)에 두 성분이
자동 포함되도록 집계 함수에 더한다.

### 3.3 파라미터 연결

CLI 또는 입력 파일에서 다음을 연결한다.

| 파라미터 | 의미 | 기본/범위 |
|---|---|---|
| `chl_mg_m3` | 엽록소 농도 | 입력 |
| `eap_mie_path` | EAP 그룹 .mie 경로 | pico/nano/micro 옵션 |
| `a_d440_chl` | detritus 흡광 절대(440nm) | Chl covary, 지역보정 |
| `S_d_nm_inv` | detritus 흡광 기울기 | ≤0→0.0109, 범위 0.0024-0.017 |

`a_d440_chl`은 Chl covary로 결정한다(방향 B 절대 스케일). 정확한 covary 계수는
지역 보정 대상이다(SPEC 문서 5.3 참조).

---

## 4. 위상함수(Mueller) 병합 — 별도 단계

현재 `rt_iop_organic`은 스칼라 IOP(a/b/bb)만 산출한다. 편광 대기보정의 핵심인
P11/P12/P33 Mueller 병합은 OCRT의 기존 phase moments 인프라
(`rt_iop_ccrr_phase_moments_load`)에 `.mie`의 P11/P12/P33를 연결하는 별도 단계다.

병합은 산란 가중 평균이다(방향 B 절대 스케일 연결, 결정 1·2).

```
P_total(λ,μ) = Σ_i [ b_i(λ) · P_i(λ,μ) ] / b_total(λ)
```

각 성분의 b_i는 3.2에서 산출된 값을 쓴다. EAP·detritus·TSM·순수해수의 P_i를
b 가중으로 합친다. 이 단계는 벡터 솔버(--vector, --surface ocean 경로)와 연동된다.

---

## 5. 옵션 게이트 (권장)

수치 솔버 옵션은 `OCRT_ADVANCED` 뒤에, 진단 옵션은 `OCRT_DEBUG` 뒤에 게이트한다.
detritus `S_d` 옵션은 일반 사용자에게 노출하되, 범위 밖 값은 경고를 출력한다
(구현에 포함됨). 기본 편광 모드는 벡터로 유지한다.

---

## 6. 검증

통합 후 커밋 게이트를 통과해야 한다.
- 직전 버전과 비트 동일 확인(유기 성분 미사용 시)
- Tier-0 앵커(rrs0minus) 재현
- STRICT 빌드 유지
- EAP a*(443)≈0.034(Diatoms), detritus 흡광 지수감쇠(S_d=0.0109) 확인

독립 검증(test_organic)에서 확인된 값:
- Huot bbp(550,Chl=1)=2.267e-3, f_ph(Chl=1)=0.035
- EAP Diatoms a(443,Chl=1)=0.0342, bb/b(550)=0.0036
- detritus a(443,Chl=1,a_d440=0.01)=0.0097, bb/b(550)=0.0059
