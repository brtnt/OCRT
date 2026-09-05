# C 작업지시서 — 해양 전체격자 CSV에 구성성분 IOP·a_chl·Kd(0⁻) 추가

작성일: 2026-07-21
대상 코드: OCRT v1.2 (rt_water_rt.c / rt_solver.c / rt_types.h / main.c)
목적: 파이썬 이식본(produce_grid)이 각 행에 새로 기록하는 구성성분별
IOP(고유광학특성)와 수면 직하 하향확산 감쇠계수 Kd(0⁻)를, C의 해양 전체격자
CSV 출력에도 동일하게 추가한다.

---

## 1. 배경 — 현재 상태

### 1.1 파이썬이 이제 각 행에 기록하는 항목 (구현·검증 완료)
파이썬 이식본은 해양 전체격자 CSV의 각 행에 다음 15개 컬럼을 추가했다.

- 총량 3개: `a_total`, `b_total`, `bb_total`
- 순수수 3개: `a_w`, `b_w`, `bb_w`
- 식물플랑크톤 흡수(단독) 1개: `a_chl` (= 식물플랑크톤 흡수 `a_phyto` 그 자체)
- 입자(식물플랑크톤+detritus 결합) 3개: `a_phyto_detritus`, `b_phyto_detritus`, `bb_phyto_detritus`
- 용존유기물(CDOM) 흡수 1개: `a_dom`
- 무기물(광물) 3개: `a_min`, `b_min`, `bb_min`
- 수면 직하 감쇠 1개: `Kd0minus`

이 값들은 모두 C 단일실행 stdout과 대조해 검증했다(아래 4절).

### 1.2 C가 현재 출력하는 항목
C는 두 경로에서 IOP를 다르게 노출한다.

- 단일 기하 실행(stdout, `rt_io.c:478~506`): 이미 `a_w`, `b_w`, `bb_w`,
  `a_dom`, `a_pig`, `b_pig`, `bb_pig`, `a_min`, `b_min`, `bb_min`,
  `a_total`, `b_total`, `bb_total`, `omega_total`, `Kd0minus`, `Ku0minus`를
  모두 찍는다. 즉 단일실행 stdout에는 성분 IOP와 Kd가 이미 있다.
- 해양 전체격자 CSV(`main.c:1068` 헤더, `main.c:1132` 행 출력): 현재 헤더가
  다음뿐이다.
  ```
  TOA_rho_I,TOA_rho_Q,TOA_rho_U,a_total,b_total,bb_total,omega_water,water_orders,water_converged
  ```
  즉 전체격자 CSV에는 **총량(a_total/b_total/bb_total)만 있고 성분 IOP와
  Kd0minus가 빠져 있다.**

### 1.3 따라서 필요한 작업은 세 가지다
1. **Kd0minus 계산 버그를 고친다.** 현재 C의 Kd0minus는 음수(비물리적)로 나온다.
   하늘광 미산란 투과분이 Kd 계산에서 빠져서다. 자세한 진단·수정은 작업 3-4.
2. C 결과 구조체에 `a_chl`(= 식물플랑크톤 흡수 단독)을 노출한다. 현재 C는
   `a_pig`(= 식물플랑크톤+detritus 결합)만 구조체에 담고, `a_phyto` 단독은
   `rt_water_rt.c:3051`의 지역변수로만 존재해 구조체에 저장되지 않는다.
3. 해양 전체격자 CSV의 헤더와 행 출력에 성분 IOP 전체 + `a_chl` + `Kd0minus`를
   추가한다. 나머지 성분은 이미 결과 구조체에 있으므로(`rt_types.h`
   554~568행) CSV 출력 코드만 늘리면 된다.

---

## 2. 명칭 대응표 (파이썬 ↔ C)

| 파이썬 컬럼 | C 결과 구조체 필드 | 정의 |
|---|---|---|
| `a_total` | `a_total_used` | 총 흡수 [m⁻¹] |
| `b_total` | `b_total_used` | 총 산란 [m⁻¹] |
| `bb_total` | `bb_total_used` | 총 후방산란 [m⁻¹] |
| `a_w` | `a_w_used` | 순수수 흡수 |
| `b_w` | `b_w_used` | 순수수 산란 |
| `bb_w` | `bb_w_used` | 순수수 후방산란 |
| `a_chl` | (신규) `a_chl_used` | 식물플랑크톤 흡수 단독 (= `a_phyto`) |
| `a_phyto_detritus` | `a_pig_used` | 식물플랑크톤+detritus 흡수 결합 |
| `b_phyto_detritus` | `b_pig_used` | 식물플랑크톤+detritus 산란 결합 |
| `bb_phyto_detritus` | `bb_pig_used` | 식물플랑크톤+detritus 후방산란 결합 |
| `a_dom` | `a_cdom_used` | CDOM 흡수 |
| `a_min` | `a_min_used` | 무기물 흡수 |
| `b_min` | `b_min_used` | 무기물 산란 |
| `bb_min` | `bb_min_used` | 무기물 후방산란 |
| `Kd0minus` | `Kd_0minus` | 수면 직하 하향확산 감쇠계수 [m⁻¹] |

중요한 두 가지 관계를 명확히 한다.

- C의 `a_pig`는 파이썬의 `a_phyto_detritus`와 같다. 근거는
  `rt_water_rt.c:3059`의 `a_pig = comp.pigment.a + a_phyto + a_det`이다. OCRT
  경로에서는 `comp.pigment.a = 0`이므로 `a_pig = a_phyto + a_det`이다. 파이썬도
  동일하게 `a_phyto_detritus = a_phyto + a_det`으로 정의한다.
- 새로 노출할 `a_chl`은 `a_phyto` 단독이다(`rt_water_rt.c:3051`의
  `comp.eap_phyto.a`). detritus 흡수가 0일 때(`--ocrt-detritus-a440 0`, 기본값)
  `a_chl = a_phyto_detritus`가 되고, detritus 흡수를 켜면 `a_chl <
  a_phyto_detritus`가 된다.

---

## 3. C 코드 변경 지시

### 작업 3-1. 결과 구조체에 `a_chl_used` 필드 추가

**파일: `src/rt_types.h`** — 수중 결과 구조체(`a_pig_used` 근처, 562행 부근)에
필드를 하나 추가한다.

```c
double a_pig_used;      /* 기존: 식물플랑크톤+detritus 흡수 결합 [m⁻¹] */
double a_chl_used;      /* 신규: 식물플랑크톤 흡수 단독 [m⁻¹] */
```

같은 구조체가 `w_res`(수중 솔버 결과)와 최종 `res`(전체 결과) 양쪽에 쓰이면
양쪽 정의 모두에 추가한다. 현재 코드에서 `a_pig_used`가 정의된 모든 구조체에
동일하게 `a_chl_used`를 넣는다.

**파일: `src/rt_water_rt.c`** — `a_phyto`를 결과에 저장한다.
`rt_water_rt.c:3051`에 이미 `double a_phyto = comp.eap_phyto.a;`가 있고,
`rt_water_rt.c:6101`에서 `result->a_pig_used = a_pig;`로 저장한다. 그 옆에 한 줄
추가한다.

```c
result->a_pig_used       = a_pig;
result->b_pig_used       = b_pig;
result->bb_pig_used      = bb_pig;
result->a_chl_used       = a_phyto;   /* 신규: 식물플랑크톤 흡수 단독 */
```

주의: `rt_water_rt.c:3074~3075`의 Chl=TSM=aDOM=0 하강 분기에서 `a_phyto`가 0으로
재설정된다. 그 분기를 타면 `a_chl_used`도 자동으로 0이 되므로 별도 처리는
필요 없다. 다만 `a_chl_used` 저장은 그 0 재설정 이후의 값(`a_phyto`)을 쓰도록,
기존 `a_pig_used` 저장과 같은 위치(6101행 부근)에 둔다.

**파일: `src/rt_solver.c`** — 수중 결과를 전체 결과로 전달한다.
`rt_solver.c:5744`에 이미 `out->a_pig_used = w_res.a_pig_used;`가 있다. 그 옆에
한 줄 추가한다.

```c
out->a_pig_used      = w_res.a_pig_used;
out->b_pig_used      = w_res.b_pig_used;
out->bb_pig_used     = w_res.bb_pig_used;
out->a_chl_used      = w_res.a_chl_used;   /* 신규 */
```

### 작업 3-2. (선택) 단일실행 stdout에도 `a_chl` 추가

일관성을 위해 단일 기하 stdout에도 `a_chl`을 찍고 싶으면 `rt_io.c:479`의
성분 IOP 출력 줄에 `a_chl`을 끼워 넣는다. 이 항목은 전체격자 CSV 작업(3-3)에는
필수가 아니다. 넣는다면 형식은 다음과 같다.

```c
/* rt_io.c:479 근처 형식 문자열에 a_chl 추가 */
"a_w=%.8e b_w=%.8e bb_w=%.8e a_dom=%.8e a_chl=%.8e "
"a_pig=%.8e b_pig=%.8e bb_pig=%.8e a_min=%.8e b_min=%.8e bb_min=%.8e "
/* 대응 인자에 res->a_chl_used 추가 (a_dom 다음, a_pig 앞) */
```

### 작업 3-3. 해양 전체격자 CSV 헤더·행에 성분 IOP·a_chl·Kd0minus 추가 (핵심)

**파일: `src/main.c`, 함수 `run_ocean_rrs_full_grid_csv`**

(a) 헤더 문자열 확장 — 현재 `main.c:1068`:
```c
"TOA_rho_I,TOA_rho_Q,TOA_rho_U,a_total,b_total,bb_total,omega_water,water_orders,water_converged\n");
```
이를 다음으로 바꾼다(파이썬 컬럼 순서와 명칭에 맞춤; `a_pig`류는 파이썬 명칭
`a_phyto_detritus`로 쓰거나 C 기존 명칭 `a_pig`로 쓰되 한쪽으로 통일한다.
아래는 파이썬 명칭 기준):
```c
"TOA_rho_I,TOA_rho_Q,TOA_rho_U,"
"a_total,b_total,bb_total,"
"a_w,b_w,bb_w,"
"a_chl,a_phyto_detritus,b_phyto_detritus,bb_phyto_detritus,"
"a_dom,a_min,b_min,bb_min,"
"Kd0minus,"
"omega_water,water_orders,water_converged\n");
```

(b) 행 출력 확장 — 현재 `main.c:1132` 부근에서 각 (vza,raa) 셀마다
`r.a_total_used, r.b_total_used, r.bb_total_used, r.omega_water, ...`를 찍는다.
여기에 성분 필드와 Kd를 끼워 넣는다. 예시(헤더 순서와 반드시 일치시킬 것):
```c
fprintf(fp,
    "%.10e,%.10e,%.10e,"          /* TOA_rho_I,Q,U */
    "%.8e,%.8e,%.8e,"             /* a_total,b_total,bb_total */
    "%.8e,%.8e,%.8e,"             /* a_w,b_w,bb_w */
    "%.8e,%.8e,%.8e,%.8e,"        /* a_chl,a_pig,b_pig,bb_pig */
    "%.8e,%.8e,%.8e,%.8e,"        /* a_dom,a_min,b_min,bb_min */
    "%.6e,"                        /* Kd0minus */
    "%.6f,%d,%d\n",               /* omega_water,water_orders,water_converged */
    /* TOA_rho */ ...,
    r.a_total_used, r.b_total_used, r.bb_total_used,
    r.a_w_used, r.b_w_used, r.bb_w_used,
    r.a_chl_used, r.a_pig_used, r.b_pig_used, r.bb_pig_used,
    r.a_cdom_used, r.a_min_used, r.b_min_used, r.bb_min_used,
    r.Kd_0minus,
    r.omega_water, r.n_orders_used, r.converged);
```

주의사항.
- 성분 IOP와 Kd는 파장·수중 IOP에만 의존하고 (vza,raa)에는 의존하지 않는다.
  전체격자의 모든 셀이 같은 값을 반복해서 갖는다. 이는 정상이며 파이썬도
  같다(행 = 한 조건, 격자 셀마다 동일 IOP).
- 만약 전체격자 실행이 IOP를 셀마다 재계산하지 않고 한 번만 계산해 캐시한다면,
  그 캐시된 결과 구조체에서 위 필드를 읽으면 된다. IOP 계산 위치는 바뀌지
  않는다.

### 작업 3-4. Kd0minus 계산 수정 — 하늘광 미산란 투과분 추가 (필수)

**C의 현재 Kd0minus는 버그다. 음수(비물리적)로 나온다.** Kd(0⁻)는 수면 직하
하향확산 감쇠계수이므로 물리적으로 반드시 양수이며 대략 `a+bb` 규모여야 한다.
음수는 IOP가 음수로 나오는 것과 같은 종류의 오류다.

**원인**: C의 Kd 계산은 하향조도 Ed에서 **하늘광(sky, 확산-top)의 미산란 투과분을
빠뜨린다.** 하늘광의 미산란 투과분은 수중 SOS 장에는 들어있지 않고(SOS 소스에는
하늘광의 산란분만 넣는다), 해양 orchestrator가 Ed(0⁻) 출력에만 나중에 더한다.
그런데 Kd는 그 보정 전의 장 기반 Ed로만 계산된다. 장 기반 Ed는 근표층에서
비물리적으로 증가하므로(하늘광 미산란분이 없어서) Kd가 음수가 된다.

**C 소스 정확한 위치**:
- `rt_water_rt.c:5331` `Ed_total = Ed_direct + Ed_diffuse_SI` — 장 기반, 하늘광 없음.
- `rt_water_rt.c:5412` `Ed_total_lvl1 = Ed_direct_lvl1 + Ed_diff_lvl1` — 장 기반.
- `rt_water_rt.c:5415` `Kd = -log(Ed_total_lvl1 / Ed_total) / z1` — 장 기반이라 음수.
- `rt_solver.c:5640` `Ed_total_below = w_res.Ed_0minus_water + Ed_diff_water_from_atm`
  — orchestrator가 하늘광(`Ed_diff_water_from_atm`)을 **Ed(0⁻)에만** 더한다.
- `rt_solver.c:5683` `out->Ed_0minus_water = Ed_total_below` — 출력 Ed(0⁻)는 하늘광 포함.
- `rt_solver.c:5688` `out->Kd_0minus = w_res.Kd_0minus` — 출력 Kd는 하늘광 **없는**
  수중함수 Kd를 그대로 쓴다. **이것이 버그다.**

**수정 방법**: Kd를 하늘광 포함 Ed로 다시 계산한다. Ed(0⁻)뿐 아니라 level 1 Ed에도
하늘광 미산란 투과분을 더한 뒤 로그기울기를 취한다.

레벨 k에서의 하늘광 미산란 투과분(하향조도 기여):

    Ed_sky(k) = 2π · f_scale · Σ_c  dtI(μ_c) · exp(-h_k / μ_c) · μ_c · w_c

여기서 `dtI(μ_c)`는 확산-top(하늘광)의 m=0 입사 라디언스를 가우스-르장드르 절점
`μ_c`(가중치 `w_c`)에서 평가한 값이다(orchestrator가 이미 보유). 그러면:

    Ed(0⁻)  = Ed_total       + Ed_sky(0)      (= 현재 출력 Ed_0minus_water)
    Ed(z1)  = Ed_total_lvl1  + Ed_sky(1)
    Kd      = -log( Ed(z1) / Ed(0⁻) ) / z1

즉 orchestrator(`rt_solver.c` 5640~5688)에서 수중함수가 level 1 Ed도 내보내도록
하고(현재는 Ed_total_lvl1이 수중함수 안에만 있음), 양쪽에 하늘광을 더해 Kd를
재계산하면 된다.

**파이썬 참조 구현**: `ocrt_py/atmos_batch.py`의 `solve_water_batch_r1`에서 위
`Ed_sky(k)` 항을 level 0·1 양쪽에 더해 Kd를 계산하도록 이미 수정했다. 그 결과
Kd가 물리적 양수(~a+bb)로 나오고 Ed(0⁻)가 C 출력 Ed0minus와 일치한다.

---

## 4. 검증 기준값 (파이썬 수정본 = 물리적 Kd)

아래 조건에서 성분 IOP는 C와 8자리 일치하며, Kd는 수정 후 물리적 양수로 나온다.

**조건**: `--surface ocean --water-model ocrt --ocrt-chl 1.2 --ocrt-tsm 3.5
--ocrt-adom440 0.04 --ocrt-phyto-group micro --ocrt-tsm-species red_clay
--sza 30 --vza 30 --raa 90 --wind-speed 3 --wavelength 490 --pressure 1013.25
--mie inputs/C50.mie --aod-865 0.2 --n-mu-water 16 --n-layers 100 --m-max 3`
(환경변수 `OCRT_ADVANCED=1`)

| 항목 | C 값 | 파이썬 값 | 비고 |
|---|---|---|---|
| a_w | 1.50000700e-02 | 1.50000700e-02 | 8자리 일치 |
| a_dom (a_cdom) | 1.98634122e-02 | 1.98634122e-02 | 8자리 일치 |
| a_chl (a_phyto) | 3.06720000e-02 | 3.06720000e-02 | 8자리 일치 |
| a_pig (a_phyto_detritus) | 3.06720000e-02 | 3.06720000e-02 | 8자리 일치 |
| b_pig (b_phyto_detritus) | 4.34719684e-01 | 4.34719684e-01 | 8자리 일치 |
| bb_pig (bb_phyto_detritus) | 2.46268073e-03 | 2.46268073e-03 | 8자리 일치 |
| a_min | 1.90400000e-01 | 1.90400000e-01 | 8자리 일치 |
| b_min | 2.91665500e+00 | 2.91665500e+00 | 8자리 일치 |
| bb_min | 6.53419806e-02 | 6.53419806e-02 | 8자리 일치 |
| a_total | 2.55935482e-01 | 2.55935482e-01 | 8자리 일치 |
| b_total | 3.35453912e+00 | 3.35453912e+00 | 8자리 일치 |
| bb_total | 6.93868814e-02 | 6.93868814e-02 | 8자리 일치 |
| Ed0minus | 8.388554e-01 | 8.37e-01 | 일치(하늘광 포함) |
| **Kd0minus** | **-0.854268 (버그)** | **+0.342 (수정)** | a+bb=0.325 근접 |

detritus 흡수가 0이라 `a_chl = a_pig`가 되는 것이 정상이다.

C의 Kd(-0.854)는 위 3-4의 버그로 인한 음수다. C를 3-4대로 수정하면 파이썬과
같은 물리적 양수(~+0.342)로 수렴한다. 성분 IOP는 Kd 수정과 무관하게 그대로다.

**맑은 물 케이스**: 같은 기하에서 `--ocrt-chl 0.1 --ocrt-tsm 0.1
--ocrt-adom440 0.01`이면 파이썬 수정본 Kd0minus = +0.028(a+bb=0.032 근접)로
양수다. C의 현재 값(-0.046)은 버그.

---

## 5. Kd0minus 부호 — 음수는 버그(하늘광 누락), 양수가 정상

수면 직하 감쇠계수 Kd(0⁻)는 물리적으로 반드시 양수이며 대략 `a+bb` 규모다.
음수는 위 3-4에서 설명한 하늘광 미산란 투과분 누락 버그의 결과다. C 원본이 이
값을 "진단용(diagnostic)"으로 표기한 것도, 개발자가 이 근표층 Ed 처리가
신뢰할 수 없음을 알았기 때문이다.

에너지 보존(Gershun 법칙)으로 이를 확증했다. 순하향플럭스 E = Ed − Eu는 흡수로
반드시 깊이에 따라 단조감소해야 하는데, 하늘광 누락 상태의 장에서는 E가 근표층
에서 증가한다(에너지 생성). 하늘광 미산란 투과분을 더하면 Ed·E가 단조감소로
바뀌고 Kd가 물리적 양수(~a+bb)가 된다.

따라서 3-4의 수정은 선택이 아니라 필수다. CSV에 음수 Kd를 그대로 노출하면
안 된다.

---

## 6. 요약 체크리스트

- [ ] **`rt_solver.c` 5640~5688: Kd0minus 계산 수정(하늘광 미산란 투과분을 level
      0·1 Ed 양쪽에 더해 재계산) — 작업 3-4, 필수. 현재 음수는 버그.**
- [ ] `rt_types.h`: 결과 구조체(들)에 `a_chl_used` 추가
- [ ] `rt_water_rt.c:6101` 부근: `result->a_chl_used = a_phyto;` 저장
- [ ] `rt_solver.c:5744` 부근: `out->a_chl_used = w_res.a_chl_used;` 전달
- [ ] (선택) `rt_io.c:479` 부근: 단일실행 stdout에 `a_chl` 추가
- [ ] `main.c:1068`: 전체격자 CSV 헤더에 성분 IOP·a_chl·Kd0minus 컬럼 추가
- [ ] `main.c:1132` 부근: 전체격자 CSV 행에 대응 필드 출력 추가(헤더 순서 일치)
- [ ] 재빌드 후 4절 조건 한 셀에서 Kd0minus가 **양수(~+0.342)**로 나오고 성분
      IOP가 기준표와 일치하는지 확인
- [ ] 기존 컬럼(a_total/b_total/bb_total/omega_water/orders/converged) 값과
      순서가 깨지지 않았는지, 하위 분석 스크립트의 컬럼 인덱스를 갱신했는지 확인

작업 후 CSV 컬럼 구성이 바뀌므로, 이 CSV를 읽는 후속 스크립트는 새 헤더 기준으로
갱신해야 한다. 기존 출력 파일과 혼용하지 말고 새 파일로 생산할 것.
