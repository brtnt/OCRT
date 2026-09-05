# OCRT v1.2 `--water-m-max` mode-bound 결함 수정 지시서

**Issue ID:** OCRT-RT-C3D-EXTBOTTOM-001  
**작성일:** 2026-07-20  
**대상 기준 소스:** `OCRT_v1.2_SNF_PLUS_ACCURT_AHMAD2010_80_AEROSOL_MODELS_2026-07-19`  
**우선순위:** **P0 / release blocker**  
**결함 유형:** silent numerical corruption, heap out-of-bounds read 가능성, 잘못된 `conv=1` 보고

---

## 1. 수정 목적

대기 Fourier 상한과 수중 Fourier 상한이 다를 때, 수중 water-leaving bottom-source 배열보다 더 높은 대기 Fourier mode를 읽는 결함을 제거한다.

대표적으로 aerosol-on 조건은 대기 `fourier_m_max=16`으로 자동 설정된다. 이 상태에서 사용자가 `--water-m-max 2`, `4`, `8`처럼 더 작은 값을 지정하면, 수중장 자체는 유한하게 계산되지만 TOA I/Q/U가 폭주하거나 불안정해진다.

이 결함은 종료코드 0과 `conv=1`을 반환하므로 단순 실패가 아니라 **조용한 결과 오염**이다.

---

## 2. 관측된 증상

공통 진단 조건은 atmosphere `n_mu=48`, water `n_mu=48`, aerosol-on, coupled ocean이다.

| `water_m_max` | TOA I | TOA Q | TOA U | 판정 |
|---:|---:|---:|---:|---|
| 2 | -1.355e+244 | 8.295e+245 | 3.505e+245 | CORRUPT_TOA |
| 4 | 5.392e+244 | 2.806e+243 | 3.171e+228 | CORRUPT_TOA |
| 8 | 0.096691082 | 0.032430329 | 1.447475770 | UNSTABLE |
| 16 | 0.096672745 | 0.032480902 | -0.011013802 | NORMAL |
| 30 | 0.096672745 | 0.032480902 | -0.011013802 | NORMAL |
| 48 | 0.096672745 | 0.032480902 | -0.011013802 | NORMAL |
| 64 | 0.096672745 | 0.032480902 | -0.011013802 | NORMAL |

핵심 특징:

- `m=2`, `m=4`: TOA가 약 `10^244–10^245`로 폭주한다.
- `m=8`: 값은 유한하지만 TOA U가 비물리적으로 커지고 반복 안정성이 없다.
- `m>=16`: 정상 범위로 복귀한다.
- 수중 `Rrs/rrs` 장은 같은 실행에서 유한하다.
- 실행은 `PASS`, `orders=7`, `conv=1`로 종료한다.

![관측된 stability defect](evidence/nmu48_water_fourier_stability_defect.png)

---

## 3. 코드상 직접 원인

### 3.1 대기 mode 수

`src/main.c:2992–3005`에서 aerosol-on이면 사용자가 `--m-max`를 지정하지 않은 경우 대기 `opts.fourier_m_max=16`이 된다.

### 3.2 C3d direct-water bottom source의 실제 배열 크기

`src/rt_solver.c:4541–4548`:

```c
int mmax = opts->fourier_m_max;
if (mmax > w_res.m_max_filled) mmax = w_res.m_max_filled;
size_t na = (size_t)(mmax + 1) * (size_t)n_mu_atm;
double *wlI = (double*)calloc(na, sizeof(double));
double *wlQ = (double*)calloc(na, sizeof(double));
double *wlU = (double*)calloc(na, sizeof(double));
```

즉 배열은 `min(atm_m_max, water_m_max)+1`개 mode만 가진다.

### 3.3 대기 pass-2는 원래 atmospheric `m_max`를 유지

`src/rt_solver.c:4580–4586`:

```c
rt_options_t opts_p2 = *opts;
...
opts_p2.ext_bottom_per_m_I = wlI;
opts_p2.ext_bottom_per_m_Q = wlQ;
opts_p2.ext_bottom_per_m_U = wlU;
```

`opts_p2.fourier_m_max`는 줄이지 않으므로 aerosol-on에서는 16이다.

### 3.4 대기 solver가 배열 끝을 넘어 mode slice를 생성

`src/rt_solver.c:675`는 `m=0..m_max`를 순회하고, `740–743`은 mode-bound 정보 없이 다음 주소를 만든다.

```c
sos_opt.ext_bottom_I = &opts->ext_bottom_per_m_I[(size_t)m * (size_t)n_mu];
sos_opt.ext_bottom_Q = &opts->ext_bottom_per_m_Q[(size_t)m * (size_t)n_mu];
sos_opt.ext_bottom_U = &opts->ext_bottom_per_m_U[(size_t)m * (size_t)n_mu];
```

따라서 `m > water source m_max`이면 heap 범위 밖을 읽는다.

### 3.5 threshold가 16인 이유

- aerosol-on atmospheric default: `m_max=16`
- `water_m_max < 16`: pass-2가 short array 뒤의 메모리를 읽음
- `water_m_max >= 16`: 배열이 atmosphere loop 전체를 포함하므로 증상이 사라짐

관측 threshold와 코드 구조가 정확히 일치한다.

### 3.6 동일 파일의 pass-3에는 이미 올바른 패턴이 존재

`src/rt_solver.c:5547–5569`의 skylight pass는 다음과 같이 full atmospheric size를 `calloc`한다.

```c
int mmax_atm = opts->fourier_m_max;
size_t na = (size_t)(mmax_atm + 1) * (size_t)n_mu_atm;
double *slI = (double*)calloc(na, sizeof(double));
```

주석도 `m>=1`을 0으로 유지해 `OOB·garbage`를 방지한다고 명시한다. Direct-water C3d pass에 같은 보호가 누락된 상태다.

---

## 4. 필수 수정 사항

### 4.1 즉시 수정: C3d 배열을 atmospheric mode 수로 zero-pad

`rt_solver.c`의 direct-water C3d pass에서 아래 두 상한을 분리한다.

```c
const int atm_mmax = opts->fourier_m_max;
const int water_src_mmax =
    (w_res.m_max_filled < atm_mmax) ? w_res.m_max_filled : atm_mmax;

size_t na = (size_t)(atm_mmax + 1) * (size_t)n_mu_atm;
double *wlI = calloc(na, sizeof(*wlI));
double *wlQ = calloc(na, sizeof(*wlQ));
double *wlU = calloc(na, sizeof(*wlU));

/* Only 0..water_src_mmax are populated. Higher atmospheric modes remain zero. */
rt_air_water_couple_water_to_atm(
    ...,
    water_src_mmax,
    ...,
    wlI, wlQ, wlU);
```

이 변경은 pass-3의 기존 구현 방식과 동일하다.

### 4.2 API hardening: external bottom-source의 shape metadata 추가

현재 `rt_options_t`에는 포인터만 있고 mode 수와 stride 정보가 없다. 다음 중 하나를 적용한다.

#### 권장 구조

```c
typedef struct {
    const double *I;
    const double *Q;
    const double *U;
    int n_mu;
    int m_max;
} rt_fourier_bottom_source_t;
```

`rt_options_t`에는 이 구조체를 넣는다.

#### 최소 변경 구조

```c
const double *ext_bottom_per_m_I;
const double *ext_bottom_per_m_Q;
const double *ext_bottom_per_m_U;
int ext_bottom_n_mu;
int ext_bottom_m_max;
```

기본값은 `ext_bottom_n_mu=0`, `ext_bottom_m_max=-1`로 한다.

### 4.3 atmosphere solver에서 bounds를 강제

mode loop마다 포인터를 먼저 NULL로 초기화하고, 유효한 mode에서만 slice를 연결한다.

```c
sos_opt.ext_bottom_I = NULL;
sos_opt.ext_bottom_Q = NULL;
sos_opt.ext_bottom_U = NULL;

if (have_external_bottom_source) {
    if (opts->ext_bottom_n_mu != n_mu || opts->ext_bottom_m_max < 0) {
        return RT_ERR_BOTTOM_SOURCE_SHAPE;
    }
    if (m <= opts->ext_bottom_m_max) {
        size_t off = (size_t)m * (size_t)n_mu;
        sos_opt.ext_bottom_I = opts->ext_bottom_per_m_I + off;
        sos_opt.ext_bottom_Q = opts->ext_bottom_per_m_Q + off;
        sos_opt.ext_bottom_U = opts->ext_bottom_per_m_U + off;
    }
}
```

`bottom_source_only=1`이고 source mode 범위를 넘으면 해당 mode는 zero source로 처리한다. Fourier mode는 서로 독립이므로 이 동작이 물리적으로 맞다.

### 4.4 모든 caller에 shape를 설정

최소 다음 경로를 전수 확인한다.

- C3d direct-water pass (`opts_p2`)
- C3-full skylight pass (`opts_p3`)
- S7/S7b atmosphere grid cache
- D3 operator/diagnostic paths
- 향후 `ext_bottom_per_m_*`를 설정하는 모든 call site

### 4.5 S7b source fingerprint의 shape 일치

현재 S7b fingerprint는 local `mmax+1`까지만 hash한다. 실제 solver가 소비하는 mode 범위와 fingerprint 범위가 동일해야 한다.

- full atmospheric zero-padding을 사용할 경우 `(atm_mmax+1)*n_mu_atm` 전체를 hash
- metadata-bound 방식을 사용할 경우 `ext_bottom_m_max`와 `ext_bottom_n_mu`를 cache key에 포함

### 4.6 임시 release guard

정식 fix가 병합되기 전 배포 브랜치에서는 다음 조합을 즉시 거부해 silent corruption을 막는다.

```text
coupled ocean + atmosphere present
AND explicit water_m_max < atmospheric fourier_m_max
```

이 guard는 임시조치다. 최종 해결을 `water_m_max >= atm_m_max` 강제나 CLI 무시로 대체하면 안 된다.

---

## 5. 금지되는 수정 방식

다음은 결함을 숨길 뿐 해결하지 못한다.

1. `--water-m-max`를 무조건 16 이상으로 올림
2. 사용자가 지정한 low mode를 조용히 무시
3. TOA 출력에 `isfinite()`만 추가
4. 폭주값만 0으로 clamp
5. stable `m>=16` 경로의 물리 결과를 변경

낮은 water Fourier cap은 진단·수렴 시험에서 유효한 설정이어야 하므로 shape mismatch를 정상 처리해야 한다.

---

## 6. 필수 회귀시험

### 6.1 sanitizer test

ASan/UBSan 빌드로 aerosol-on coupled-ocean case를 실행한다.

```text
atmosphere n_mu = 48
water n_mu = 48
atmosphere m_max = 16 (aerosol auto)
water_m_max = 2, 4, 8, 16, 30, 48, 64
```

요구사항:

- heap-buffer-overflow 0건
- invalid read/write 0건
- UAF/double-free 0건
- 모든 실행 return code 0 또는 명시적인 입력오류

### 6.2 numerical stability matrix

각 `water_m_max`에 대해 다음을 기록한다.

```text
TOA I/Q/U
Lu(0+)/Lu(0−) I/Q/U
Rrs(0+)/rrs(0−) I/Q/U
orders, conv
wall time
```

판정:

- 모든 출력 finite
- `abs(TOA I/Q/U) < 100`의 broad corruption guard 통과
- `m=2,4,8`에서 `10^N` 폭주 없음
- `m>=16`은 수정 전 정상 baseline과 bit-identical 또는 `max_abs<=1e-12`

### 6.3 zero-mode padding unit test

대기 `m_max=16`, external bottom source `m_max=2`인 bottom-only solve를 직접 구성한다.

- mode 0..2: 입력 source가 전달됨
- mode 3..16: exact zero source
- mode 3..16 결과가 zero
- sanitizer clean

### 6.4 기존 production regression

다음은 변경되면 안 된다.

- auto water `m_max=30` full-vector 대표 cell
- Item 6 golden 5-cell
- Item 6 전체 162-cell baseline
- Rayleigh-only pure-ocean 6-band regression
- atmosphere-off water solve
- S7/S7b cache ON/OFF 결과

### 6.5 cache determinism

같은 case를 cold/warm 및 반복 실행해:

- output hash 동일
- cache key 동일
- low-mode 실행 후 high-mode 실행 및 역순 실행 결과 동일

---

## 7. 완료 판정 기준

| 항목 | 수용 기준 |
|---|---|
| Low `water_m_max` TOA | finite, 비폭주 |
| ASan/UBSan | 0 errors |
| Stable `m>=16` 결과 | bit-identical 또는 `max_abs<=1e-12` |
| Auto `m=30` Item 6 | 162/162 baseline 유지 |
| `conv` 신뢰성 | corrupt output에 `conv=1` 발생하지 않음 |
| Cache | 실행순서 독립, hash 재현 |
| 성능 | 정상 `m=30` 경로 wall time 회귀 1% 이하 |
| 문서 | CHANGELOG 및 CLI 설명 갱신 |

---

## 8. 메인 업데이트 세션 제출물

수정 세션은 다음을 반환한다.

1. Git commit hash 및 patch
2. 수정 소스 archive
3. clean release build와 ASan/UBSan build 로그
4. `water_m_max=2,4,8,16,30,48,64` 결과 CSV
5. 수정 전/후 stability 산포도
6. Item 6 golden 5 및 162-cell 회귀 결과
7. binary SHA-256와 source tree hash
8. 변경된 source line 목록

---

## 9. 메인 세션에 전달할 짧은 지시문

> OCRT v1.2의 coupled-ocean C3d atmosphere bottom-source pass에서 `opts->fourier_m_max`보다 짧은 `wlI/Q/U` mode 배열을 전달해 heap OOB read와 silent TOA corruption이 발생한다. Aerosol-on atmospheric m_max=16에서 explicit `--water-m-max 2/4`는 TOA가 10^244–10^245로 폭주하고 `m=8`은 불안정하며, `m>=16`은 정상이다. `rt_solver.c:4541–4548`의 direct-water arrays를 full atmospheric mode count로 zero-pad하고, `rt_options_t` external bottom-source에 `n_mu/m_max` metadata를 추가해 `rt_solve_case_pol_impl`에서 bounds를 강제하라. Pass-3 `5547–5569`의 full-size zero-padding 구현을 기준으로 하며, CLI 최소값 강제나 출력 clamp로 우회하지 말 것. ASan/UBSan, m=(2, 4, 8, 16, 30, 48, 64), Item 6 162-cell 회귀를 완료해야 병합한다.

---

## 10. 근거 파일

- `evidence/nmu48_water_fourier_stability.csv`
- `evidence/nmu48_water_fourier_stability_defect.png`
- `evidence/source_code_evidence.txt`
