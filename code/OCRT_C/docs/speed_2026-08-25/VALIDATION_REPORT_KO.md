# OCRT S1/S2 속도 패치 적용·무회귀 검증 보고서

- 작성일: 2026-08-25
- 기준 트리: `MIGRATION_PKG_2026-08-19/01_OCRT_C`
- 선행 통합: EAP runtime-disabled overlay + T_wa incoming-U Q-fold fix + FR631 771-column reader compatibility
- 대상 패치: air-side Cox–Munk all-m Fourier kernel cache, 실제 `m_max` build cap, air–water interpolation sort-map memoization
- EAP 상태: production runtime OFF, 명시 selector fail-loud

## 1. 통합 판정

**PASS.** 첨부 속도 패치를 기존 EAP-disabled/T_wa-Q-fold 정본 위에 hunk 단위로 적용했다. 첨부된 완성본 파일을 무조건 덮어쓰지 않았으며, 다음 충돌을 보존·해소했다.

1. `rt_solver.c`: 기존 Windows compatibility include를 유지했다.
2. `rt_water_rt.c`: 기존 Windows compatibility include와 EAP scattering fail-loud guard를 유지했다.
3. `surface.c`: T_wa FKC hit와 direct path 모두에서 incoming-U column `IU/QU` fold가 한 번만 적용되는 구조를 유지했다.
4. FR631 line buffer와 AHN 771-column capacity 패치를 유지했다.
5. multi-geometry one-run과 `water_cold=1`, `replay=0` 경로를 변경하지 않았다.

첨부 완성본과 통합본 비교 결과:

- `surface.c`, `surface.h`, `rt_air_water_coupling.c`: SHA-256 완전 일치
- `rt_solver.c`: 통합본에 `rt_windows_compat.h` include 한 줄 추가 유지
- `rt_water_rt.c`: 통합본에 `rt_windows_compat.h`와 EAP fail-loud guard 유지

## 2. 적용 내용

### S1a — air-side rough-surface all-m FKC

`surface_coxmunk_fourier_kernel()`가 Fourier mode마다 동일한 per-azimuth Cox–Munk Mueller kernel을 다시 적분하지 않고, 한 번의 azimuth pass에서 필요한 모든 mode를 투영한다.

- thread-local 8-slot LRU
- exact scalar/grid byte comparison key
- `OCRT_AIR_FKC_OFF=1` authoritative rollback path
- m>0 `IU/QU` sign fold 유지

### S1b — 실제 Fourier 상한 전달

`surface_fkc_set_build_mmax()`를 다음 세 m-loop 앞에 배치했다.

- atmosphere air-side surface kernel
- water-side `R_ww`
- rough `T_wa`

compile ceiling 64 전체를 무조건 만들지 않고 실제 loop 상한까지만 생성한다.

### S2 — air–water interpolation sort-map memoization

`aw_interp_on_unsorted()`가 동일한 direction grid를 매 component/mode/target마다 insertion sort와 dedup하지 않도록 thread-local map을 재사용한다.

- input grid는 byte comparison으로 확인
- 보간에 들어가는 값과 연산 순서는 기존 경로와 동일
- 320개 초과 grid는 기존 heap fallback 유지

## 3. Build

GNU C 14.2.0, CMake/Ninja release build를 사용했다.

- `-O3`
- `-ffp-contract=fast`
- `-fassociative-math`
- `-DOCRT_FAST_KERNELS`
- OpenMP link
- build warning: 0

속도 패치 Linux binary SHA-256:

```text
b12d1423c719005614b6e14f2f5d4a54878d8e0641d113657cea44ae96d76109
```

## 4. 무회귀 결과

### 4.1 Pure water 및 Red clay full-grid

조건:

```text
555 nm, SZA 40°, wind 3 m/s, n_water 1.34
n_mu_water 48, water tau cap 16.6358629872763
VZA 0–60°/5°, RAA 0–345°/15°
312 geometry/case
PSSA OFF, gases OFF, EAP OFF
```

| Case | 비교 | 결과 |
|---|---|---|
| Pure water | baseline vs speed-patched | full CSV byte-identical |
| Red clay 0.5 g m⁻³ | baseline vs speed-patched | full CSV byte-identical |
| Red clay | baseline vs speed-patched + `OCRT_AIR_FKC_OFF=1` | full CSV byte-identical |

Rrs I/Q/U 624개 paired sample의 최대 절대차는 전부 0이다. 산포도와 exact source CSV를 함께 보존했다.

### 4.2 Wind 0

speed-patched FKC ON/OFF 결과가 byte-identical이다.

```text
SHA-256 = 89c5a6aa4edf7136c21293d40d6740afaa0c932d2a98afb6c1477515812d8b5c
```

### 4.3 PSSA smoke

조건: 555 nm, SZA 75°, wind 5 m/s, black Fresnel ocean, 1,296 geometry.

baseline과 speed-patched CSV가 byte-identical이다.

```text
SHA-256 = 552524adeab28f535d90b18fa22cc4991b2787332d826f84950778019c15a843
```

### 4.4 T_wa spin-2 gate

VZA 2–40°에서:

```text
T_I, T_Q2, T_U2 = 0.999–1.005
```

FKC ON/OFF text output도 byte-identical이다.

### 4.5 Q-fold 작업지시서 재검증

새 speed-patched binary로 작업지시서 검증기를 재실행했다.

| Case | I median/max | Q median/max | U median/max | 판정 |
|---|---|---|---|---|
| Pure water | +0.240 / 0.298 | −0.029 / 0.127 | ≈0 / 0.091 | PASS |
| Red clay | +0.646 / 0.889 | −0.012 / 0.085 | ≈0 / 0.091 | PASS |

6/6 PASS다.

### 4.6 FR631 reader

Canonical Red clay 파일:

```text
n_ang=631
n_wl=771
n_phase_wl=771
status=PASS
```

### 4.7 EAP fail-loud

명시적 `--ocrt-phyto-group micro` 요청은 exit code 2로 종료했다. 일반 Chl 경로로 fallback하거나 EAP Mie를 load하지 않았다.

## 5. 속도 결과

### 5.1 정례/회귀 profile

단일 반복, 같은 container, `OMP_NUM_THREADS=1`:

| Profile | Baseline | Speed-patched | 가속 |
|---|---:|---:|---:|
| Pure water, 312 geometry | 11.10 s | 5.66 s | 1.96× |
| Red clay, 312 geometry | 17.44 s | 11.06 s | 1.58× |
| PSSA BFO, 1,296 geometry | 0.53 s | 0.27 s | 1.96× |

### 5.2 공정 B2a profile

조건:

```text
555 nm, SZA 50°, wind 5 m/s, Red clay 0.5 g m⁻³
PSSA OFF, gases OFF, EAP OFF
atmosphere/water positive Gauss nodes 48/48
atmosphere/water layers 26/640
order caps 100/100
Fourier m_max 32
VZA 18 × RAA 72 = 1,296 geometry
OMP_NUM_THREADS=1
```

Speed-patched OCRT 3회:

```text
27.21 s, 26.77 s, 24.60 s
median = 26.77 s
```

동일 profile의 unpatched baseline은 420 s timeout 내 완료되지 않았다. 따라서 이 환경에서 속도향상은 **15.69× 초과의 하한**만 제시하며, 정확한 baseline 배율로 간주하지 않는다. 첨부 문서의 3.25×와 차이가 큰 것은 CPU/build/runtime data 상태가 다르고 baseline이 timeout됐기 때문이다.

### 5.3 OSOAA 비교의 정확한 범위

같은 host, 같은 red-clay fair profile에서:

| OSOAA 상태 | OSOAA | OCRT median | OCRT/OSOAA | 판정 |
|---|---:|---:|---:|---|
| `SURF_MATR` fresh 생성 포함 | 33.47 s | 26.77 s | 0.800 | OCRT 1.25× 빠름 |
| 기존 `SURF_MATR` 재사용 | 20.99 s | 26.77 s | 1.275 | OCRT가 1.28× 느림 |

따라서 **“OSOAA보다 빨라졌다”는 주장은 fresh surface-matrix end-to-end 조건에서는 확인됐지만, OSOAA warm cache 재사용 조건까지 일반화하면 아직 성립하지 않는다.** 공식 판정은 native Windows에서 fresh와 warm을 분리해 반복 median으로 다시 수행해야 한다.

두 번째 OSOAA fresh 반복은 surface `TWA` matrix 재생성 도중 환경 timeout으로 중단되어, OSOAA 수치는 현재 각 상태 1회 진단값이다. 이를 공식 median으로 표현하지 않는다.

## 6. 구현 순서 갱신

1. **완료 — 이번 S1/S2 patch**
   - air-side all-m FKC
   - actual m_max cap
   - interpolation sort-map memo
   - Q-fold/EAP/FR631 무회귀
2. **다음 — native Windows fair speed gate + warm OSOAA gap 해결**
   - OCRT persistent/versioned surface-operator cache 검토
   - T_aw 전용 all-m builder
   - T_wa view-node batch key
3. output-only exact near-nadir projection으로 coupled-call 중복 제거
4. contracted `T_aw/T_wa` operator cache와 direct indexing
5. 5-subsystem R/T operator cache 및 Jacobian partial rebuild

층분리 R/T/Jacobian 구조는 유지하되, 이번 patch와 warm-cache gap 해결 뒤에 구현한다.

## 7. 승인 판정

| Gate | 판정 |
|---|---|
| speed patch merge | PASS |
| EAP-disabled 정책 보존 | PASS |
| T_wa Q-fold 보존 | PASS |
| FR631 631×771 reader | PASS |
| build warning 0 | PASS |
| baseline–optimized I/Q/U 산포도 | PASS, exact |
| fresh OSOAA speed diagnostic | PASS |
| warm OSOAA speed diagnostic | FAIL — 후속 최적화 필요 |
| native Windows official gate | PENDING |

이번 patch는 production C tree에 채택한다. 다만 장기 편광 campaign의 최종 speed release gate는 native Windows warm/fresh 결과와 I/Q/U OSOAA 산포도까지 봉인한 뒤 승인한다.
