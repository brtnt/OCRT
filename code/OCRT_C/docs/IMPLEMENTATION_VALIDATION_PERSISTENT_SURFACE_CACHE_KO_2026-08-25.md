# OCRT versioned persistent Cox–Munk/interface operator cache 구현·검증 보고서

## 1. 판정

2026-08-25 C 정본에 versioned persistent Cox–Munk/interface operator cache를 통합했다.

```text
구현                         PASS
clean full CMake build         PASS
CTest                         2/2 PASS
cache off/build/required      full CSV byte-identical
Rrs/rrs/TOA I/Q/U scatter     exact 1:1, 624 rows
전체 numeric column regression max |Δ| = 0
T_wa spin-2 regression        PASS
required cache miss           fail-loud PASS
version/corruption rejection  fail-loud PASS
batch build-mode rejection    PASS
4-thread row independence     PASS
4-process row independence    PASS
EAP runtime policy            OFF/fail-loud 유지
OS-specific official gate      폐기; platform-neutral same-host 계약으로 대체
추가 속도 개선                동결
```

이번 구현 이후 output-only exact-view projection, contracted interface operator, subsystem R/T/Jacobian reuse 등 추가 속도 개선은 진행하지 않는다.

## 2. 기준 source

입력 정본:

```text
OCRT_C_SPEEDPATCHED_EAP_DISABLED_2026-08-25.tar.gz
```

기존 정본에서 보존한 사항:

- T_wa incoming-U Q-fold 수정
- FR631 reader compatibility
- air-side all-mode Fourier kernel cache
- 실제 m_max build cap
- air–water interpolation sort/dedup memoization
- EAP production runtime OFF와 명시 selector fail-loud

## 3. 구현 내용

### 3.1 persistent cache 대상

| operator | persistent payload |
|---|---|
| `R_aa` | pre-fold all-mode slab |
| `R_ww` | raw all-mode slab |
| `T_aw` | incoming-U fold가 적용된 final per-mode matrix |
| `T_wa` | T_wa Q-fold 적용 전 raw all-mode slab |

`R_aa`와 `T_wa`는 기존 public wrapper가 fold를 정확히 한 번 적용한다. persistent load가 fold의 중복 적용이나 누락을 만들지 않는다.

### 3.2 mode

```text
OCRT_SURFACE_PERSIST_CACHE_MODE=off|readonly|required|build
OCRT_SURFACE_PERSIST_CACHE_DIR=<directory>
OCRT_SURFACE_PERSIST_CACHE_TRACE=0|1
```

- `build`: batch 이전에 cache를 생성하는 명시적 priming run
- `required`: production row가 immutable cache만 읽음; miss/stale/corrupt 즉시 실패
- `readonly`: optional hit; miss는 해당 row가 local memory에서 계산하고 파일은 쓰지 않음
- `off`: 기존 in-process cache만 사용

### 3.3 versioning 및 무결성

현재 cache version:

```text
format   1
ABI      2026082501
physics  2026082501
Q-fold   1
```

key에는 operator, Fourier mode 범위, quadrature 크기, n_phi, Cox–Munk sigma type, Q convention, wind, n_water, exact `mu_o/mu_i` binary arrays를 포함한다. 파일은 dual key hash와 dual payload hash를 갖는다. 파일을 열 때 header와 모든 방향격자를 bit-for-bit로 확인한다.

### 3.4 게시 방식

`build` mode는 명시적 single-prebuilder 단계이다. Cache core는 ISO C file I/O와 `rename`만 사용하며 OS별 process-ID, directory, lock API를 호출하지 않는다. lock, poll, sleep, retry-wait는 없다. Interrupted/partial file은 header와 payload hash 검증에서 거부된다.


### 3.5 build-system integrity repair

기존 정본의 `test_fr631_mie_reader` target은 `mie_io.c`가 참조하는 `io_utils.c`와 `numerics.c`가 link source에 없어 전체 CMake build가 실패하는 상태였다. 이를 CMakeLists와 Makefile에서 보완하고, FR631 reader CTest 입력을 실제 631×771 canonical `Red_clay_AHN.mie`로 교정했다. `ocrt_nadir_smoke`에는 필요한 `OCRT_ADVANCED=1` test environment를 명시했다. 이는 RT 물리나 속도 알고리즘 변경이 아니라 release build/test 무결성 보완이다.

최종 결과:

```text
cmake --build <build-dir>   PASS
ctest --output-on-failure  2/2 PASS
```

## 4. simulation row 독립성

사용자 요구사항에 따라 행 간 선후행 의존을 명시적으로 금지했다.

### 4.1 build mode의 batch 금지

`--batch-full-grid` 진입 전에 persistent cache mode를 확인한다. `build`이면 계산 시작 전에 exit code 2로 종료한다.

따라서 다음 구조는 허용되지 않는다.

```text
row 1이 cache 생성
→ row 2…N이 row 1 완료를 기다림
```

운영 구조는 다음과 같다.

```text
single-prebuilder 단계
→ immutable cache 완성
→ 독립 process들이 순서와 무관하게 required-mode로 동시 읽기
```

### 4.2 p2b hot cache lock 제거

기존 `rt_air_water_coupling.c`의 `g_p2b_rows`는 process-global mutable cache였고 `critical(p2b_rowcache)`로 보호됐다. 이를 다음과 같이 변경했다.

```text
process-global + OpenMP critical
→ _Thread_local worker-private cache
```

lookup과 insert의 OpenMP critical을 제거했다. 어느 row도 다른 row의 p2b build를 기다리지 않는다.

### 4.3 persistent file 공유의 의미

여러 row/process가 같은 파일을 읽을 수 있으나 파일은 immutable이다. 이는 계산 결과나 mutable state 공유가 아니다. cache miss in `required` mode는 다른 process를 기다리지 않고 해당 실행을 즉시 실패시킨다.

### 4.4 동적 검증

4개 row를 다음 세 방식으로 실행했다.

```text
serial, OMP_NUM_THREADS=1
one process / 4 OpenMP workers
four independent processes, each OMP_NUM_THREADS=1
```

결과:

- 4개 row 전체 CSV가 세 방식에서 byte-identical
- 4-process latest start: 약 0.0019 s
- 4-process first finish: 약 6.3305 s
- 모든 process가 첫 process 종료 전에 시작하여 실제 overlap 확인
- `critical(p2b_rowcache)` source 잔존 0건

## 5. 수치 무회귀 검증

### 5.1 조건

```text
wavelength     555 nm
SZA            40°
wind           3 m/s
n_water        1.34
n_mu_water     48
water tau cap  16.6358629872763
VZA            0–60°, 5°
RAA            0–345°, 15°
geometry       312/case
cases          pure water, Red clay TSM 0.5 g m-3
PSSA           OFF
absorbing gas  OFF
EAP            OFF
```

### 5.2 full-file equality

각 case에서 다음 세 파일이 byte-identical이다.

```text
persistent cache OFF
persistent cache BUILD
persistent cache REQUIRED
```

Pure-water SHA-256:

```text
6bbf5356f7132b83ba85ca8054d1beb608e87da1b8b605bb3d8df2b72401cb9f
```

Red-clay SHA-256:

```text
9e9390209e2410eb8ebcf58a86531e1667c84ffb80192f4099177e0e75fb540c
```

### 5.3 산포도

두 case 624개 geometry를 결합해 다음 9개 산포도를 생성했다.

```text
Rrs I/Q/U
rrs I/Q/U
TOA rho I/Q/U
```

모든 항목:

```text
max absolute difference = 0
RMSE                    = 0
slope                   = 1
R²                      = 1
verdict                  = PASS
```

CSV의 59개 numeric column 전체를 검사한 236개 case/comparison 조합도 max absolute difference와 RMSE가 모두 0이다.

## 6. T_wa Q-fold 및 cache convention 검증

합성 spin-2 test를 cache off/build/required로 실행했다. 세 log는 byte-identical하다.

VZA 2–40°에서:

```text
T_I  = 0.999–1.005
T_Q2 = 0.999–1.005
T_U2 = 0.999–1.005
```

persistent T_wa payload는 pre-fold로 저장되고 load 후 public wrapper가 Q-fold를 한 번 적용한다. T_aw payload는 기존 direct 함수의 final folded 결과를 저장한다.

## 7. fail-loud Gate

| 시험 | exit code | 판정 |
|---|---:|---|
| `required` + empty directory | 3 | PASS |
| header physics-version byte 손상 | 3 | PASS |
| `build` + `--batch-full-grid` | 2 | PASS |
| explicit EAP selector | 2 | PASS |

손상 cache는 `stale-or-corrupt`로 거부됐으며 fallback 결과를 승인하지 않는다.

## 8. cache inventory

검증 profile에서 생성된 예시 cache:

| operator | files | total |
|---|---:|---:|
| R_aa | 4 | 0.99 MiB |
| R_ww | 2 | 7.39 MiB |
| T_aw | 6 | 0.96 MiB |
| T_wa | 74 | 39.71 MiB |
| 합계 | 86 | 약 49.1 MiB |

이는 validation example이다. 실제 campaign cache 개수는 unique wind, n_water, quadrature, Fourier limit 및 view grid 조합에 따라 달라진다.

## 9. Same-host 진단 timing

동일 container에서 cache I/O 효과를 확인한 진단값이다. 특정 OS를 공식 Gate로 지정하지 않으며, 모든 성능 비교는 동일 host/build/options/cache-state 조건에서 해석한다.

| case | off | build | required | off/required |
|---|---:|---:|---:|---:|
| pure water | 3.10 s | 3.12 s | 0.66 s | 4.70× |
| Red clay 0.5 | 6.34 s | 5.39 s | 5.30 s | 1.20× |

Red-clay에서는 water RT 자체가 지배해 persistent surface cache 이득이 제한적이다. 이번 사용자의 지시에 따라 이 결과 이후 추가 속도 구조 변경은 하지 않는다.

## 10. Windows 상태

소스는 `_WIN32` 경로, `_getpid`, Windows path separator를 지원하며 PowerShell prebuild script를 포함한다. 그러나 현재 세션 실행환경은 Linux이므로 native Windows build/wall-time Gate를 실제 통과했다고 주장하지 않는다.

Windows에서 수행할 작업은 성능 패치가 아니라 다음 release 확인이다.

1. 최신 source native build
2. cache prebuild
3. `required` mode로 공정 fresh/warm 반복
4. I/Q/U 산포도 및 exact source CSV 재생성

## 11. 속도 개선 동결

다음 항목은 별도 사용자 지시 전까지 보류한다.

- output-only exact-view projection
- contracted T_aw/T_wa cache 재설계
- 5-block subsystem R/T operator
- Jacobian partial rebuild
- 기타 solver 구조 최적화

현재 release의 병렬 실행 계약은 `16 independent processes × OMP_NUM_THREADS=1`, immutable persistent cache read-only, row 간 wait/lock 없음이다.


## 10. 2026-08-26 플랫폼 중립 롤백 검증

- native-Windows 전용 gate·runner·PowerShell helper를 최종 정본에서 제외했다.
- `surface_persistent_cache.c`에서 `_WIN32`, POSIX `unistd`, process-ID, `stat` 기반 directory probe를 제거했다.
- persistent-cache core는 ISO C file I/O만 사용한다.
- generic helper `scripts/prebuild_surface_cache.py`로 prebuild 절차를 통일했다.
- cache format/key/payload version은 유지하여 기존 valid cache와 수치 계약을 바꾸지 않았다.
- clean build, CTest, cache off/build/required regression, row-independence를 재실행해 PASS를 확인했다.
- native-Windows 전용 결과를 release blocker로 삼는 규칙은 폐기했다.
