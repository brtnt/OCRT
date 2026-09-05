# OCRT native-Windows 롤백 및 플랫폼 중립 최종 정본 검증 보고서

**일자:** 2026-08-26  
**최종 source:** `OCRT_C_FINAL_PLATFORM_NEUTRAL_PERF_2026-08-26`

## 1. 최종 판정

```text
native-Windows 전용 Gate/runner/workflow     ROLLED BACK
특정 OS 필수 실행·승인 규칙                  REMOVED
persistent-cache core 신규 OS dependency     0
PowerShell script in final source             0
Windows executable in final source            0
clean CMake build                             PASS
CTest                                          2/2 PASS
previous canonical vs final full CSV          BYTE-IDENTICAL
Rrs I/Q/U scatter                              exact 1:1, 624 rows
cache format old vs final                      5/5 BYTE-IDENTICAL
cache off/build/required                       BYTE-IDENTICAL
row independence                               PASS
EAP runtime                                    OFF / fail-loud
추가 속도개선                                  FROZEN
```

최종 정본은 특정 운영체제의 실행·파일·process API를 성능 경로에 추가하지 않는다. 기존 조건부 compatibility shim은 이전 정본과 byte-identical하게 유지했으며, 이는 OS별 물리/성능 분기가 아니라 동일 소스 이식성을 위한 기존 최소 계층이다.

## 2. 롤백 항목

다음 산출물과 정책은 최종 정본에서 사용하지 않는다.

- native-Windows fresh/warm speed gate
- Windows runtime 조립·재빌드 runner
- Windows-only PowerShell prebuild helper
- “Windows 결과만 공식 timing”이라는 승인 규칙
- 특정 OS 이름을 campaign 시작 조건으로 사용하는 규칙

`metrics/DEPRECATED_NATIVE_WINDOWS_ARTIFACTS.csv`에 비정본 산출물을 기록했다.

## 3. 유지한 최종 성능 패치

1. T_wa incoming-U Q-fold
2. air-side Cox–Munk all-mode FKC
3. true Fourier `m_max` build cap
4. air–water interpolation sort/dedup map memoization
5. versioned persistent `R_aa/R_ww/T_aw/T_wa` cache
6. `p2b_row` worker-local cache 및 row 간 OpenMP critical 제거
7. FR631 reader
8. EAP production OFF/fail-loud

추가 solver 최적화는 중단했다.

## 4. Persistent cache 플랫폼 중립화

수정 파일:

```text
src/shared/surface_persistent_cache.c
src/shared/surface_persistent_cache.h
```

제거한 dependency:

```text
_WIN32 branch
direct.h / process.h
unistd.h / sys/stat.h
_getpid / getpid
OS-specific path separator macro
Windows-only PowerShell helper
```

현재 cache core가 사용하는 외부 기능은 ISO C file I/O와 standard C types뿐이다. Cache directory 생성은 solver 외부에서 처리한다. Generic helper는 다음이다.

```text
scripts/prebuild_surface_cache.py
```

Cache format/key/payload version은 바꾸지 않았다.

## 5. Build 및 test

```text
CMake Release configure   PASS
CMake build               PASS
CTest                     2/2 PASS
FR631 reader              PASS
OCRT nadir smoke          PASS
```

로그:

```text
logs/platform_neutral_cmake_configure.log
logs/platform_neutral_cmake_build.log
logs/platform_neutral_ctest.log
```

## 6. 수치 무회귀

조건:

```text
555 nm
SZA 40°
wind 3 m/s
n_water 1.34
n_mu_water 48
VZA 0–60°, 5°
RAA 0–345°, 15°
312 geometry/case
pure water + Red clay TSM 0.5 g m-3
PSSA OFF
absorbing gases OFF
EAP OFF
```

비교:

```text
previous persistent-cache canonical binary
vs
platform-neutral final binary
```

두 case의 full CSV는 byte-identical이다. 624개 geometry의 Rrs I/Q/U 결과:

| 성분 | n | max abs diff | RMSE | slope | intercept | R² |
|---|---:|---:|---:|---:|---:|---:|
| I | 624 | 0 | 0 | 1 | 0 | 1 |
| Q | 624 | 0 | 0 | 1 | 0 | 1 |
| U | 624 | 0 | 0 | 1 | 0 | 1 |

산포도와 exact source CSV는 `figures/`, `source/`, `metrics/`에 함께 보존했다.

## 7. Cache format 회귀

동일 profile에서 이전 구현과 최종 구현이 생성한 5개 persistent-cache binary를 비교했다.

```text
filename              identical
byte count            identical
SHA-256               identical
result                 5/5 PASS
```

즉, 플랫폼 중립화는 cache ABI, key, payload 또는 수치 결과를 변경하지 않았다.

## 8. 행 독립성

재실행 결과:

- serial vs OpenMP 4 workers: byte-identical
- serial vs 4 independent processes: byte-identical
- OpenMP rows 실제 overlap: PASS
- independent processes 실제 overlap: PASS
- `critical(p2b_rowcache)` 잔존: 0
- batch 내부 build mode: fail-loud
- required cache miss: fail-loud

행들은 다른 행의 계산 완료나 cache 생성을 기다리지 않는다.

## 9. 플랫폼 중립 성능 비교 계약

향후 OCRT와 reference 속도를 다시 비교할 때 OS를 고정하지 않는다. 다음만 동일하게 맞춘다.

```text
same host
same compiler optimization class
same physics and numerical resolution
same process/thread count
same output workload
same fresh/warm cache definition
repeated median
I/Q/U regression artifacts
```

이번 정본에서는 추가 속도개선과 OS-specific gate를 모두 중단한다.
