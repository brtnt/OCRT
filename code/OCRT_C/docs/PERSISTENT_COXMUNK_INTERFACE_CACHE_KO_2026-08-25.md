# OCRT versioned persistent Cox–Munk/interface operator cache

## 1. 적용 범위

이 문서는 T_wa incoming-U Q-fold와 S1/S2 속도 패치가 반영된 OCRT C 정본에 추가한 마지막 성능개선 항목을 정의한다.

| 연산자 | 저장 단위 | 저장 규약 |
|---|---|---|
| 대기측 거친 해면 반사 `R_aa` | 0…m_max all-mode slab | public wrapper의 U-column fold 전 raw 계수 |
| 해수측 내부반사 `R_ww` | 0…m_max all-mode slab | raw 계수 |
| 공기→해수 투과 `T_aw` | Fourier mode별 완성 행렬 | incoming-U fold 적용 완료 |
| 해수→공기 투과 `T_wa` | 0…m_max all-mode slab | T_wa Q-fold 전 raw 계수 |

output-only exact-view projection, contracted interface operator, atmosphere/ocean subsystem R/T cache, Jacobian 재사용은 이번 정본에 포함하지 않는다.

## 2. 플랫폼 중립 구현 계약

- 소스는 특정 운영체제의 파일·process·locking API를 호출하지 않는다.
- persistent-cache core는 ISO C `fopen/fread/fwrite/fflush/fclose/remove/rename`만 사용한다.
- cache directory는 실행 전 외부에서 생성한다.
- path, process ID, mutex, lock file, filesystem polling에 대한 OS별 분기가 없다.
- 운영체제 이름은 성능·수치 acceptance 조건이 아니다.

기존 조건부 compatibility shim은 동일 source의 이식성을 위한 최소 계층이며, persistent-cache 물리·성능 경로는 OS별로 분기하지 않는다.

## 3. 행 독립성 계약

시뮬레이션 행은 서로의 계산 결과에 의존하지 않는다.

- cache 생성은 simulation batch 이전의 명시적 **단일 prebuilder** 단계이다.
- production batch는 immutable cache 파일을 읽기만 한다.
- lock file, mutex, OpenMP critical, polling, sleep, retry-wait를 사용하지 않는다.
- `--batch-full-grid` 안에서는 `build` mode를 fail-loud로 거부한다.
- `required` mode에서 miss/stale/corrupt이면 즉시 실패하며 다른 row를 기다리지 않는다.
- `p2b_row` hot cache는 `_Thread_local`이다.
- process 시작·완료 순서는 수치 결과에 영향을 주지 않는다.

## 4. 실행 mode

```text
OCRT_SURFACE_PERSIST_CACHE_MODE=off|readonly|required|build
OCRT_SURFACE_PERSIST_CACHE_DIR=<existing directory>
OCRT_SURFACE_PERSIST_CACHE_TRACE=0|1
```

| mode | 동작 | batch 허용 |
|---|---|---|
| `off` | persistent cache 미사용 | 허용 |
| `readonly` | exact hit만 load; miss는 각 row가 독립 계산하고 파일은 쓰지 않음 | 허용 |
| `required` | exact cache만 load; miss/stale/corrupt 즉시 실패 | 권장 |
| `build` | 별도 single-prebuilder run에서 cache 생성 | `--batch-full-grid` 금지 |

## 5. 파일 버전과 cache key

```text
ocrt_surface_v1_<operator>_<dual-64-bit-key-hash>.bin
```

Header/key 구성:

```text
format version
operator ABI version
surface physics version
T_wa Q-fold convention version
operator kind
Fourier mode first/count
n_o, n_i, n_phi
sigma type, Q convention
wind speed, n_water
mu_o/mu_i exact binary64 arrays
payload element count
key dual hash
payload dual hash
```

현재 version:

```text
format   = 1
ABI      = 2026082501
physics  = 2026082501
Q-fold   = 1
```

## 6. 플랫폼 중립 prebuild

```bash
python scripts/prebuild_surface_cache.py \
  --ocrt <OCRT_EXECUTABLE> \
  --cache-dir <CACHE_DIRECTORY> \
  -- <OCRT_ARGUMENTS>
```

도구는 cache directory를 생성하고 다음 환경을 child process에만 설정한다.

```text
OMP_NUM_THREADS=1
OCRT_SURFACE_PERSIST_CACHE_MODE=build
OCRT_SURFACE_PERSIST_CACHE_DIR=<CACHE_DIRECTORY>
```

완료 후 `SURFACE_CACHE_SHA256.csv`를 생성한다.

Production rows:

```text
OMP_NUM_THREADS=1
OCRT_SURFACE_PERSIST_CACHE_MODE=required
OCRT_SURFACE_PERSIST_CACHE_DIR=<CACHE_DIRECTORY>
```

## 7. 성능 비교 원칙

OCRT와 reference를 비교할 때는 어느 지원 OS에서도 다음을 동일하게 맞춘다.

- 동일 host
- 동일 compiler optimization class
- 동일 numerical/physical options
- 동일 layer/quadrature/order/Fourier settings
- 동일 output geometry/workload
- 동일 fresh 또는 warm cache definition
- 반복 median

특정 OS에서만 얻은 결과를 공식 정본으로 지정하지 않는다.

## 8. 검증 Gate

1. cache off/build/required full CSV 회귀
2. Rrs/rrs/TOA I/Q/U 1:1 산포도
3. T_wa spin-2 cache mode parity
4. version/header/payload corruption fail-loud
5. `build` + `--batch-full-grid` rejection
6. serial/OpenMP/process 결과 동일
7. process 실제 overlap
8. `critical(p2b_rowcache)` 부재
9. EAP runtime OFF/fail-loud
10. persistent-cache source의 OS-specific dependency 부재

모든 검증자료는 exact source CSV, metrics, figures, commands, manifest를 원자적으로 함께 보존한다.
