# OCRT C 최종 플랫폼 중립 성능개선 정본 — 2026-08-26

이 tree가 현재 C 정본이다.

## 포함

- T_wa Q-fold fix
- FR631 reader compatibility
- EAP runtime OFF/fail-loud
- all-mode Cox–Munk FKC
- true `m_max` cap
- interpolation map memoization
- versioned persistent `R_aa/R_ww/T_aw/T_wa` cache
- worker-local row cache 및 row 간 lock 제거

## 플랫폼 정책

- 특정 OS는 필수 실행환경이나 공식 성능 Gate가 아니다.
- 소스의 물리·성능 경로를 OS별로 분기하지 않는다.
- persistent-cache core는 ISO C file I/O만 사용한다.
- native-Windows gate와 Windows-only helper는 최종 정본에서 제외했다.
- 성능 비교는 동일 host/build/options/cache-state 조건으로 수행한다.

## Cache prebuild

```bash
python scripts/prebuild_surface_cache.py \
  --ocrt <OCRT_EXECUTABLE> \
  --cache-dir <CACHE_DIRECTORY> \
  -- <OCRT_ARGUMENTS>
```

Simulation rows:

```text
OCRT_SURFACE_PERSIST_CACHE_MODE=required
OCRT_SURFACE_PERSIST_CACHE_DIR=<CACHE_DIRECTORY>
OMP_NUM_THREADS=1
```

각 row는 다른 row의 계산이나 cache 생성을 기다리지 않는다.

## 동결

위 성능개선까지만 최종 적용한다. 추가 solver 최적화는 별도 승인 전까지 진행하지 않는다.

상세 결정: `docs/NATIVE_WINDOWS_ROLLBACK_PLATFORM_NEUTRAL_FINAL_2026-08-26_KO.md`
