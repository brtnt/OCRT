# OCRT C final platform-neutral performance tree — 2026-08-26

This tree is the current C reference tree.

## Included

- T_wa Q-fold fix
- FR631 reader compatibility
- EAP runtime OFF / fail-loud
- all-mode Cox–Munk Fourier kernel cache (FKC)
- true `m_max` cap
- interpolation map memoization
- versioned persistent `R_aa/R_ww/T_aw/T_wa` cache
- worker-local row cache without locks between rows

## Platform policy

- No particular operating system is a required environment or an official performance gate.
- The physics and performance paths of the source do not branch by operating system.
- The persistent-cache core uses ISO C file I/O only.
- Native-Windows gates and Windows-only helpers were removed from the final tree.
- Performance comparisons are made on the same host with the same build, options and cache state.

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

No row waits for another row's computation or cache build.

## Freeze

Only the performance changes listed above are applied. Further solver optimization waits for a
separate approval.

Decision record: `docs/NATIVE_WINDOWS_ROLLBACK_PLATFORM_NEUTRAL_FINAL_2026-08-26_KO.md`

---

# OCRT C 최종 플랫폼 중립 성능개선 정본 — 2026-08-26

이 트리가 현재 C 정본이다.

## 포함

- T_wa Q-fold 수정
- FR631 reader 호환
- EAP 런타임 OFF / fail-loud
- 전 모드 Cox–Munk 푸리에 커널 캐시(FKC)
- 실제 `m_max` 상한
- 보간 맵 메모이제이션
- 버전 관리되는 영속 `R_aa/R_ww/T_aw/T_wa` 캐시
- worker 국소 행 캐시, 행 사이 lock 제거

## 플랫폼 정책

- 특정 운영체제는 필수 실행환경도 공식 성능 게이트도 아니다.
- 소스의 물리·성능 경로를 운영체제별로 분기하지 않는다.
- 영속 캐시 코어는 ISO C 파일 입출력만 쓴다.
- native-Windows 게이트와 Windows 전용 보조 도구는 최종 정본에서 제외했다.
- 성능 비교는 같은 호스트에서 같은 빌드·옵션·캐시 상태로 수행한다.

## 캐시 미리 만들기

```bash
python scripts/prebuild_surface_cache.py \
  --ocrt <OCRT_EXECUTABLE> \
  --cache-dir <CACHE_DIRECTORY> \
  -- <OCRT_ARGUMENTS>
```

시뮬레이션 행:

```text
OCRT_SURFACE_PERSIST_CACHE_MODE=required
OCRT_SURFACE_PERSIST_CACHE_DIR=<CACHE_DIRECTORY>
OMP_NUM_THREADS=1
```

각 행은 다른 행의 계산이나 캐시 생성을 기다리지 않는다.

## 동결

위 성능개선까지만 적용한다. 추가 솔버 최적화는 별도 승인 전까지 진행하지 않는다.

결정 기록: `docs/NATIVE_WINDOWS_ROLLBACK_PLATFORM_NEUTRAL_FINAL_2026-08-26_KO.md`
