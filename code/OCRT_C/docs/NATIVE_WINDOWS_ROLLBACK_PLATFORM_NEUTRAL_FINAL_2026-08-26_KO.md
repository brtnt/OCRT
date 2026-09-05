# OCRT native-Windows 작업 롤백 및 플랫폼 중립 최종화

**일자:** 2026-08-26  
**최종 정본:** `OCRT_C_FINAL_PLATFORM_NEUTRAL_PERF_2026-08-26`

## 1. 결정

특정 운영체제를 공식 실행환경이나 성능 승인환경으로 강제하지 않는다. 이전 단계에서 준비한 native-Windows 전용 speed-gate, PowerShell 실행기, Windows runtime 조립 절차 및 “Windows에서만 공식 판정” 규칙은 모두 폐기한다.

소스의 물리·성능 경로는 운영체제에 따라 달라지지 않아야 한다. 성능 비교는 어느 지원 플랫폼에서도 가능하지만 다음 조건을 동일하게 맞춰야 한다.

- 동일 호스트와 동일 CPU 상태
- 동일 compiler optimization class
- 동일 입력·물리 옵션
- 동일 quadrature, layer, scattering-order, Fourier 설정
- 동일 thread/process 설정
- 동일 fresh/warm cache 상태
- 동일 출력 workload

운영체제 이름 자체는 acceptance 조건이 아니다.

## 2. 롤백 범위

최종 정본에서 제외한 항목:

- `OCRT_NATIVE_WINDOWS_FAIR_WARM_GATE_V2_2026-08-25`
- `OCRT_STAGE2_WINDOWS_NATIVE_GATE_KIT_2026-08-23`
- Windows-only runtime setup/build/gate 절차
- `PREBUILD_SURFACE_CACHE_WINDOWS.ps1`
- native-Windows 결과가 있어야만 campaign을 실행할 수 있다는 문서 규칙
- Windows 전용 wall-time을 공식 정본으로 지정하는 판정

이전 산출물은 역사적 감사자료일 뿐이며 최종 정본의 실행·승인 절차로 사용하지 않는다.

## 3. 유지한 성능·정합 패치

다음 항목까지를 최종 성능개선 범위로 고정한다.

1. T_wa incoming-U Q-fold 수정
2. air-side Cox–Munk all-mode Fourier kernel cache
3. 실제 Fourier `m_max` build cap
4. air–water interpolation sort/dedup map memoization
5. versioned persistent Cox–Munk/interface operator cache
6. `p2b_row` cache worker-local화 및 계산 경로의 row 간 lock 제거
7. EAP production runtime OFF와 명시 selector fail-loud
8. FR631 631-angle × 771-wavelength reader compatibility

추가 구조 최적화는 별도 승인 전까지 동결한다.

## 4. 플랫폼 중립 persistent-cache 구현

`src/shared/surface_persistent_cache.c`는 다음 원칙으로 정리했다.

- ISO C file I/O만 사용
- `_WIN32`, POSIX `unistd`, process-ID API, OS lock API 사용 없음
- cache directory 생성은 solver 밖에서 수행
- build는 명시적 단일 prebuilder 단계
- simulation row는 `required` 또는 `readonly`로 immutable 파일만 읽음
- row 간 lock, wait, poll, retry 없음
- cache miss in `required` mode는 즉시 실패
- cache file format, key, payload hash는 기존 version 1 계약 유지

보조 prebuild 도구는 `scripts/prebuild_surface_cache.py`로 통일했다. 이 도구는 특정 shell이나 운영체제를 요구하지 않는다.

## 5. 기존 호환 shim의 취급

기존 source tree에 있던 조건부 C-runtime compatibility shim은 제거하지 않았다. 이는 특정 OS 최적화 경로가 아니라 동일한 소스를 여러 compiler/runtime에서 빌드하기 위한 최소 이식성 계층이다. 이번 변경에서 해당 shim의 기능을 확장하지 않았고, 물리·성능 알고리즘을 OS별로 분기하지 않았다.

## 6. 병렬 실행 계약

- 각 simulation row는 완전히 독립적이어야 한다.
- 장기 campaign은 독립 process fan-out을 사용한다.
- 각 process는 `OMP_NUM_THREADS=1`을 기본으로 한다.
- cache는 batch 전에 단일 prebuilder로 생성한다.
- production row가 cache를 생성하거나 다른 row의 완료를 기다리는 구조는 금지한다.
- 결과 병합은 모든 독립 계산이 끝난 뒤 수행한다.

## 7. 최종 검증 요구

- clean CMake/Make build
- cache off/build/required 결과 회귀
- Rrs/rrs/TOA I/Q/U 산포도
- exact source CSV, metrics, manifest 동봉
- T_wa spin-2 gate
- cache corruption/miss fail-loud
- serial/OpenMP/process row independence
- EAP OFF/fail-loud
- source의 신규 OS-specific dependency 0건
