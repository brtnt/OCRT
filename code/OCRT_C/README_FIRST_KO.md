# OCRT 속도향상 최종 코드 패키지 — Mie payload 제외

**패키지:** `OCRT_SPEED_OPTIMIZED_PLATFORM_NEUTRAL_NO_MIE_2026-08-26`  
**기준 정본:** `OCRT_C_FINAL_PLATFORM_NEUTRAL_PERF_2026-08-26`  
**기준 archive SHA-256:** `0debb9b92f652ee6a9621b28ff1ce40cbeb0d566b8f0080b3c9c1c8580cab629`

## 1. 패키지 범위

현재 승인된 C 성능개선 정본을 그대로 포함한다.

- T_wa incoming-U Q-fold 수정
- air-side Cox–Munk all-mode Fourier kernel cache
- 실제 Fourier `m_max` build cap
- air–water interpolation sort/dedup map memoization
- versioned persistent `R_aa/R_ww/T_aw/T_wa` interface-operator cache
- worker-local `p2b_row` cache 및 row 간 lock 제거
- persistent-cache core의 플랫폼 중립 ISO C file I/O
- FR631 631-angle × 771-wavelength reader compatibility
- EAP production runtime OFF 및 명시 selector fail-loud

추가 solver 최적화는 동결 상태다.

## 2. Mie 자료 제외

요청에 따라 **패키지 내부 `.mie` 파일 수는 0개**다.

- 현재 production에서 사용할 수 있는 canonical non-EAP 자료: 외부 181개
- EAP canonical 자료: 17개, archival only; production에 설치·로드하지 않음
- 최신 canonical 구조: 330–1100 nm, 1 nm, 771 wavelength columns, FR631 631 angles, P11/P12/P33

정확한 경로·크기·SHA-256:

```text
_RELEASE/external_mie/EXTERNAL_MIE_RUNTIME_REQUIRED_181.csv
_RELEASE/external_mie/EAP_MIE_ARCHIVAL_ONLY_17.csv
_RELEASE/external_mie/MIE_CANONICAL_INVENTORY_198.csv
_RELEASE/external_mie/MIE_198_FULL_AUDIT.csv
```

외부 자료 설치 확인:

```bash
python _RELEASE/tools/verify_external_mie.py --ocrt-root .
```

검증기는 181개 non-EAP 파일의 size/SHA-256을 확인하고, production EAP `.mie`가 설치돼 있으면 실패한다.

## 3. Build

Mie 파일은 compile 단계에 필요하지 않다.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

본 패키지를 Mie 없이 별도 build한 결과 `cmake --build`가 PASS했고, `ocrt_nadir_smoke` 1/1이 PASS했다.
`fr631_mie_reader`와 aerosol/TSM/detritus RT test는 외부 Mie 설치 후 수행한다. Build log는 `_RELEASE/build_audit/`에 있다.

## 4. Persistent surface cache와 행 독립성

Cache는 simulation batch 전에 별도 prebuild한다.

```bash
python scripts/prebuild_surface_cache.py   --ocrt <OCRT_EXECUTABLE>   --cache-dir <CACHE_DIRECTORY>   -- <OCRT_ARGUMENTS>
```

각 simulation row:

```text
OCRT_SURFACE_PERSIST_CACHE_MODE=required
OCRT_SURFACE_PERSIST_CACHE_DIR=<CACHE_DIRECTORY>
OMP_NUM_THREADS=1
```

각 row는 다른 row의 계산 또는 cache 생성을 기다리지 않는다. `required` mode miss는 즉시 실패한다.

## 5. 검증 상태

기준 정본에서 다음을 확인했다.

- clean release build PASS
- CTest 2/2 PASS — Mie 포함 기준 정본
- cache off/build/required full CSV byte-identical
- pure-water/Red-clay 624 geometry Rrs I/Q/U 산포도 exact 1:1
- serial/OpenMP/process row 결과 byte-identical 및 실제 overlap 확인
- EAP selector fail-loud

이 code-only 패키지의 source/build/test/script 파일은 기준 정본과 byte-identical하며,
`_RELEASE/SOURCE_CODE_EQUIVALENCE_SHA256.csv`에서 확인할 수 있다.

## 6. 패키지 자체 검증

```bash
./verify_package.sh
```

또는:

```bash
python _RELEASE/tools/verify_package.py
```
