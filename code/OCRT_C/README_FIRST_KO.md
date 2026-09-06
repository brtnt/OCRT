# OCRT speed-optimized final code package — Mie payload excluded

**Package:** `OCRT_SPEED_OPTIMIZED_PLATFORM_NEUTRAL_NO_MIE_2026-08-26`
**Reference tree:** `OCRT_C_FINAL_PLATFORM_NEUTRAL_PERF_2026-08-26`
**Reference archive SHA-256:** `0debb9b92f652ee6a9621b28ff1ce40cbeb0d566b8f0080b3c9c1c8580cab629`

## 1. Package scope

The package contains the approved C performance tree as it is.

- T_wa incoming-U Q-fold fix
- air-side Cox–Munk all-mode Fourier kernel cache
- true Fourier `m_max` build cap
- memoized air–water interpolation sort/dedup map
- versioned persistent `R_aa/R_ww/T_aw/T_wa` interface-operator cache
- worker-local `p2b_row` cache without locks between rows
- platform-neutral ISO C file I/O in the persistent-cache core
- FR631 reader compatibility (631 angles × 771 wavelengths)
- EAP production runtime OFF; an explicit selector fails loudly

Further solver optimization is frozen.

## 2. Mie data excluded

By request, the package contains **no `.mie` file**.

- Canonical non-EAP data usable in production: 181 external files
- EAP canonical data: 17 files, archival only; not installed or loaded in production
- Latest canonical layout: 330–1100 nm at 1 nm (771 wavelength columns), FR631 631 angles,
  P11/P12/P33

Exact paths, sizes and SHA-256:

```text
_RELEASE/external_mie/EXTERNAL_MIE_RUNTIME_REQUIRED_181.csv
_RELEASE/external_mie/EAP_MIE_ARCHIVAL_ONLY_17.csv
_RELEASE/external_mie/MIE_CANONICAL_INVENTORY_198.csv
_RELEASE/external_mie/MIE_198_FULL_AUDIT.csv
```

Check of the external data installation:

```bash
python _RELEASE/tools/verify_external_mie.py --ocrt-root .
```

The checker verifies size and SHA-256 of the 181 non-EAP files and fails if a production EAP
`.mie` file is installed. (In the public repository the 181 tables are distributed as GitHub
Release `data-v1`; see the root `README.md`.)

## 3. Build

Mie files are not needed at compile time.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

A separate build of this package without Mie data passed `cmake --build`, and `ocrt_nadir_smoke`
passed 1/1. The `fr631_mie_reader` test and the aerosol/mineral/detritus RT tests run after the
external Mie data are installed. Build logs are in `_RELEASE/build_audit/`.

## 4. Persistent surface cache and row independence

The cache is built separately before a simulation batch.

```bash
python scripts/prebuild_surface_cache.py   --ocrt <OCRT_EXECUTABLE>   --cache-dir <CACHE_DIRECTORY>   -- <OCRT_ARGUMENTS>
```

Each simulation row:

```text
OCRT_SURFACE_PERSIST_CACHE_MODE=required
OCRT_SURFACE_PERSIST_CACHE_DIR=<CACHE_DIRECTORY>
OMP_NUM_THREADS=1
```

No row waits for another row's computation or cache build. A cache miss in `required` mode fails
immediately.

## 5. Validation status

Confirmed on the reference tree:

- clean release build PASS
- CTest 2/2 PASS (reference tree with Mie data)
- cache off / build / required: full CSV byte-identical
- pure-water / red-clay 624-geometry Rrs I/Q/U scatter: exact 1:1
- serial / OpenMP / process rows: byte-identical results and real overlap confirmed
- EAP selector fail-loud

The source, build, test and script files of this code-only package are byte-identical to the
reference tree; see `_RELEASE/SOURCE_CODE_EQUIVALENCE_SHA256.csv`.

## 6. Package self-check

```bash
./verify_package.sh
```

or:

```bash
python _RELEASE/tools/verify_package.py
```

---

# OCRT 속도향상 최종 코드 패키지 — Mie 자료 제외

**패키지:** `OCRT_SPEED_OPTIMIZED_PLATFORM_NEUTRAL_NO_MIE_2026-08-26`
**기준 정본:** `OCRT_C_FINAL_PLATFORM_NEUTRAL_PERF_2026-08-26`
**기준 archive SHA-256:** `0debb9b92f652ee6a9621b28ff1ce40cbeb0d566b8f0080b3c9c1c8580cab629`

## 1. 패키지 범위

승인된 C 성능개선 정본을 그대로 포함한다.

- T_wa 입사-U Q-fold 수정
- 공기 쪽 Cox–Munk 전 모드 푸리에 커널 캐시
- 실제 푸리에 `m_max` 빌드 상한
- 공기–물 보간 정렬/중복제거 맵 메모이제이션
- 버전 관리되는 영속 `R_aa/R_ww/T_aw/T_wa` 계면 연산자 캐시
- worker 국소 `p2b_row` 캐시, 행 사이 lock 제거
- 영속 캐시 코어의 플랫폼 중립 ISO C 파일 입출력
- FR631 reader 호환(631 각도 × 771 파장)
- EAP 생산 런타임 OFF, 명시 선택 시 fail-loud

추가 솔버 최적화는 동결 상태다.

## 2. Mie 자료 제외

요청에 따라 패키지 안에 **`.mie` 파일은 없다**.

- 생산에서 쓸 수 있는 표준 non-EAP 자료: 외부 181개
- EAP 표준 자료: 17개, 보관용. 생산에 설치·로드하지 않음
- 최신 표준 구조: 330–1100 nm, 1 nm(771 파장 열), FR631 631 각도, P11/P12/P33

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

검증기는 181개 non-EAP 파일의 크기와 SHA-256 을 확인하고, 생산 EAP `.mie` 가 설치되어 있으면 실패한다.
(공개 저장소에서는 181개 표를 GitHub Release `data-v1` 로 배포한다. 루트 `README.md` 참조.)

## 3. 빌드

Mie 파일은 컴파일 단계에 필요하지 않다.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

이 패키지를 Mie 없이 별도로 빌드한 결과 `cmake --build` 가 통과했고 `ocrt_nadir_smoke` 1/1 이 통과했다.
`fr631_mie_reader` 와 에어로졸/광물/쇄설물 RT 테스트는 외부 Mie 설치 후 수행한다. 빌드 로그는
`_RELEASE/build_audit/` 에 있다.

## 4. 영속 표면 캐시와 행 독립성

캐시는 시뮬레이션 배치 전에 별도로 미리 만든다.

```bash
python scripts/prebuild_surface_cache.py   --ocrt <OCRT_EXECUTABLE>   --cache-dir <CACHE_DIRECTORY>   -- <OCRT_ARGUMENTS>
```

각 시뮬레이션 행:

```text
OCRT_SURFACE_PERSIST_CACHE_MODE=required
OCRT_SURFACE_PERSIST_CACHE_DIR=<CACHE_DIRECTORY>
OMP_NUM_THREADS=1
```

각 행은 다른 행의 계산이나 캐시 생성을 기다리지 않는다. `required` 모드에서 캐시 miss 는 즉시 실패한다.

## 5. 검증 상태

기준 정본에서 다음을 확인했다.

- clean release build PASS
- CTest 2/2 PASS — Mie 포함 기준 정본
- 캐시 off/build/required 전체 CSV 바이트 동일
- 순수해수/Red-clay 624 기하 Rrs I/Q/U 산포도 정확히 1:1
- serial/OpenMP/process 행 결과 바이트 동일, 실제 겹침 실행 확인
- EAP 선택 fail-loud

이 코드 전용 패키지의 source/build/test/script 파일은 기준 정본과 바이트 동일하며,
`_RELEASE/SOURCE_CODE_EQUIVALENCE_SHA256.csv` 에서 확인할 수 있다.

## 6. 패키지 자체 검증

```bash
./verify_package.sh
```

또는:

```bash
python _RELEASE/tools/verify_package.py
```
