# OCRT 정본 교체 검증·반영 기록 — 2026-08-26

대상: `OCRT_SPEED_OPTIMIZED_PLATFORM_NEUTRAL_NO_MIE_2026-08-26.zip` (sha256 b19b44c2bf3b288c…, 정본 트리 `OCRT_C_FINAL_PLATFORM_NEUTRAL_PERF_2026-08-26`). 조건: **기존과 결과값 불변, 속도만 개선** 확인 시 정식 반영. 판정: **통과 — 반영 완료**.

## 1. 패키지 내용 (기존 대비 변경: src 17파일 + 신규 2모듈)

기존 반영분(T_wa Q-fold, air FKC, m_max cap, 보간 메모)에 더해: versioned persistent R_aa/R_ww/T_aw/T_wa 계면연산자 디스크 캐시(신규 `surface_persistent_cache.c/h`, ISO C I/O, 기본 **off**), p2b row 캐시 worker-local화(행간 락 제거), FR631 reader 호환, EAP runtime OFF fail-loud, in-tree Windows compat(`rt_windows_compat.h` — 외부 compat 파일 불필요). Mie payload 0개(외부 설치형).

## 2. 결과 불변 검증 (샌드박스, 전 항목 통과)

| 게이트 | 결과 |
|---|---|
| B1 정례 red_clay fullgrid (기본 모드) | 현행 트리와 **비트 동일** |
| B2a 공정-근사 1296기하 (기본 모드) | **비트 동일** |
| PSSA 스모크 SHA | 기본경로 `2ec2592b…`·레거시경로(OCRT_AIR_FKC_OFF=1) `98e4b243…` — 기록된 기준 SHA **완전 재현** |
| wind=0 (평면 경로) | **비트 동일** |
| persistent cache REQUIRED(사전빌드 후) | 무캐시 산출과 **비트 동일** (B1·B2a 모두) |
| 캐시 교차설정 가드 | wind5 요청 × wind3 캐시 → content-addressed **miss + fail-loud** (오재사용 없음) |
| OMP 1 vs 2 (p2b 락 제거 검증) | **비트 동일** |
| 공식 하니스 | pure PASS +0.245% / rayleigh_bfo PASS 0.162% / red_clay κ-별건 −0.258·+0.659% — 현행과 자릿수 동일 |
| T_wa 단위 게이트 | (T_I,T_Q2,T_U2)=(1,1,1) PASS |
| Q-fix 검증기 (NEW-BUILD) | 6/6 PASS |
| 외부 mie 검증기 (`verify_external_mie.py`) | 181/181 present, bad_hash 0, production EAP 0 → **PASS** |

## 3. 속도 (동일 컨테이너 연속 측정, OMP 1)

| 케이스 | 현행 | 신판 cold | 신판 warm(pcache) |
|---|---|---|---|
| B1 정례 | 7.1–7.2 s | 7.0–7.2 s (동일) | **3.8 s (≈1.9×)** |
| B2a 공정-근사 | 17.7 s | 17.6 s (동일) | **8.5 s (≈2.1×)** |

기본(cold) 속도는 현행과 동일(코드 경로 비트 동일이므로 당연). 신규 이득은 **사전빌드 캐시 재사용**에서 나온다 — 동일 (풍속·격자·n·SZA) 반복 실행(캠페인)에서 커널 빌드 전체가 상각된다. 사용법:

```
python scripts/prebuild_surface_cache.py --ocrt <ocrt> --cache-dir <DIR> -- <run args>
env OCRT_SURFACE_PERSIST_CACHE_MODE=required OCRT_SURFACE_PERSIST_CACHE_DIR=<DIR> ocrt <run args>
```

미설정 시 off(현행과 완전 동일 동작). readonly=있으면 사용, required=없으면 fail-loud.

## 4. 반영 내역

- 샌드박스 런타임 `01_OCRT_C`: 정본 파일 동기화 + 재빌드 — 패키지 빌드와 **바이너리 비트 동일** 확인.
- 로컬 `code\OCRT_C`: 정본 zip(package\ 보관본)에서 1,289파일 직접 전개(inputs·build 제외), 대표 파일 SHA 대조 일치. `build\ocrt.exe`(=v3, sha b9bea815…)·`ocrt_x86-64-v3.exe`·`ocrt_x86-64-baseline.exe`(739af90b…) 재컴파일 커밋 — **in-tree compat로 외부 ocrt_win_compat.c 불필요해짐**.
- EAP 정책: production 트리의 EAP mie 17개를 `package\EAP_MIE_ARCHIVAL_ONLY_17\`로 이동(패키지 정책: EAP는 archival only, runtime OFF fail-loud). 외부 mie 검증기 PASS 상태.
- 정본 zip은 `package\`에 보관(버전 관리).

## 5. 참고

- 이 교체로 기존 문서의 소스 SHA 기준이 갱신된다: surface.c 9c41f187…(구) → 7b295238…(정본), 신규 surface_persistent_cache.c dd836890….
- 비트 기준: PSSA 스모크 ON/OFF SHA 쌍(2ec2592b…/98e4b243…)은 정본에서도 유효 — 회귀 비교 기준 유지.
- 배치/campaign 러너에서 pcache를 쓸 때는 문서대로 **행 시작 전 단일 prebuild** 후 required 모드로 실행한다(행간 대기·락 없음).
