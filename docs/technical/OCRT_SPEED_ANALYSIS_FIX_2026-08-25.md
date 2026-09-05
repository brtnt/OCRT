# OCRT 속도 분석·개선 보고 및 후속 로드맵 — 2026-08-25

첨부 패키지 `OCRT_STAGE2_FAIR_SPEED_PARITY_20260823.zip`(속도 관련 내용만 참조; 결과값 불참조)의 공정 비교 계약을 기준으로 현재 트리(MIGRATION_PKG_2026-08-19, T_wa Q-fix 반영판)를 실측·프로파일링하고, 지배 병목 2건을 수정해 **공정 프로파일 3.25×, 정례 케이스 2.53× 가속**을 얻었다. 동일 머신·근사 공정 조건에서 OCRT가 OSOAA보다 느리던 관계(1.76×)가 **1.85× 우위로 반전**됐다.

## 1. 첨부 패키지 검토 (속도 관련 요지)

- 구판 4.991×(OCRT 84.1s vs OSOAA 16.9s) 판정은 조건 불일치(PSSA ON, 대기 400층, 수중 quadrature 96, 차수·기체·캐시 상태 상이)로 **철회 상태(NON_COMPARABLE_DIAGNOSTIC)** — 병목 참고자료로만 유효.
- 공정 계약: 공통 물리(555 nm/SZA 50/wind 5/red_clay 0.5/M80C AOD865 0.1/PSSA OFF/기체 0) + 공통 해상도(Gauss 48, 대기 26층/해수 640층, 차수 cap 100, particle Fourier m=32, 출력 18×72=1296기하) + OMP 1 + fresh/warm 분리. OSOAA 내부 rough-surface OS_NS=96 vs OCRT m≤32는 `STRICT_INTERNAL_FOURIER_PARITY=false`로 명시 유지.
- OSOAA 공통해상도 probe 17.21s(패키지 측 Linux). 패키지의 ADVANCED 노브 패치(--n-layers-water 등)는 구판용이며, 현재 트리는 동일 지점을 `--debug-water-n-layers/--debug-water-max-orders`(OCRT_DEBUG=1)로 재현 가능 — 신규 패치 불필요.

## 2. 수정 전 실측 (이 샌드박스, 2코어, OMP_NUM_THREADS=1)

| 케이스 | 구성 | wall |
|---|---|---|
| B1 정례(cross-check red_clay) | 555/SZA40/w3, n_mu_water 48, 자동층, 13×24기하 | 14.18 s |
| B2a 공정-근사(에어로졸 제외) | 계약 §2·§4 재현(48/48, 26/640층, 차수100, m32, 18×72) | **45.98 s** |
| B2b = B2a + M80C AOD 0.1 | 동일 | 48.14 s |
| OSOAA 공정-근사 대응(무에어로졸) | NbGauss 48/Mie 16/IGmax 100/NT 26·640, 동일 red_clay·기하 조건 | **26.1–26.3 s** (표면행렬 fresh 확인) |
| OMP 2 (B2a) | — | 44.9 s (이득 2.4% — 단일 케이스 경로 사실상 직렬) |

→ 수정 전 OCRT/OSOAA = **1.75×** (구판 4.99×보다 이미 개선된 상태였음). 참고: fullgrid는 native coupled LUT(1회 해석+재구성) 구조로 기하당 재해석 문제는 이미 없음(`calls=3 water_cold=1`).

## 3. 프로파일 (gprof, B2a): 병목은 물리 코어가 아니라 표면 커널 생성이다

| 순위 | 함수 | 비중 | 호출수 | 진단 |
|---|---|---|---|---|
| 1 | `surface_R_coxmunk_trig` | **43.0%** | 1.62억 | 공기측 거친면 반사 커널의 φ-적분을 **푸리에 모드 m마다 재수행** (33모드 × 동일 기하·Fresnel) |
| 2 | `aw_interp_on_unsorted` | **16.5%** | 2,460만 | 매 호출마다 동일 방향격자를 삽입정렬+dedup (O(n²)/호출) |
| 3 | `rt_sos_operator_apply_vector` | 7.9% | 414 | 실제 RT 물리 코어 — 전체의 8%뿐 |
| 4 | `ocrt_fkc_serve` | 6.8% | 2,183 | all-m 커널 캐시 서빙(복사)·빌드 |
| 5 | air 커널 드라이버 + T_aw/T_wa trig | ~12% | — | 표면 커널 계열 합계 ~62% |

핵심 발견: v1.09-opt에 **all-m 일괄 커널 캐시(FKC, OSOAA SURF_MATR 유사 설계)가 이미 존재**하고 R_ww·T_wa에는 적용돼 있으나, **최대 비용인 공기측 R 커널만 미적용**(per-m 재적분 + 별도 per-m 캐시)이었다. 또 FKC는 요청 m_max와 무관하게 65모드 전부를 빌드했다.

## 4. 적용 수정 (v1.11-speed S1·S2, 5파일, diff 312줄 동봉)

- **S1a — 공기측 R 커널 FKC 통합** (`shared/surface.c`): `surface_coxmunk_fourier_kernel`이 8-슬롯 LRU FKC(키: 격자·풍속·n·nphi·규약 바이트 대조)에서 서빙. m>0 IU/QU 부호접기는 기존 규칙 그대로 서빙 후 적용. 회귀 스위치 **`OCRT_AIR_FKC_OFF=1`= 구경로 비트 재현**. 대기 m-루프가 grid-grid/grid-solar 키를 교대 요청하므로 다중 슬롯이 재빌드 스래시를 방지.
- **S1b — FKC 빌드 모드 상한 힌트** (`surface_fkc_set_build_mmax`, 3개 m-루프 앞 설정): 65모드 고정 빌드 → 실제 m_max로 캡. R_ww·T_wa에도 공통 이득. Rayleigh-only(m_max=2) 케이스의 과빌드 회귀 방지.
- **S2 — `aw_interp_on_unsorted` 정렬 메모화** (`rt_air_water_coupling.c`): 동일 격자 내용(바이트 대조)의 정렬/dedup 인덱스맵을 스레드-로컬 메모화, 값은 맵 경유 수집 — 보간 입력이 직접 경로와 비트 동일.

## 5. 수정 후 실측

| 케이스 | 수정 전 | 수정 후 | 배율 |
|---|---|---|---|
| B2a 공정-근사 | 45.98 s | **14.13 s** | **3.25×** |
| B2b (+M80C) | 48.14 s | **14.87 s** | 3.24× |
| B1 정례 | 14.18 s | **5.61 s** | 2.53× |
| PSSA 대기 스모크(555/SZA75/w5) | 0.39 s | 0.15 s | 2.6× |
| **동일 머신 OSOAA 대비 (B2a)** | 1.76× 느림 | **1.85× 빠름** (14.1 vs 26.1 s) | — |

주의: OSOAA 대비 수치는 이 샌드박스(Linux, 동일 CPU, 단일스레드, 표면행렬 fresh 동등)의 공정-근사 조건 값이다. 계약상 공식 판정은 native Windows + 패키지 러너(PARITY_STATUS 확인 + 1,296행 IQU 산포도 봉인)로 수행해야 하며 이는 로컬 실행 항목이다(§7).

## 6. 정합성 검증 (전 항목 수행 완료)

| 게이트 | 결과 |
|---|---|
| 메모·힌트 비트투명성 (B1을 `OCRT_AIR_FKC_OFF=1`로 실행 vs 수정 전 출력) | **비트 동일** — 편차는 air FKC 단독 기원임을 격리 증명 |
| wind=0 (평면 경로) | FKC on/off **비트 동일** |
| B2a 수치 편차 | 13/1296행, 절대치 max 1e-16(Rrs_U)·1e-17(ρ_Q) — 마지막-ULP 급, 수중장·I열 비트 동일 |
| B1 / B2b 편차 | 7/312행 max 1e-17 / 13/1296행 max 1e-14 — 동급 |
| PSSA 스모크 | ρ_I **비트 동일**, Q/U 6행 1e-22. 레거시 SHA `98e4b2433be300eb…`은 `OCRT_AIR_FKC_OFF=1`로 **자릿수까지 재현**; 신규 기본경로 SHA `2ec2592be03ec88b…`(비트 기준 재설정) |
| 공식 하니스(run_cross_check_fr631.sh) | pure PASS +0.245% · rayleigh_bfo PASS 0.162% · red_clay κ-별건 −0.258/+0.659% — **세 수치 모두 수정 전과 인쇄 자릿수까지 동일** |
| T_wa 단위 게이트 | (T_I,T_Q2,T_U2)=(1,1,1) PASS |
| Q-fix 검증기(작업지시서 동봉) NEW-BUILD | **6/6 PASS**, 통계 동봉 기준과 동일 (pure Q 0.127, red_clay Q 0.085 등) |

ULP 편차의 성격: `-fassociative-math` 빌드에서 φ-적분 누적 루프의 구조가 바뀌면 컴파일러 벡터화 재결합이 달라져 마지막 비트가 흔들릴 수 있다(물리 게이트 0.15~0.6% 대비 11~13자릿수 아래). 완전 비트 재현이 필요한 회귀 비교는 `OCRT_AIR_FKC_OFF=1` 경로를 쓰면 된다.

## 7. 남은 병목과 로드맵 (수정 후 프로파일 근거, 우선순위순)

수정 후 1위는 물리 코어 `rt_sos_operator_apply_vector`(22.2%)로 정상화됐다. 남은 오버헤드성 항목:

- **S3** `aw_interp` 잔여 1.02s(11.9%): 정렬은 제거됐고 남은 것은 호출당 바이트대조+수집. 타깃별 (구간, 가중) 사전계산을 (격자,타깃)당 1회로 끌어올려 성분·모드 간 재사용(≈33×)하면 ~0.1s로. 난도 하.
- **S4** T_aw 계열 1.2s(14%): T_aw는 generic FKC 부적합 이력(2026-08-16, II 오차 기록)이 있어 **전용** all-m 배치 빌더가 필요(자체 memo와 병존). 난도 중.
- **S5** T_wa 뷰노드 일괄화 0.68s: 현재 뷰 μ당 별도 키로 all-m 재빌드(18회) — 뷰노드 배열을 한 키로 묶으면 1회. 난도 중.
- **S6** FKC accumulate SIMD/OpenMP + 포인터 서빙 1.3s: 빌드 내적(모드×9 madd)의 병렬화·복사 제거. 난도 중.
- **S7** 단일 케이스 인트라 병렬화: OMP2 이득이 2.4%뿐 — sos_apply·층 스윕·모드 루프의 OpenMP화가 멀티코어 스케일링의 관건(16코어 로컬에서 배수 기대). 배치 러너는 이미 행 병렬이므로 캠페인은 현재도 코어 스케일함. 난도 상.
- **S8** OSOAA식 표면행렬 디스크 캐시: 동일 (wind, n, 격자, SZA) 재사용 캠페인에서 커널 빌드 자체를 상각. FKC 블록의 파일 직렬화로 구현 용이. 난도 하~중.

추정 여지: S3–S6로 단일스레드 ~10–11s 수준(추가 ~25%), S7로 멀티코어 배수.

- **공식 Windows 공정 게이트(사용자 실행 항목)**: 첨부 패키지 러너(RUN_FAIR_SPEED_GATE_WINDOWS_ANACONDA.ps1)를 새 exe로 native 실행해 PARITY_STATUS=PASS + IQU 산포도 봉인과 함께 공식 속도비를 산출할 것. 구판용 패치 적용 단계는 건너뛰고 현재 트리의 debug 노브 대응을 쓰면 된다(§1).

## 8. 반영 위치

- 소스(수정 5파일): `code/OCRT_C/src/shared/surface.c`(SHA 9c41f187…)·`surface.h`·`src/rt_air_water_coupling.c`(ea3327f9…)·`src/rt_solver.c`·`src/rt_water_rt.c`. diff: `docs/technical/speed_patch_2026-08-25.diff`(312줄, 5파일).
- Windows exe 재빌드: `code/OCRT_C/build/ocrt.exe`(x86-64-v3, SHA 94a07316…)·`ocrt_x86-64-baseline.exe`(9e23316f…) — T_wa Q-fix + speed 패치 포함.
- 샌드박스 런타임 재빌드 완료(binary dcce3c10…). 벤치 산출물: `~/ocrt/speedbench/`.
- 첨부 패키지 사본: 샌드박스 `~/ocrt/speed_pkg/`(구판용 — 속도 계약·러너 참조용으로만 유지).
