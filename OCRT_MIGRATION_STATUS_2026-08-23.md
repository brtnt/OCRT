# OCRT 마이그레이션 상태 문서 (2026-08-23)

이 문서가 D:\Claud_Cowork\OCRT의 **유일한 최상위 기준 문서**다. 새 작업 세션은 이 문서만 읽으면 현재 상태에서 이어갈 수 있다.

프로젝트 목표: **OCRT의 Python & C 페어 구현 개발. 구현 정확도 레퍼런스는 OSOAA.**

---

## 1. 폴더 구조 (2026-08-23 재조직 완료)

```text
D:\Claud_Cowork\OCRT\
├─ OCRT_MIGRATION_STATUS_2026-08-23.md   ← 이 문서 (root에는 이것만 유지)
├─ code\
│  ├─ OCRT_C\            C 소스 최신 패치 상태. build\ocrt.exe = Windows 실행파일
│  ├─ OCRT_Python\       pyOCRT v1.2 (2026-08-14 트리 + 08-20/08-21 패치 3파일 적용)
│  ├─ OSOAA\             OSOAA_OCRT_REFERENCE (OSOAA 소스+OCRT 패치+manifest)
│  ├─ campaign_runner\   편광 PSSA 16-worker 캠페인 러너 (2026-08-22 패키지의 러너부)
│  ├─ run_tools\         과거 run 스크립트 (2026-08-04 패키지)
│  ├─ validation_05\     MIGRATION_PKG의 05_VALIDATION (bit_baselines·harness·artifacts)
│  └─ OCRT_C_pkg_tests\  MIGRATION_PKG의 tests
├─ package\              ★ 모든 원본 아카이브의 버전관리 저장소 (내용 변경 금지)
│  ├─ PART00–PART11 zip 12종 (FR631 분광자료 인계 세트 2026-08-20)
│  ├─ OCRT_POLARIZATION_PSSA_16CORE_CORE_2026-08-22.zip (코어: base 아카이브+러너)
│  ├─ OCRT_paper_migration_2026-08-04.tar.gz
│  └─ raw_campaign\      OSOAA underwater 레퍼런스 데이터 v1/v2, workorder code (2026-07-26)
├─ docs\
│  ├─ paper\             01_paper(논문 초안 docx·빌드스크립트) + figures
│  ├─ technical\         fr631_handoff(계약문서·매니페스트), workorders_2026-07,
│  │                     state_docs(04_STATE_DOCS), polarization_legacy_2026-08-04,
│  │                     legacy_build_2026-08-04, README_MIGRATION_2026-08-04.md
│  └─ validation\        bugreports, sandbox_setup(샌드박스 설치검증 JSON)
└─ _xfer\                전송용 임시 청크 — **이 폴더와 안의 파일 전부 삭제해도 된다**
```

루트의 원본 압축파일은 전부 `package\`로 이동했다. 압축파일 삭제 시 `package\`는 유일본이므로 삭제 대상이 아니다.

## 2. 코드 버전 계보 (code\가 담고 있는 상태)

적용 순서: ① MIGRATION_PKG_2026-08-19 기본 트리 → ② truncation·validation-artifact 정책 오버레이(2026-08-21 zip) → ③ FR631 소비자 패치(PART00 2026-08-20: rt_aerosol.c, rt_aerosol_runtime.c, rt_iop_ahn_mineral.c 771열) → ④ Python bb_b_ratio FR631 이식(아래 §6-b).

- C 델타 21파일 + non-Mie 입력 19파일이 stock 대비 변경분이다(전체 목록: docs\validation\sandbox_setup 옆 `sandbox_delta.tar.gz`의 내용 = _xfer에 있음, 삭제 전 필요시 확인).
- `code\OCRT_C\build\` 의 리눅스 바이너리들(ocrt, ocrt_merged 등)은 **패키지 동봉 stock 빌드**라서 최신 패치 미반영이다. 리눅스/WSL에서는 `scripts\build_release_v1.2.sh`로 재빌드할 것.
- Windows 실행파일은 최신 패치 소스에서 크로스컴파일했다:
  - `build\ocrt.exe` : -march=x86-64-v3 (AVX2), -static, OpenMP 포함 → 이 PC 권장
  - `build\ocrt_x86-64-baseline.exe` : 구형 CPU 안전판
  - 컴파일 방법 재현: `build\build_win.bat` + `build\ocrt_win_compat.c` (MSYS2 MinGW64)
  - 주의: exe는 샌드박스(리눅스)에서 wine 없이 실행 검증을 못 했다. 첫 실행으로 `ocrt.exe --version` 확인 필요.

## 3. 데이터 상태

- 계약: 330–1100 nm 1 nm(771), FR631 정확 각도 631점(180→0°), P11/P12/P33 771열, 소비 보간은 파장·각도 모두 선형.
- canonical Mie 198종의 **로컬 설치는 하지 않았다**(약 4.94 GB 중복 방지). 원본은 `package\PART01–10`이며 SHA-256 전수 검증 완료.
- non-Mie 런타임 데이터(gas xsec, water_iop, tsm_ahn, psi_T 등 19파일)는 `code\OCRT_C\inputs\`에 설치했다.
- 로컬 전체 런타임(198종 포함)이 필요하면 WSL에서:
  `code\campaign_runner\runner\SETUP_AND_BUILD_WSL.ps1` 실행, parts 폴더로 `package\`를 지정. (기존 검증 이력: 2026-08-22 PASS)

## 4. 샌드박스(클라우드 개발환경) 상태 — 완전 구축됨

경로 `~/ocrt/` (이 Cowork 세션 유지 중 유효):

- `runtime/MIGRATION_PKG_2026-08-19/` : 공식 설치기 `setup_runtime.py`로 구축한 **완전 런타임**. 217파일 설치, canonical 198 전수, C release 빌드(march=cascadelake) + PSSA full-grid 스모크 PASS.
- **스모크 CSV SHA-256이 2026-08-22 레퍼런스 런과 비트단위 일치** (`98e4b243…c723`) — 샌드박스가 레퍼런스 환경을 정확 재현함을 뜻한다.
- `parts/` 아카이브 13종(SHA 검증), `archives/`, `polar_pkg/`(러너), `handoff_part00/`, `winbuild/`(exe), gcc 13.3 + mingw-w64 + gfortran + Python 3.11(numpy 2.4) 설치 완료.
- OSOAA 소스는 `runtime_src/MIGRATION_PKG_2026-08-19/02_OSOAA/`에 있고 **아직 샌드박스에서 빌드하지 않았다**(다음 단계).

## 5. 검증 요약 (2026-08-23 수행)

| 항목 | 결과 |
|---|---|
| 아카이브 SHA-256 (코어2 + PART00–10) | 전수 일치 (REQUIRED_FILES.json 대조) |
| 설치기 source audit (198 canonical) | PASS (failure 0) |
| C 빌드 + PSSA 스모크 1296행 | PASS, 레퍼런스와 비트 일치 |
| Python FR631 리더 (5 패밀리 대표: SNF·Ahmad paper·AccuRT·EAP·AHN) | 631×771×771 계약·내림차순 각도·cache_key PASS |
| 이식된 bb_b_ratio | 실행 정상, 값 물리적 타당 (aerosol ~0.08, EAP 0.006, mineral 0.004) |
| Windows exe 링크·정적화 | OK (v3 1.96MB / baseline 1.95MB) |

## 6. 미해결·주의 항목 (다음 세션 우선 확인)

- (a) **Python 위상 파장보간 불일치**: `ocrt_py/aerosol.py`의 `phase_matrix_at_wavelength_pchip()`이 실제 PCHIP을 수행한다. C쪽은 08-20 패치 이후 전부 선형(래퍼 `..._pchip`도 linear 위임)이다. 정수 nm 조회는 노드 정확 일치라 GOCI 밴드 계산엔 영향 없음. 분수 파장에서만 미세 차이. C와 동일하게 linear 위임으로 정리할 것.
- (b) **bb_b_ratio 이식 결정 기록**: 런타임 python constituent.py가 구형(파장 PCHIP·각도 PCHIP 200노드·id() 캐시)으로 남아 있어 08-20 핸드오프 구현(선형 파장·native FR631 trapezoid·안정 캐시키)을 이식했다. 이식 시 08-14 트리의 fail-loud 게이트 2줄(require_wavelength, require_query_in_table)은 **유지**했다 — 핸드오프 원본에는 없던 줄이므로 검토 요망. 원본 백업: 샌드박스 `…/ocrt_py/constituent.py.pre_fr631_graft.bak`.
- (c) **핸드오프 리더 테스트 원본은 import 불일치로 실행 불가**(`phase_matrix_at_wavelength_linear` 부재 — 08-21 오버레이에서 함수명 변경). 동등 검증은 인라인 테스트로 수행했고 PASS. (a) 정리 시 테스트도 갱신할 것.
- (d) **8,045-run historical campaign matrix 미복원**: `campaign_runner/matrices/HISTORICAL_RUN_MATRIX_REQUIRED.csv`가 비어 있으면 full campaign은 fail-loud로 중단된다(설계 의도). 복원 경로: `runner/inspect_historical_package.py --archive package/OCRT_paper_migration_2026-08-04.tar.gz` → `build_matrix_from_run_commands.py`.
- (e) **EAP 정책**: 이 판본은 elastic OCRT·EAP 기본 비활성이다. EAP 17종 Mie는 데이터로는 존재(package). 이전에 EAP 제거를 결정한 바 있으므로, C `rt_iop_organic` 쪽 EAP 필수 로드 요구가 남아있는지 ocean 실행 시 확인 필요.
- (f) **Python 실행 요구사항**: numpy ≥ 2.0 (`np.trapezoid` 사용). 사용자 로컬은 Python 3.13 pycache 흔적 있음 — 로컬 가상환경 구성 시 유의.
- (g) exe 실행 스모크는 Windows에서 직접 1회 필요 (§2).

## 7. 다음 단계 (재개 시)

1. 샌드박스에서 OSOAA 빌드(gfortran 준비됨) → OSOAA-OCRT 비교 파이프라인 재가동 (docs\technical\workorders_2026-07\DOC_OSOAA_reference_generation_2026-07-26.md 참조)
2. §6(a) Python 선형보간 정리 + 리더 테스트 갱신
3. historical matrix 복원 → PSSA full campaign
4. 검증 목표 재확인: OSOAA 대비 ±0.5% (GOCI 6밴드), P0 = 고 ω rrs(0−) 잔차 개선

## 8. 작업 이력 (2026-08-23)

로컬→샌드박스 1.72GB 청크 릴레이(8MB×~215, 디바이스 브리지 호출당 8MB 한계 우회) → SHA 전수 검증 → 공식 설치기로 샌드박스 런타임 구축·스모크 비트일치 확인 → mingw-w64 Windows exe 크로스컴파일 → 로컬 폴더 재조직(code/package/docs) 및 패치 델타 적용 → 본 문서 작성. 이후 세션 일시정지(사용자 이동).

---

## 9. 세션 2 결과 (2026-08-23 오후, 샌드박스 단독)

### 9.1 OSOAA 샌드박스 빌드·검증 완료

- gfortran 13.3으로 빌드. **내 빌드 = 동봉 640층 바이너리와 수치 완전 동일**(스모크 케이스 101행 rel diff 0). `-march` 무관.
- 스모크 기준 CSV(`reference_values_Z09_20260728.csv`)는 **80층 시절 stale**임을 확정: `exe80`(=CTE_NT_SEA 80 보관본, TODO#4의 원복물)만 통과하고 모든 640 빌드는 rel ~1.4e-4로 "FAIL"한다. → 640 기준 `reference_values_Z09_NT640_20260823.csv`를 **추가**(구판 보존, check_smoke 11/11 <1e-9 PASS). 스모크는 이제 640 기준으로 돌릴 것.

### 9.2 OCRT-OSOAA 교차 하니스(05_VALIDATION/harness) 재가동 — red_clay −70% 사건 규명·해소

- 최초 실행: pure_d200 **PASS(+0.005%)**, rayleigh_bfo **PASS(0.162%)**, red_clay_tsm555 **FAIL(−70%)**.
- 원인 확정: 설치기가 stock Red_clay_AHN.mie(구판: 16열 sparse, bb/b(555)=0.02168)를 **FR631 canonical**(신판: 771열, 해수배경 forward Mie 재계산, bb/b(555)=0.01000, 전방피크 4배 예리)로 교체 — 저장 OSOAA 참조(`runs/nt640_osoaa_tsm555`)는 구판 위상으로 생성된 stale 기준이었다. 코드 회귀 아님(direct/moment 커널 동일 결과로 커널 가설도 기각).
- 저장 기준의 생성 조건을 완전 복원(재현 관문: stock 위상 재생성이 저장본과 Flux·radiance 전 자릿수 일치, c(555)=0.457744·depth 36.343m·MOT 0.0935485·NbGauss48/100 확정) 후 **canonical 위상으로 신규 기준 생성**: `runs/nt640_osoaa_tsm555_fr631/`(생성 스크립트 `05_VALIDATION/tools/regen_osoaa_tsm555_reference.sh`).
- 신규 기준 대비: **rrs 중앙 −0.258%, Rrs 중앙 +0.642%** (직전 −70%에서 수렴; OSOAA 신/구 기준 radiance 비 0.305가 위상 교체 효과를 정량 확인).
- 잔여 −0.26%의 소재 판별: OSOAA L=200→600 절단 무이동, CTE_NT_SEA 640→1280 무이동 → **κ-계열 합성 잔차 범주**(기존 12′ 결론 "≤0.41%, 논문 기술 방침·추가 환원 중단 권고"와 정합). 단 게이트(rrs 0.15%/Rrs 0.60%)는 그대로 두었으므로 공식 기록은 `run_cross_check_fr631.sh` 기준 2/3 PASS + red_clay 여유초과 FAIL이다. 게이트 재조정 여부는 사용자 판단 사항.

### 9.3 Python 파장보간 linear 정리 (§6a 해소)

- `ocrt_py/aerosol.py`에 `phase_matrix_at_wavelength_linear`(FR631 계약: 구간선형, 정수 nm 정확 조회) 신설, 기존 `..._pchip`는 C `mie_io.c` 패턴대로 linear 위임 래퍼로 전환(DOC-REF 주석 포함). 정수 nm 정확조회·분수 파장 bracketing 검증 통과.
- **공식 핸드오프 리더 테스트(`test_fr631_python_reader.py`) 5개 패밀리 전 항목 PASS(exit 0)** — §6(c)도 함께 해소.

### 9.4 Historical 8,045-run 행렬 복구 — 대부분 해소

`inspect_historical_package.py`가 paper_migration에서 `runs_summary.csv`(8,045행, 전 행 status=ok) 추출: **aer_model 24종, AOD865 7값(0.02~0.6), 밴드 412/555/667/748/865, 풍속 3/5/9, 수체상태 65종(chl/tsm/cdom), branch N/S/A/NS, 케이스별 코드버전·바이너리 SHA·통계** — blocker 문서의 미확보 항목 대부분이 파라미터 수준에서 복구됐다. 미복구는 **run별 정확한 command line/advanced 옵션**뿐. 파라미터→명령행 재구성은 현행 정책(PSSA ON 등) 주입이 필요한 **정책 결정 사항**이라 사용자 승인 대기로 남긴다. 추출물: 샌드박스 `polar_pkg/.../historical_import/` (+ atm_reference.csv.gz 42.8MB, cases_thin.csv.gz, underwater_validation.csv).

### 9.5 로컬 동기화 대기물 (재접속 시 커밋)

세션 종료 시점에 로컬 디바이스 오프라인이면 다음이 샌드박스에 대기한다: `~/ocrt/sandbox_delta2_2026-08-23.tar.gz`(14KB: aerosol.py 정리본, fr631 하니스 2종, regen 스크립트, NT640 스모크 기준 CSV). 적용 위치: code/OCRT_Python/ocrt_py/, code/validation_05/harness·tools/, code/OSOAA/verification/. 신규 OSOAA 기준 run(55MB)은 대용량이라 미포함 — regen 스크립트로 WSL 재생성 가능(재현 관문 검증됨).

### 9.6 갱신된 미해결 목록 (§6 대체)

- (신규-1) 교차 하니스 red_clay 잔여 −0.258%/+0.642% — κ-잔차 판정. 게이트(0.15/0.60) 유지 중이므로 게이트 재조정 또는 잔차 추가 환원 여부 결정 필요. P0(고 ω rrs 잔차)와 함께 다룰 것.
- (신규-2) historical matrix: 파라미터 복구 완료, 명령행 재구성 정책 승인 대기 (§9.4).
- (신규-3) OSOAA 스모크 스크립트가 기본으로 구판 CSV를 참조 — NT640 CSV로 전환하는 한 줄 수정 대기(원본 스크립트 보존 정책 때문에 미수정).
- (유지) §6(b) bb_b_ratio 이식 게이트 2줄 검토, §6(e) EAP 정책 확인, §6(g) Windows exe 실행 1회 확인, §6(f) numpy≥2.0.
- (해소) §6(a)·(c) → §9.3. red_clay −70% → §9.2. OSOAA 빌드 → §9.1.

---

## 10. 정본 교체 (2026-08-26): OCRT_C_FINAL_PLATFORM_NEUTRAL_PERF

- **code\OCRT_C는 2026-08-26 정본**(`package\OCRT_SPEED_OPTIMIZED_PLATFORM_NEUTRAL_NO_MIE_2026-08-26.zip`, sha b19b44c2…)으로 교체됐다(1,289파일 전개, inputs·build 제외 후 exe 재컴파일 커밋). 결과 불변 검증 전 항목 통과 기록: `docs\technical\OCRT_CANONICAL_UPDATE_VERIFICATION_2026-08-26.md`.
- 포함: T_wa Q-fold fix, air-side all-m FKC + m_max cap, 보간 메모화(이상 기존 반영분), **persistent 계면연산자 디스크 캐시(기본 off; prebuild→required 운용)**, p2b 행간 락 제거, FR631 reader 호환, EAP OFF fail-loud, in-tree Windows compat.
- 기본 모드 산출은 직전 판과 **비트 동일**(PSSA 스모크 SHA 2ec2592b…/레거시 98e4b243… 재현). warm pcache에서 추가 ≈2×.
- **§2의 Windows 빌드 방법 갱신**: 외부 `ocrt_win_compat.c` 불필요(in-tree `rt_windows_compat.h`). exe: `build\ocrt.exe`=v3(sha b9bea815…), baseline 739af90b…. Windows 실행 1회 스모크는 여전히 필요(§6g).
- **EAP mie 17종은 production 트리에서 `package\EAP_MIE_ARCHIVAL_ONLY_17\`로 이동**(패키지 정책: archival only). 외부 mie 검증기 181/181 PASS.
- 샌드박스 런타임도 동일 정본으로 동기화·재빌드 완료(바이너리 비트 동일).

## 11. Mie 자료 규격 통일 (2026-08-26)

`code\` 전체 .mie를 **FR631 정본 규격(631 각도 × 771 파장, 330–1100 nm 1 nm)** 으로 통일했다. 최종: **363개 전부 규격 일치, 비규격 0**.

- **`code\OCRT_C\inputs\` (181개)는 이미 정본이었다** — 패키지 매니페스트(`_RELEASE\external_mie\EXTERNAL_MIE_RUNTIME_REQUIRED_181.csv`) 대비 size 181/181 일치, SHA-256 표본(Red_clay 포함 4건) 일치. 교체 불필요.
- **`code\OCRT_Python\data\`가 구형이었다** — 199개 전부 361-각도 계열(구 22/101/771 파장 혼재, 합 60 MB). C inputs의 정본 181개로 전량 교체(4.21 GB), manifest size 181/181 일치, Red_clay SHA-256 일치 확인.
- **비규격 24개 삭제**: `inputs\deprecated\M50C.mie`(83각×20파장; 정본 `inputs\M50C.mie` 별도 존재), `water_iop\candidates\...CANDIDATE_REJECTED_L200.mie`(C·Py 2개), `validation\mie_candidates\...rejected_20260816\`(4개), `OCRT_Python\data\water_iop\eap\`(구 EAP 17개). 삭제 전 전량을 `package\LEGACY_NON_FR631_MIE_2026-08-26\`에 원본 경로 구조로 보존했다(약 7 MB) — 완전 제거를 원하면 이 폴더를 지우면 된다.
- EAP 정책 유지: production 트리에 EAP 자료 없음(정본 FR631 EAP 17종은 `package\EAP_MIE_ARCHIVAL_ONLY_17\`).
- 영향: `task\` 검증 기준 CSV는 C inputs 기반이라 무효화되지 않는다(입력 불변). Python 페어는 이제 C와 동일 자료를 읽는다.

## 12. TSM mineral a*/b* Mie-graft 개정 반영 (2026-08-29)

전달 패키지 `package\OCRT_TSM_astar_bstar_MIEgraft_2026-08-29_KST.tar.gz`의 8파일(4광물 × a*/b*)을 검증·반영했다. 상세: `docs\technical\OCRT_TSM_ASTAR_BSTAR_UPDATE_RECORD_2026-08-29.md`.

- 개정 범위: 330–400 nm(구 OLS→Mie graft), 750–1100 nm(.mie bulk 유래→Mie graft), 400–450·700–750 nm blend 창 신설. **450–700 nm와 >1100 nm는 bit-identical 불변.**
- 검증: 패키지 gate 8/8 PASS, SHA256SUMS 34/34, `data_original`이 우리 기존 설치본과 8/8 일치, 450–700 nm 재검산 상대차 0.0, red_clay 555 nm a*=0.027810·b*=0.764760 불변.
- **OCRT 실행 검증: 555 nm 정례 케이스 산출이 기준과 비트 동일** → `task\` 기준 CSV·교차검증 게이트·속도 벤치 전부 유효 유지. 865 nm 대조군은 b_total +2.02%로 신 파일 소비 확인.
- 반영: `code\OCRT_C\inputs\tsm_ahn\`·`code\OCRT_Python\data\tsm_ahn\` 각 8파일(SHA 8/8 일치), 원본 tar.gz는 `package\`에 보관.
- **`.mie` bulk 컬럼 re-mirroring은 미수행이나 production 영향 없음** — `rt_iop_ahn_mineral.c` 확인 결과 `.txt`가 유일한 production 출처이고 `.mie` bulk는 `.txt` 부재 시 fallback 전용이다.
- 후속: 330–399 nm 채움 방식 최종 결정, `.mie` re-mirroring, 문서 정정 3건(H5.0 spline·zero-crossing 수치 폐기, miepython 원인 규명, Junge ξ=4 부적합), R_eff를 위상함수·bb/b에 사용 금지.

## 13. S_CDOM / S_detritus 제어 옵션 (2026-08-29)

상세: `docs\technical\OCRT_CONSTITUENT_SLOPE_OPTIONS_2026-08-29.md`, diff: `docs\technical\py_constituent_slope_options_2026-08-29.diff`.

- **C는 기존부터 지원**: `--ocrt-adom-slope`(기본 0.014 nm⁻¹), `--ocrt-detritus-slope`(0.0109), 짝으로 `--ocrt-adom440`·`--ocrt-detritus-a440`. `--help` 문서화·게이트 없음·배치 CSV 열 지원. detritus는 `--ocrt-chl>0` 필요.
- **Python 페어에 옵션 추가**(3파일): `ocrt_solve.py`에 `--adom-slope`(별칭 `--cdom-slope`)·`--detritus-a440`·`--detritus-slope` + grid CSV 선택 열(C 철자 수용) + 출력 열 기록, `batch_driver.py`에 슬로프 전달·혼합 슬로프 폴백, `produce_grid.py` 열 전달. 기본값은 C와 동일.
- 검증: C↔Py IOP 차 ≤2.5e-7, 슬로프 응답비 상대차 0.0004%(Rrs)·0.0002%(rrs), **기본값 실행은 패치 전과 자릿수 동일(무회귀)**.
- 후속: `evaluate_batch`의 슬로프 벡터화(대량 스윕 시), C detritus의 Chl>0 제약 재검토 여부.

## 14. 구면 보정 교체: PSSA → IPSS 단일화 + VZA 지상앵커 (2026-09-05)

근거: Zhai & Hu (2022), JQSRT 282, 108132. 상세: `code\OCRT_C\validation\ipss_2026-09-05\OCRT_IPSS_REPLACEMENT_RECORD_2026-09-05.md`.

- **IPSS 신규 구현**: `src/rt_ipss.{c,h}` — 파라메트릭 레이-셸 추적(지표 차폐 판정 포함), 구면/천저/평면평행 3규약 단일산란 적분, Eq.(7) κ = I₁,ss/I₁,pp, 정확 구면 직달빔 광학두께. κ는 (I,Q,U)를 같은 배로 곱하는 스칼라 재척도이며 SOS 해는 평면평행을 유지한다.
- **구 PSSA 완전 삭제**(사용자 지시): `src/rt_pssa.{c,h}`, `tests/test_pssa_geometry.c`, `scripts/test_pssa_internal.sh` 제거, `rt_options_t.pssa_mode`·`rt_pssa_mode_t`·`--pssa-mode` 제거, `rt_atm_t`의 `pssa_*` 필드와 `rt_surface_boundary.c`의 PSSA 분기 제거. `--pssa`는 on/off 스위치 하나만 남았고 `--pssa-mode`를 주면 exit 2. 삭제 파일은 `validation/ipss_2026-09-05/removed_legacy/`에 보존.
- **VZA 앵커 결함 수정**: 초판은 `--vza`를 논문 §3 규약대로 TOA 천정각으로 해석했는데, 그 앵커에서는 태양이 천정에 있어도 −2.9 %(VZA 55)~−11 %(VZA 70)의 SZA 무관 편차가 남고 대기를 지표로 붕괴시켜도 사라지지 않는다. 지상 화소 기준(위성 L1B senz)으로 고쳐 모든 평면평행 극한에서 보정량이 0으로 수렴하게 했다. 게이트 G-I3가 이 극한을 상시 검사한다.
- **영향**: `--pssa` 산출이 이전 IPSS 대비 ρ_I +2.6~+3.1 % 변한다(이전 값이 틀렸다). `--pssa` 없는 산출은 **비트 불변**(10케이스 배터리, 비-PSSA 스모크 SHA `4ced4f18d50e4ac9` 동일). `T_dir_dn`도 불변.
- **구면 보정이 실제로 필요한 구간**: VZA 55 / RAA 90 / M80C AOD 0.1 / BFO w5 기준으로 무보정 오차가 SZA 65°에서 0.06 %, 70° −0.22 %, 80° −2.27 %, 85° −12.13 %(555 nm). **SZA 75° 이상 슬롯이 실질 영향권**이다. 구 legacy PSSA의 잔차는 이 단면에서 −0.4~+0.9 %지만, 전 방위·VZA 격자에서는 SZA 85°에 −1.4~+5.1 %이고 주평면을 가로지르며 부호가 바뀐다(legacy가 가질 수 없는 자유도).
- **게이트**: `scripts/run_ipss_gates.sh` — G-I0(기하 vs Python 오라클) / G-I1(논문 T9 앵커) / G-I1b(구적 수렴) / G-I2(평면평행·방위·가드) / **G-I3(앵커 극한)** + `--pssa-mode` 거부 + 소스 트리 legacy 심볼 부재 + `--pssa` on≠off + 비-PSSA 불변 + 결합해양 fail-loud. 샌드박스·로컬 양쪽 ALL PASS.
- **결합 해양(`--water-model`) + `--pssa`는 fail-loud**(exit 1). Phase I 범위는 대기 단독/흑해양이며 legacy 되돌림 경로가 없다.
- **비용**: 1,260기하 fullgrid, n_layers=40, OMP 1 — off 0.045 s → on 0.112 s.
- **빌드 주의**: 이 PC의 CPU(Core Ultra 9 275HX)는 AVX-512 미지원이라 Makefile 기본값 `-march=cascadelake` 산출물은 즉시 Illegal instruction으로 죽는다. `-march=x86-64-v3`로 빌드할 것. Windows 재빌드 절차는 `code\OCRT_C\REBUILD_WIN_v1.11_KO.txt`.
- 후속: Windows `.exe` 재빌드(사용자 실행 필요), 과거 `--pssa` 캠페인 산출물 재실행 대상 판정(SZA 75° 이상), Phase II 결합 해양 Eq.(7) 적용 정책.

### 14.1 v1.11.1 — 혼합 대기 κ 위상함수 가중 수정 (2026-09-05 재조사)

상세: 기록 문서 §11.6. 변경 파일 `src/rt_ipss.{c,h}`, `src/rt_solver.c`, `tests/test_ipss_gates.c`; `SHA256SUMS.txt` 재생성(1,724건).

- **결함**: v1.11 의 κ 는 층 원천을 β_j·ω_j(위상 미가중)로 두었다. Eq.(7)의 I₁ 은 실제 단일산란 radiance 이므로 위상함수는 Rayleigh/에어로졸 혼합비가 고도에 무관할 때만 상쇄된다. 운영 대기(에어로졸 2 km / Rayleigh ~8 km)에서는 저고도 에어로졸 층이 P_R/P_A ≈ 8–10 배 과대 가중됐다. RTSOS 는 층별 위상행렬로 I₁ 을 만든다(논문 정의).
- **수정**: 층별 산란분율(`ydel/xdel`)과 P_R(Θ)·P_A(Θ)(FR631 값표 우선, Legendre 모멘트 대체)로 w_src,j = w_R,j P_R + w_A,j P_A 를 기하마다 구성해 구면·평면평행 두 적분에 적용. Rayleigh 단독 프로파일은 비트 동일, 비-PSSA 경로 불변(SHA `4ced4f18d50e4ac9`).
- **영향(VZA 55 / RAA 90 / M80C AOD 0.1 / BFO w5, IPSS/무보정 − 1)**: 저 SZA 평탄부 412/555/865 = −0.159/−0.101/−0.087 % → −0.191/−0.169/−0.176 %; SZA 85 = +9.36/+13.80/+2.56 % → +8.77/+11.91/+1.89 %. 즉 §14 의 "무보정 오차 SZA 85° −12.13 %(555)" 는 −10.64 % 로 정정.
- **게이트 추가**: G-I4a(혼합비 고도-일정 시 미가중과 동일, 2.4e-15), G-I4b(2성분 합성 프로파일 κ vs 독립 scipy 적분, |Δκ| ≤ 7.5e-6). 전 게이트 PASS.
- **진단**: `OCRT_IPSS_PROFILE_DUMP=<file>` 환경변수로 κ 입력 프로파일·P11 표 덤프; `validation/ipss_2026-09-05/phasew/`.
- **후속**: Windows `.exe` 재빌드 필요(`code\OCRT_C\REBUILD_WIN_v1.11.1_KO.txt`). 에어로졸 대기에서 Eq.(7) 근사 자체의 정확도(전체 구면 RT 대비)는 여전히 미검증 외삽.
