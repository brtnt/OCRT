# OCRT C — integrated S1/S2 speed patch (2026-08-25)

This tree is the C reference candidate that keeps the following together.

- EAP species scattering: production runtime OFF; an explicit selector fails loudly
- T_wa incoming-U column Q-fold fix
- FR631 reader compatibility (631 angles × 771 wavelengths)
- air-side Cox–Munk all-m Fourier kernel cache
- true Fourier `m_max` build cap
- memoized air–water interpolation sort/dedup map

## EAP data policy

`inputs/water_iop/eap/*.mie` are excluded from this production source archive. The EAP generator,
sources and provenance documents remain as reproduction material, but the production runtime
does not install or load EAP tables.

## Validation

- clean GNU C release build: PASS, 0 warnings
- pure water / red clay, baseline vs patched full CSV: byte-identical
- Rrs I/Q/U scatter: exact 1:1
- T_wa spin-2 gate: PASS
- Q-fold work order: 6/6 PASS
- explicit EAP selector: fail-loud, exit code 2
- FR631 reader: 631 × 771 PASS

## Scope of the speed statement

In a same-host diagnostic, OCRT was faster than the reference model even when the reference
model's fresh surface-matrix generation was included. No operating system is an official gate.
Later comparisons use the same host, the same build class, the same numerical settings, the same
definition of fresh/warm cache and the median of repeated runs, and they are sealed together with
the I/Q/U scatter plots and the exact source CSV.

## Final stage of the persistent surface-operator cache (2026-08-25)

- adds a versioned, content-addressed persistent cache of `R_aa`, `R_ww`, `T_aw`, `T_wa`
- operation: a separate `build` pre-run, then simulation rows use `required` read-only mode
- `build` mode is refused (fail-loud) under `--batch-full-grid`
- `p2b_row` cache made worker-local; `critical(p2b_rowcache)` removed
- no lock, poll or wait; no ordering between simulation rows
- pure water / red clay: cache off / build / required full CSV byte-identical
- Rrs / rrs / TOA I/Q/U scatter: exact 1:1
- serial, 4-thread and 4-process rows: byte-identical results and real execution overlap
- after this stage, further speed work (output-only exact view, contracted operator, subsystem
  R/T/Jacobian) is on hold

Contract and usage: `docs/PERSISTENT_COXMUNK_INTERFACE_CACHE_KO_2026-08-25.md`

## Platform-neutral finalization (2026-08-26)

- native-Windows-only gates, runners and PowerShell procedures removed
- persistent-cache core reduced to ISO C file I/O
- Windows-only prebuild helper removed; `scripts/prebuild_surface_cache.py` is the single tool
- same-host / same-options comparison contract independent of the operating system
- further speed optimization frozen

---

# OCRT C — 2026-08-25 S1/S2 속도 패치 통합본

이 트리는 다음을 함께 보존한 C 정본 후보다.

- EAP 종별 산란: 생산 런타임 OFF, 명시 선택 시 fail-loud
- T_wa 입사-U 열 Q-fold 수정
- FR631 reader 호환(631 각도 × 771 파장)
- 공기 쪽 Cox–Munk 전 모드 푸리에 커널 캐시
- 실제 푸리에 `m_max` 빌드 상한
- 공기–물 보간 정렬/중복제거 맵 메모이제이션

## EAP 자료 정책

`inputs/water_iop/eap/*.mie` 는 이 생산 소스 아카이브에서 제외했다. EAP 생성기·원자료·출처 문서는 재현
자료로 남지만 생산 런타임은 EAP 표를 설치·로드하지 않는다.

## 검증

- clean GNU C release build: PASS, 경고 0
- 순수해수 / Red clay 기준본–패치본 전체 CSV: 바이트 동일
- Rrs I/Q/U 산포도: 정확히 1:1
- T_wa spin-2 게이트: PASS
- Q-fold 작업지시: 6/6 PASS
- 명시 EAP 선택: fail-loud, 종료 코드 2
- FR631 reader: 631 × 771 PASS

## 속도 판정 범위

같은 호스트 진단에서 OCRT 는 참조 모델의 fresh surface-matrix 생성 시간을 포함한 조건보다도 빨랐다.
특정 운영체제를 공식 게이트로 지정하지 않는다. 이후 비교는 같은 호스트, 같은 빌드 종류, 같은 수치 설정,
같은 fresh/warm 캐시 정의와 반복 실행의 중앙값으로 수행하며, I/Q/U 산포도와 정확한 원본 CSV 를 함께
봉인한다.

## 2026-08-25 영속 표면 연산자 캐시 최종 단계

- 버전 관리·내용 주소 방식의 `R_aa`, `R_ww`, `T_aw`, `T_wa` 영속 캐시 추가
- 운영 방식: 별도 `build` 사전 실행 후 시뮬레이션 행은 `required` 읽기 전용 사용
- `build` 모드는 `--batch-full-grid` 에서 fail-loud 로 거부
- `p2b_row` 캐시를 worker 국소로 바꾸고 `critical(p2b_rowcache)` 제거
- lock/poll/wait 없음, 시뮬레이션 행 사이 선후 관계 없음
- 순수해수/Red-clay 캐시 off/build/required 전체 CSV 바이트 동일
- Rrs/rrs/TOA I/Q/U 산포도 정확히 1:1
- serial, 4-thread, 4-process 행 결과 바이트 동일, 실제 실행 겹침 확인
- 이 단계 이후 추가 속도 개선(출력 전용 정확 뷰, 축약 연산자, 하위계 R/T/Jacobian)은 보류

상세 계약과 실행법: `docs/PERSISTENT_COXMUNK_INTERFACE_CACHE_KO_2026-08-25.md`

## 2026-08-26 플랫폼 중립 최종화

- native-Windows 전용 게이트·실행기·PowerShell 절차 폐기
- 영속 캐시 코어를 ISO C 파일 입출력으로 정리
- Windows 전용 사전 빌드 도구를 제거하고 `scripts/prebuild_surface_cache.py` 로 통일
- 운영체제와 무관한 같은 호스트/같은 옵션 비교 계약 채택
- 추가 속도 최적화 동결
