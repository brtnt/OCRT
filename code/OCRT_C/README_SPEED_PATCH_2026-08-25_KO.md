# OCRT C — 2026-08-25 S1/S2 속도 패치 통합본

이 트리는 다음을 함께 보존한 C 정본 후보이다.

- EAP species scattering production runtime OFF 및 명시 selector fail-loud
- T_wa incoming-U column Q-fold 수정
- FR631 631-angle × 771-wavelength reader compatibility
- air-side Cox–Munk all-m Fourier kernel cache
- 실제 Fourier `m_max` build cap
- air–water interpolation sort/dedup map memoization

## EAP 자료 정책

`inputs/water_iop/eap/*.mie`는 이 production source archive에서 제외했다. EAP generator/source/provenance 문서는 재현 자료로 남지만 production runtime은 EAP table을 설치·로드하지 않는다.

## 검증

- clean GNU C release build: PASS, warning 0
- pure water / Red clay baseline–patched full CSV: byte-identical
- Rrs I/Q/U scatter: exact 1:1
- T_wa spin-2 gate: PASS
- Q-fold work-order: 6/6 PASS
- explicit EAP selector: fail-loud, exit code 2
- FR631 reader: 631 × 771 PASS

## 속도 판정 범위

동일 host 진단에서 OCRT는 OSOAA의 fresh surface-matrix 생성 포함 조건보다 빨랐다. 특정 운영체제를 공식 Gate로 지정하지 않는다. 이후 비교는 동일 host, 동일 build class, 동일 수치설정, 동일 fresh/warm cache 정의와 반복 median으로 수행하며 I/Q/U 산포도와 exact source CSV를 함께 봉인한다.

## 2026-08-25 persistent surface-operator cache 최종 단계

- versioned/content-addressed `R_aa`, `R_ww`, `T_aw`, `T_wa` persistent cache 추가
- 운영 mode: 별도 `build` pre-run 후 simulation rows는 `required` read-only 사용
- `build` mode는 `--batch-full-grid`에서 fail-loud 거부
- `p2b_row` cache를 worker-local로 변경하고 `critical(p2b_rowcache)` 제거
- lock/poll/wait 없음; simulation row 간 계산 선후관계 없음
- pure water/Red-clay cache off/build/required full CSV byte-identical
- Rrs/rrs/TOA I/Q/U 산포도 exact 1:1
- serial, 4-thread, 4-process row 결과 byte-identical 및 실제 실행 overlap 확인
- 이 단계 이후 추가 속도 개선(output-only exact view, contracted operator, subsystem R/T/Jacobian)은 보류

상세 계약 및 실행법: `docs/PERSISTENT_COXMUNK_INTERFACE_CACHE_KO_2026-08-25.md`

## 2026-08-26 플랫폼 중립 최종화

- native-Windows 전용 gate·runner·PowerShell 절차 폐기
- persistent-cache core를 ISO C file I/O로 정리
- Windows-only prebuild helper를 제거하고 `scripts/prebuild_surface_cache.py`로 통일
- 운영체제와 무관한 same-host/same-options 비교 계약 채택
- 추가 속도 최적화 동결
