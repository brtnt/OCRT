# OCRT 330–1100 nm spectral data integration — read first

Date: 2026-08-14
Target: the paired OCRT C / Python implementation
(The external reference model kept its 2026-07-28 configuration; it was not part of this spectral
extension.)

## 1. Final versions

- C: `OCRT-v1.2-2026-08-14-KST-spectral-data-330-1100-component-gated`
- Python: `pyOCRT-v1.2-2026-08-14-spectral-data-330-1100-component-gated`

## 2. Main result

The handoff contained 197 paired runtime data files for C / Python. Following the handoff's own
instruction, one detritus file was conditional and the other 196 were default installation
targets.

Result of this integration:

- 196 default files: applied to C and Python byte-identically
- pure-water `a_w, b_w`: no numerical change, because the supplied data equal the existing
  canonical data down to the SHA-256
- `psi_T`: complete 330–1100 nm table applied
- default phytoplankton `a_ph*`: complete 330–1100 nm table applied; the 800–1100 nm rows are
  exactly 0
- AHN mineral TSM: 8 scalar tables + 4 vector Mie tables applied
- atmospheric aerosol Mie: 176 tables applied
- gas absorption, 6 gases: 330–1100 nm, 40 levels, 1 nm tables applied
- EAP, 17 species: production scattering stays disabled
- Raman: excluded
- external reference model: not modified

## 3. Important limitation — Chl-linked detritus

The 330–1100 nm detritus file of the handoff was a **conditional candidate**. The default
installation tool was designed to exclude it.

A matched-input RT test with the candidate inserted directly into the current OCRT `L = 200`
coefficient path gave:

- negative reconstructed P11: minimum about `−8.070e-02`
- `rt_solver_sos_pol` failure even in the reduced SOS
- no valid output

The candidate was therefore not adopted for production.

Final policy:

- active detritus: the validated existing 350–850 nm data
- 330–1100 candidate: kept under `water_iop/candidates/`
- a Chl > 0 constituent run that requests a wavelength outside 350–850 nm fails loudly in both C
  and Python
- the Chl absorption table itself covers 330–1100 nm, but because of the Chl-linked scattering
  the full Chl constituent run range is currently 350–850 nm

This does not hide a failure; it applies the conditional-adoption rule of the handoff as written.

## 4. Official range interpretation

330–1100 nm is not one global number. It is a **component-gated contract** that holds when every
selected component satisfies its data range.

Usable over the full range:

- pure water
- CDOM
- the 4 AHN mineral TSM species
- atmospheric Rayleigh / aerosol
- six-gas absorption
- black / flat / black-Fresnel / coupled-ocean RT, provided the selected water constituents
  satisfy the range

Limits:

- Chl > 0 OCRT constituent: 350–850 nm, because of the active detritus
- EAP species scattering: disabled; native range 350–850 nm
- Raman: excluded
- legacy CCRR: non-default, outside this range contract

## 5. Reading order

1. `docs/DATA_INTEGRATION_MANIFEST_KO.md`
2. `docs/SPECTRAL_DATA_METHODS_AND_PROVENANCE_330_1100_KO.md`
3. `docs/DETRITUS_CANDIDATE_DECISION_KO.md`
4. `docs/VALIDATION_REPORT_KO.md`
5. `docs/HANDOFF_ARCHIVE_RECOVERY_NOTE_KO.md`
6. `validation/INSTALLATION_DECISION_MAP_C_PYTHON.csv`
7. `validation/validation_summary.json`

## 6. Representative validation results

- handoff mapping: 197 rows
- accepted runtime data: 196 files
- C / Python accepted data exact equality: PASS
- full-range active Mie: 180/180 pass range, finite-value and basic Mueller-bound checks
- C / Python scalar IOP maximum relative difference: `5.49e-15`
- C native coupled LUT at 330 / 1100 nm: physical outputs identical to the authoritative replay;
  only diagnostic round-off up to `1e-13`
- Python native coupled LUT at 330 / 1100 nm: exact equality; pass1 = 1, water = 1, pass2 = 1,
  per-cell solves = 0
- C fixed-IOP performance median: old 0.72 s, new 0.72 s; 0 % regression
- selected C regressions: PASS
- selected Python tests: 25 pass; final quick set 9 pass

## 7. Package integrity

Use the top-level `SHA256SUMS.txt`. The C / Python data-file mapping is defined by
`validation/ACCEPTED_DATA_SHA256SUMS.txt` and `validation/INSTALLATION_DECISION_MAP_C_PYTHON.csv`.

---

# OCRT 330–1100 nm 분광자료 통합본 — 먼저 읽기

작성일: 2026-08-14
대상: OCRT C/Python 쌍 구현
(외부 참조 모델은 2026-07-28 구성을 유지했으며 이번 분광 확장 대상이 아니다.)

## 1. 최종 버전

- C: `OCRT-v1.2-2026-08-14-KST-spectral-data-330-1100-component-gated`
- Python: `pyOCRT-v1.2-2026-08-14-spectral-data-330-1100-component-gated`

## 2. 핵심 결론

인계본에는 C/Python 용 쌍 실행 자료 197개가 들어 있었다. 인계본 자체의 지침대로 쇄설물 1개는 조건부였고,
나머지 196개는 기본 설치 대상이었다.

이번 통합 결과:

- 기본 설치 196개: C/Python 에 바이트 동일하게 적용 완료
- 순수해수 `a_w, b_w`: 기존 표준 자료와 인계 자료가 SHA-256 까지 같으므로 수치 변경 없음
- `psi_T`: 330–1100 nm 완전표 적용
- 기본 식물플랑크톤 `a_ph*`: 330–1100 nm 완전표 적용, 800–1100 nm 행은 정확히 0
- AHN 광물 TSM: 스칼라 8개 + 벡터 Mie 4개 적용
- 대기 에어로졸 Mie 176개: 적용
- 기체 흡수 6종: 330–1100 nm, 40 level, 1 nm 표 적용
- EAP 17종: 생산 산란 비활성 유지
- Raman: 제외
- 외부 참조 모델: 수정하지 않음

## 3. 중요한 제한 — Chl 연동 쇄설물

인계본의 330–1100 nm 쇄설물 파일은 **조건부 후보**였다. 기본 설치 도구도 이를 제외하도록 설계되어 있었다.

현재 OCRT 의 `L = 200` 계수 경로에 후보를 직접 넣어 같은 입력으로 RT 를 시험한 결과:

- 복원 P11 이 음수: 최소 약 `−8.070e-02`
- 축소 SOS 에서도 `rt_solver_sos_pol` 실패
- 정상 산출값 생성 실패

따라서 이 후보는 생산에 채택하지 않았다.

최종 정책:

- 활성 쇄설물: 검증된 기존 350–850 nm 자료 유지
- 330–1100 후보: `water_iop/candidates/` 에 보존
- Chl > 0 성분 실행이 350–850 nm 밖을 요청하면 C/Python 모두 fail-loud
- Chl 흡수표 자체는 330–1100 nm 를 완전히 지원하지만, Chl 연동 산란 때문에 전체 Chl 성분 실행 범위는
  현재 350–850 nm

이는 실패를 숨긴 것이 아니라 인계본의 조건부 채택 규칙을 그대로 적용한 결과다.

## 4. 공식 범위 해석

330–1100 nm 는 전역 숫자 하나가 아니라 **선택한 모든 성분이 자료 범위를 충족할 때 유효한 성분 게이트
계약**이다.

전 범위 사용 가능:

- 순수해수
- CDOM
- AHN 광물 TSM 4종
- 대기 Rayleigh/에어로졸
- 6종 기체 흡수
- 흑색/평면/흑색 Fresnel/결합 해양 RT, 단 선택한 해수 성분이 범위를 만족할 것

제한:

- Chl > 0 OCRT 성분: 활성 쇄설물 때문에 350–850 nm
- EAP 종별 산란: 비활성, 고유 범위 350–850 nm
- Raman: 제외
- 옛 CCRR: 비기본, 이번 범위 계약 밖

## 5. 문서 읽는 순서

1. `docs/DATA_INTEGRATION_MANIFEST_KO.md`
2. `docs/SPECTRAL_DATA_METHODS_AND_PROVENANCE_330_1100_KO.md`
3. `docs/DETRITUS_CANDIDATE_DECISION_KO.md`
4. `docs/VALIDATION_REPORT_KO.md`
5. `docs/HANDOFF_ARCHIVE_RECOVERY_NOTE_KO.md`
6. `validation/INSTALLATION_DECISION_MAP_C_PYTHON.csv`
7. `validation/validation_summary.json`

## 6. 대표 검증 결과

- 인계 매핑: 197행
- 승인된 실행 자료: 196 파일
- C/Python 승인 자료 정확 일치: PASS
- 전 범위 활성 Mie: 180/180 범위·유한값·기본 Mueller 한계 검사 PASS
- C/Python 스칼라 IOP 최대 상대차: `5.49e-15`
- C native 결합 LUT 330/1100 nm: 물리 출력이 정본 재실행과 동일, 진단 반올림만 `1e-13` 까지
- Python native 결합 LUT 330/1100 nm: 정확 일치, pass1=1, water=1, pass2=1, 셀별 해석=0
- C 고정 IOP 성능 중앙값: 이전 0.72 s, 현재 0.72 s, 회귀 0 %
- 선택 C 회귀: PASS
- 선택 Python 테스트: 25 통과, 최종 빠른 세트 9 통과

## 7. 패키지 무결성

최상위 `SHA256SUMS.txt` 를 쓴다. C/Python 자료 파일 매핑은 `validation/ACCEPTED_DATA_SHA256SUMS.txt` 와
`validation/INSTALLATION_DECISION_MAP_C_PYTHON.csv` 가 기준이다.
