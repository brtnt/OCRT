# OCRT 330–1100 nm 분광자료 통합본 — 먼저 읽기

작성일: 2026-08-14  
대상: OCRT C/Python paired implementation  
OSOAA: 현행 2026-07-28 reference package 유지, 이번 분광확장 대상 아님

## 1. 최종 버전

- C: `OCRT-v1.2-2026-08-14-KST-spectral-data-330-1100-component-gated`
- Python: `pyOCRT-v1.2-2026-08-14-spectral-data-330-1100-component-gated`

## 2. 핵심 결론

첨부 handoff에는 C/Python용 paired runtime 자료 197개가 들어 있었다. handoff 자체의 지침대로 detritus 1개는 조건부였으며, 나머지 196개는 기본 설치 대상이었다.

이번 통합 결과는 다음과 같다.

- 기본 설치 196개: C/Python에 byte-identical하게 적용 완료
- 순수해수 `a_w,b_w`: 기존 canonical 자료와 첨부 자료가 SHA-256까지 동일하므로 수치변경 없음
- `psi_T`: 330–1100 nm 완전표 적용
- 기본 phytoplankton `a_ph*`: 330–1100 nm 완전표 적용, 800–1100 nm는 데이터 행 자체가 정확히 0
- AHN 광물 TSM scalar 8개 + vector Mie 4개: 적용
- 대기 aerosol Mie 176개: 적용
- 기체흡수 6종: 330–1100 nm, 40 level, 1 nm 표 적용
- EAP 17종: production scattering 비활성 유지
- Raman: 제외
- OSOAA: 수정하지 않음

## 3. 중요한 제한 — Chl-linked detritus

handoff의 330–1100 nm detritus 파일은 **조건부 candidate**였다. 기본 설치도구도 이를 제외하도록 설계되어 있었다.

현재 OCRT의 `L=200` coefficient path에 candidate를 직접 넣어 matched-input RT를 시험한 결과:

- 복원 P11이 음수: 최소 약 `-8.070e-02`
- 축소 SOS에서도 `rt_solver_sos_pol` 실패
- 정상 산출값 생성 실패

따라서 이 candidate는 production에 채택하지 않았다.

최종 정책:

- active detritus: 검증된 기존 350–850 nm 자료 유지
- 330–1100 candidate: `water_iop/candidates/`에 보존
- Chl>0 constituent run이 350–850 nm 밖을 요청하면 C/Python 모두 fail-loud
- Chl absorption 표 자체는 330–1100 nm를 완전히 지원하지만, Chl-linked scattering 때문에 전체 Chl constituent run 범위는 현재 350–850 nm

이는 실패를 숨긴 것이 아니라 handoff의 조건부 채택 규칙을 그대로 적용한 결과다.

## 4. 공식 범위 해석

330–1100 nm는 전역 숫자 하나가 아니라 **선택한 모든 성분이 자료범위를 충족할 때 유효한 component-gated contract**이다.

전 범위 사용 가능:

- pure water
- CDOM
- AHN mineral TSM 4종
- atmospheric Rayleigh/aerosol
- six-gas absorption
- black/flat/black-Fresnel/ocean coupled RT, 단 선택 water constituent가 범위를 만족할 것

제한:

- Chl>0 OCRT constituent: active detritus 때문에 350–850 nm
- EAP species scattering: 비활성, native 350–850 nm
- Raman: 제외
- legacy CCRR: 비기본, 이번 범위계약 밖

## 5. 문서 읽는 순서

1. `docs/DATA_INTEGRATION_MANIFEST_KO.md`
2. `docs/SPECTRAL_DATA_METHODS_AND_PROVENANCE_330_1100_KO.md`
3. `docs/DETRITUS_CANDIDATE_DECISION_KO.md`
4. `docs/VALIDATION_REPORT_KO.md`
5. `docs/HANDOFF_ARCHIVE_RECOVERY_NOTE_KO.md`
6. `validation/INSTALLATION_DECISION_MAP_C_PYTHON.csv`
7. `validation/validation_summary.json`

## 6. 대표 검증 결과

- handoff mapping: 197 rows
- accepted runtime data: 196 files
- C/Python accepted data exact equality: PASS
- full-range active Mie: 180/180 range·finite·basic Mueller bounds PASS
- C/Python scalar IOP max relative difference: `5.49e-15`
- C native coupled LUT at 330/1100: physical outputs identical to authoritative replay; only diagnostic roundoff up to `1e-13`
- Python native coupled LUT at 330/1100: exact equality; pass1=1, water=1, pass2=1, per-cell solve=0
- C fixed-IOP performance median: old 0.72 s, new 0.72 s; 0% regression
- selected C regressions: PASS
- selected Python tests: 25 pass; final quick set 9 pass

## 7. 패키지 무결성

최상위 `SHA256SUMS.txt`를 사용한다. C/Python data file mapping은 `validation/ACCEPTED_DATA_SHA256SUMS.txt`와 `validation/INSTALLATION_DECISION_MAP_C_PYTHON.csv`가 기준이다.
