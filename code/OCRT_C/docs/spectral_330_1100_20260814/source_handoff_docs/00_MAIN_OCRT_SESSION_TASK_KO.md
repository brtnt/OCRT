# 본 OCRT 개발 세션 전달용 작업지시

이 ZIP은 OCRT의 공식 지원범위를 330-1100 nm로 확장하기 위해 이 세션에서 확정한 교체자료, 원자료·재현자산, 통합가이드, 코드점검표, acceptance test plan을 포함한다.

## 본 세션에서 수행할 작업

1. `00_README_FIRST_KO.md`와 `01_DOCUMENTS/OCRT_330_1100_INTEGRATION_GUIDE_KO.docx`를 먼저 읽는다.
2. `python 06_TOOLS/verify_handoff_payload.py`를 실행해 197개 mapping row와 paired SHA-256 검사를 통과시킨다.
3. 최신 OCRT C와 최신 paired OCRT Python을 같은 integration branch에서 고정한다.
4. `06_TOOLS/apply_data_patch.py`를 dry-run한 뒤 기본 196개 파일을 C/Python에 동시에 설치한다.
5. `Detritus_Stramski2001.mie`는 조건부 자료이다. 본 개발 세션에서 source-recipe candidate를 새 canonical로 명시적으로 채택하고 I/Q/U regression을 수행할 때만 `--include-detritus`로 설치한다.
6. `01_DOCUMENTS/CODE_CHANGE_CHECKLIST_KO.md`에 따라 330 nm range gate, exact endpoint flag, silent clamp, table capacity, gas reader, EAP 비활성 상태를 점검·수정한다.
7. `01_DOCUMENTS/OCRT_330_1100_ACCEPTANCE_TEST_PLAN_KO.md`에 따라 data/loader, IOP, phase physicality, C/Python matched-input RT, coupled spectral sweep를 수행한다.
8. EAP는 비활성 상태를 유지하고 Raman 및 O4/CIA를 이번 완료범위에 포함하지 않는다.
9. 결과물로 최신 OCRT C 패키지, paired Python 패키지, 변경내역, C/Python regression 결과, 330-1100 nm acceptance 보고서를 반환한다.

## 금지사항

- 한 구현에만 자료를 적용하지 않는다.
- 파장별 과학값을 코드의 component-specific `if` 문으로 대체하지 않는다.
- 범위 밖 값을 silent endpoint clamp로 정상값처럼 처리하지 않는다.
- AHN scalar·vector 또는 C/Python의 신·구 자료를 혼합한 partial integration을 남기지 않는다.
- EAP를 우회 활성화하지 않는다.

## 공식 범위

```text
330-1100 nm
elastic atmosphere-ocean coupled OCRT
C/Python paired implementation
EAP disabled
Raman excluded
legacy CCRR non-default
O2-O2/O4 CIA excluded
```
