# OCRT 330-1100 nm 분광자료 최종 인계 패키지

## 목적

이 패키지는 2026-08-09~2026-08-14 세션에서 준비·검증한 OCRT 분광확장용 교체자료를 최신 본 개발 세션으로 전달하기 위한 것이다. C와 Python의 대응 구현에 같은 자료를 동시에 설치하고, 코드의 파장범위·보간 계약을 정리한 후 330-1100 nm acceptance regression을 수행하는 것을 전제로 한다.

## 동결 범위

- 330-1100 nm elastic coupled ocean-atmosphere OCRT
- EAP 17종 species-specific phytoplankton phase: 비활성 유지
- Raman: 이번 완료범위에서 제외
- legacy CCRR Chl branch: 기본모델이 아니며 이번 패치 대상 아님
- OSOAA: 비교 레퍼런스이며 자료를 수정하지 않음
- C/Python: 대응되는 모든 runtime 파일은 byte-identical하게 유지

## 먼저 읽을 문서

1. `01_DOCUMENTS/OCRT_330_1100_INTEGRATION_GUIDE_KO.docx`
2. `01_DOCUMENTS/OCRT_330_1100_SPECTRAL_DATA_TECHNICAL_NOTE_KO.docx`
3. `01_DOCUMENTS/CODE_CHANGE_CHECKLIST_KO.md`
4. `01_DOCUMENTS/ACCEPTANCE_TEST_PLAN_KO.md`

## 자료 상태의 핵심

- 각 구현의 runtime payload는 총 197개이다: 기본 `a_ph*` 1, AHN scalar 8, AHN mineral Mie 4, detritus Mie 1, aerosol Mie 176, gas xsec 6, `psi_T` 1.
- 기본 설치는 조건부 detritus 1개를 제외한 196개를 C/Python에 동시에 설치한다.
- `Detritus_Stramski2001.mie`는 전 범위 파일이 준비됐지만 legacy phase와 잔차가 남아 **조건부**이다. 본 세션에서 명시적으로 채택하고 matched-input I/Q/U regression을 통과시킨 뒤 활성화해야 한다.
- AHN mineral phase는 Ahn (1990) source-informed seawater forward-Mie 복원이다. exact historical reproduction을 주장하지 않는다.
- NO2 330-349 nm는 Vandaele 220/294 K 상대분광형을 기존 canonical 350 nm 층별 값에 정합한 확장이다. O2-O2/O4 CIA는 미포함이다.
- PLOPS 기본 흡광은 단일 Diatom 스펙트럼이 아니라 34개 배양자료·6개 분류군을 class-balanced robust median으로 집계한 `a_ph*`이다.

## 검증과 설치

```bash
python 06_TOOLS/verify_handoff_payload.py
python 06_TOOLS/apply_data_patch.py --c-root <OCRT_C_ROOT> --python-root <OCRT_PY_ROOT>
python 06_TOOLS/apply_data_patch.py --c-root <OCRT_C_ROOT> --python-root <OCRT_PY_ROOT> --apply
```

기본 설치도구는 detritus를 건너뛴다. 조건부 detritus까지 설치할 때만 다음 옵션을 사용한다.

```bash
python 06_TOOLS/apply_data_patch.py --c-root <OCRT_C_ROOT> --python-root <OCRT_PY_ROOT> --apply --include-detritus
```

설치도구는 자료 파일만 복사한다. 필수 코드수정과 acceptance regression은 통합가이드에 따라 별도로 수행한다.

## 패키지 구성

- `02_REPLACEMENT_DATA/`: C/Python 목적지 구조를 그대로 반영한 paired payload
- `03_SOURCE_AND_PROVENANCE/`: 원자료, 생성 스크립트, 보고서, 그래프, 검증자료. PLOPS와 AHN scalar는 중복 대용량 archive 대신 완전한 펼침형 source/build tree를 여기에 수록한다.
- `04_REPRODUCIBILITY_ARCHIVES/`: AHN mineral Mie, detritus, aerosol, gas, `psi_T`의 원 재현 ZIP. PLOPS와 AHN scalar는 `03_SOURCE_AND_PROVENANCE/`의 펼침형 tree가 재현 기준이다.
- `05_VALIDATION_AND_MANIFESTS/`: 설치맵, 상태표, 최종 지원계약, SHA-256
- `06_TOOLS/`: 무결성 검사, dry-run/backup 설치, 코드 integration point 검색
- `07_REFERENCES/`: 참고문헌과 source-file register

## 현재 상태

`INSTALL_READY_NOT_YET_APPLIED_TO_LATEST_MAIN_OCRT`
