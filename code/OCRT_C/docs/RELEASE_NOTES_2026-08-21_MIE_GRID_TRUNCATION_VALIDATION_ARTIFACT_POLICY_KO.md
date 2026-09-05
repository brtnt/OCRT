# OCRT Mie grid·truncation·validation artifact 정책 반영

**작성일:** 2026-08-21  
**범위:** C/Python paired implementation, 기술문서, 검증 산출물 package

## 확정 사항

- exact FR631의 canonical scientific mode는 raw phase, broad truncation OFF이다.
- `0–0.005°` local cap은 시험 범위에서 무시 가능한 영향이었으나, broad truncation은 directional `Rrs/rrs`를 수 % 변화시킬 수 있다.
- recognized legacy 361 grid는 FR631 교체를 우선 경고한다. unknown non-FR631은 coarse로 단정하지 않고 unvalidated로 경고한다. 어떤 경우에도 truncation을 자동 활성화하지 않는다.
- `--ocrt-mie-truncation`과 `--n-mu-water`는 advanced option이다. 테스트 전용 신규 제어는 `--debug-*` 및 `OCRT_DEBUG=1`을 사용한다.
- 기술적 결정은 구현 위치 코드주석에 요약하고 `DOC-REF`로 상세문서를 연결한다.
- 모든 quantitative figure는 authoritative source CSV와 함께 package에 저장한다. 재처리 시 CSV, metrics, figure, map, metadata, scripts/configuration 및 SHA manifest를 같은 revision으로 atomic replacement한다.

## 구현 및 검증

- C reader에 exact FR631 / recognized legacy 361 / unknown noncanonical grid 분류와 one-time warning을 적용했다.
- Python reader에 동일 정책을 적용했다.
- CLI에서 `--ocrt-mie-truncation`과 `--n-mu-water`의 `OCRT_ADVANCED=1` gate를 명문화했다. 신규 CLI option은 추가하지 않았다.
- C grid-policy/FR631 reader test 및 Python grid-policy test를 추가했다.
- `05_VALIDATION/artifacts/`에 기존 4개 검증 campaign의 PNG와 exact source CSV, map, metadata, repro 및 SHA manifest를 보존했다.
- `05_VALIDATION/tools/audit_validation_artifacts.py`를 release artifact Gate로 사용한다.

## 상세 문서

- `OCRT_MIE_GRID_TRUNCATION_AND_VALIDATION_ARTIFACT_POLICY_2026-08-21.md`
- `OCRT_330_1100_FR631_INTEGRATION_GUIDE_KO.md`
- `OCRT_330_1100_FR631_SPECTRAL_DATA_TECHNICAL_NOTE_KO.md`
- `OCRT_330_1100_FR631_ACCEPTANCE_TEST_PLAN_KO.md`
