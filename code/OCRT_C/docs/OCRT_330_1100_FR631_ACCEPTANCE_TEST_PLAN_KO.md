# OCRT 330–1100 nm FR631 최종 Acceptance Test Plan

## A. 필수 입력

- 최신 paired OCRT C/Python
- PART00–PART10 전체 추출본
- OSOAA reference package
- 기존 350–1100 nm regression baseline

## B. 순서

1. `verify_extracted_bundle.py` PASS
2. integration patch 적용
3. C full release build warning 0
4. dry-run 후 data patch 적용
5. Mie 198 structural audit 재실행
6. C/Python reader smoke
7. 구성요소별 330/340/350/400/443/550/865/1100 nm smoke
8. matched-input C/Python I/Q/U regression
9. atmosphere-only, water-only, coupled regression
10. performance/memory audit

## C. 물리 Gate

- finite values only
- `P11>0`, `|P12|<=P11`, `|P33|<=P11`
- phase normalization ≤2e-9
- g consistency ≤1e-9
- SSA in [0,1]
- `Scatter=Extinction×SSA`
- 20°C pure-water exact baseline

## D. RT Gate

- Rayleigh I/Q/U
- aerosol M50C + Ahmad paper/AccuRT representatives
- pure water Rrs/rrs I/Q/U
- Chl-only
- four AHN minerals
- detritus
- coupled clear/turbid cases
- pass-region SZA<75°, VZA<65°

## E. Performance Gate

- no per-angle file reads
- no all-198 simultaneous resident load
- multi-geometry one-run property preserved
- highest turbidity coupled case not slower than reference solely because of repeated phase I/O

## F. Release Gate

- C/Python installed files byte-identical
- manifest retained
- deprecated aliases not silently used
- EAP default disabled
- Raman and O4 exclusions documented

## G. Mie grid·truncation Gate

- exact FR631 validation 기준은 raw/raw, broad truncation OFF이다.
- recognized legacy 361: warning + FR631 replacement preferred
- unknown non-FR631: forward-grid 및 RT convergence audit 전까지 unvalidated
- grid warning은 truncation을 자동 활성화하지 않음
- raw/raw 또는 broad/broad만 유효; mixed policy = `REFERENCE_MISMATCH`
- broad/broad는 transform, `A`, `b/omega/tau`, vector moments, source correction 및 physical depth 동일
- `--ocrt-mie-truncation`, `--n-mu-water`는 advanced option으로 기록

## H. Figure–CSV Artifact Gate

- 모든 I/Q/U·DoLP 산포도는 authoritative paired/source CSV와 함께 package에 저장
- 모든 PNG는 `PLOT_DATA_MAP.csv`에서 source CSV, subset/filter, transform 및 plotting script와 연결
- 자료 갱신 시 CSV, metrics, figures, map, metadata, SHA manifest를 동시에 교체
- orphan PNG, source-less figure, stale hash = 0
- `05_VALIDATION/tools/audit_validation_artifacts.py` PASS

