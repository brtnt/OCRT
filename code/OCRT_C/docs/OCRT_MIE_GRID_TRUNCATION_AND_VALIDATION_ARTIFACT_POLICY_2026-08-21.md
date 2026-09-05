# OCRT Mie 각도격자·truncation·검증 산출물 정책

**문서 버전:** 2026-08-21 Rev.1  
**대상:** OCRT C/Python paired implementation, OSOAA external validation, release packaging  
**결정 상태:** APPROVED

## 1. 확정된 물리 정책

### 1.1 Exact FR631

- exact FR631 631점은 canonical angle representation이다.
- production 및 과학 기준해는 **raw FR631, broad truncation OFF**이다.
- 이 결론은 broad truncation의 영향이 작아서가 아니다. 대조실험에서 broad truncation은 입자산란의 큰 비율을 변환하고 directional `Rrs/rrs`를 수 % 이상 변화시킬 수 있었다.
- `0–0.005°` local cap은 시험한 수중 Mie에서 raw 결과와 machine precision 수준으로 동일했으나, 별도 production option으로 추가하지 않는다.

### 1.2 Non-FR631

`non-FR631 = coarse`로 자동 판정하지 않는다.

- **recognized legacy 361×0.5° grid:** 전방피크가 저해상도일 수 있음을 강하게 경고한다. 우선 조치는 FR631 재생성이다. 불가능할 때만 `--ocrt-mie-truncation`을 advanced fallback으로 평가한다.
- **unknown/custom non-FR631 grid:** point count만으로 해상도를 단정하지 않는다. 최소 양의 각도, 전방간격, phase 적분 및 RT 수렴성을 검사하기 전까지 `unvalidated`로 분류한다.
- reader는 truncation을 자동 활성화하지 않는다.

### 1.3 Broad truncation advanced mode

`--ocrt-mie-truncation`은 `OCRT_ADVANCED=1`이 필요한 advanced physics option이다. Canonical FR631 default는 OFF이다.

Broad compatibility test를 수행할 때에는 OCRT와 reference model에서 다음을 모두 맞춘다.

- truncation algorithm 및 `mu1/mu2/threshold`
- coefficient `A`, forward fraction `f=A/2`
- transformed `b`, `omega`, `tau`
- `P11/P12/P33` 및 vector moments
- source-function correction
- finite-column physical depth

다음 비교는 무효이다.

- OCRT raw / OSOAA broad
- OCRT broad / OSOAA raw
- OCRT broad auto-depth / OSOAA broad fixed-depth

## 2. 옵션 분류

| 옵션 | 분류 | 정책 |
|---|---|---|
| `--ocrt-mie-truncation` | advanced | `OCRT_ADVANCED=1` 필수; default OFF |
| `--n-mu-water` | advanced | 수렴성 제어; `OCRT_ADVANCED=1` 필수 |
| 신규 테스트 전용 제어 | debug-only | 반드시 `--debug-*`; `OCRT_DEBUG=1` 강제 |

비-debug 신규 옵션 추가 또는 기본값 변경은 사전 합의 없이 수행하지 않는다.

## 3. 런타임 경고 정책

### 3.1 Exact FR631

- grid warning을 출력하지 않는다.
- raw/OFF가 canonical이라는 사실은 help와 구현 주석에 기록한다.

### 3.2 Legacy 361

한 파일·worker당 한 번 다음 내용을 경고한다.

- forward peak under-resolution 가능성
- FR631 replacement 우선
- broad truncation은 advanced fallback이며 물리량을 변경함
- 자동 활성화하지 않음

### 3.3 Unknown non-canonical

한 파일·worker당 한 번 다음을 출력한다.

- angle count
- minimum positive angle
- 20° 이하 최대 인접간격
- forward adequacy가 검증되지 않았음
- non-FR631 전체를 coarse로 단정하지 말 것

## 4. 기술문서–코드주석 동기화 정책

기술문서 전체를 소스에 복사하지 않는다. 대신 모든 물리·수치 결정은 구현 위치의 코드주석에도 반드시 요약하고 `DOC-REF`로 상세문서를 연결한다.

한 물리 patch는 다음 네 항목을 동시에 갱신해야 완료로 인정한다.

1. implementation
2. implementation-site comment
3. technical document
4. regression/acceptance test

코드주석에는 다음을 포함한다.

- 수식 또는 transform 정의
- 기본값과 활성조건
- 적용범위 및 제한
- 부작용과 reference matching 조건
- 관련 문서명/절

대규모 산포도, case별 표, 실행로그 및 전체 참고문헌은 validation artifact와 기술문서에 보존한다.

## 5. 산포도 및 CSV 보존 정책

### 5.1 CSV가 authoritative source

모든 정량 validation figure는 figure 생성에 사용한 CSV를 OCRT package에 같이 저장한다.

- figure만 저장하는 것을 금지한다.
- `PLOT_DATA_MAP.csv`에서 모든 PNG를 source CSV와 연결한다.
- plot에 적용한 subset/filter/transform을 map 또는 plotting script에 기록한다.
- paired source CSV, derived metrics CSV, figure, plotting script, run configuration 및 hash manifest를 같은 campaign directory에 저장한다.

### 5.2 데이터 업데이트

동일 validation campaign의 자료처리를 다시 수행하여 데이터가 갱신되면 다음 파일을 **한 번에 재생성하고 교체**한다.

1. authoritative paired/source CSV
2. metrics CSV
3. 모든 관련 scatter/residual figure
4. `PLOT_DATA_MAP.csv`
5. campaign metadata
6. SHA-256 manifest

CSV만 또는 그림만 부분 교체하는 것을 금지한다. package installer는 staging directory를 검증한 뒤 campaign directory를 atomic replacement한다.

### 5.3 Package layout

```text
05_VALIDATION/artifacts/
  README_KO.md
  CAMPAIGN_INDEX.csv
  <campaign_id>/
    CAMPAIGN_METADATA.json
    PLOT_DATA_MAP.csv
    *.csv
    *.png
    repro/ or plotting scripts
    SHA256SUMS.txt
    ARTIFACT_MANIFEST.json
```

### 5.4 Release Gate

다음을 통과해야 release validation artifact가 유효하다.

- 모든 PNG가 `PLOT_DATA_MAP.csv`에 존재
- map의 모든 CSV가 실제 존재
- orphan PNG 없음
- SHA-256 일치
- source-data revision과 figure revision 동일
- `audit_validation_artifacts.py` PASS

## 6. P0 TODO

1. CLI option hygiene audit: production / advanced / debug-only 재분류
2. raw FR631 reference convergence freeze
3. broad all-view coupled-LUT state 저장 및 multi-geometry one-run 복구
4. `Rrs-Q` direct-beam/interface reference-physics 정합
5. validation artifact audit를 release Gate에 통합

## 7. 근거 validation campaign

- `2026-08-21_noatm_recheck`
- `2026-08-21_forward_cap_control`
- `2026-08-21_broad_truncation_match`
- `2026-08-21_rayleigh_water_stage1`
