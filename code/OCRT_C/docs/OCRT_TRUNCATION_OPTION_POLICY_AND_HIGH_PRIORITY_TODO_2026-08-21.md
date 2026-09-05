# OCRT truncation·고급 옵션 정책 및 최우선 TODO

**결정일:** 2026-08-21

## 확정 정책

1. `--ocrt-mie-truncation`과 `--n-mu-water`는 **advanced option**이다.
2. 테스트 전용 신규 제어는 반드시 `--debug-*`로 추가하고 `OCRT_DEBUG=1`을 강제한다.
3. 비-debug 옵션의 신규 추가·기본값 변경은 사전 합의 없이 수행하지 않는다.
4. OCRT–OSOAA 비교에서 truncation 정책이 다르면 그 비교는 무효이다.
5. 모든 비교는 명시적으로 다음 중 하나로 고정한다.
   - `raw/raw`: OCRT OFF, OSOAA OFF
   - `broad/broad`: OCRT ON, OSOAA ON, 동일 변환·동일 물리수심

## Truncation 최종 결정

### Production 및 과학 기준해

- **Raw FR631**
- **Broad truncation OFF**
- Canonical truncation coefficient: `A = 0`

근거:

- 최신 FR631은 631개 명시 각도와 771개 phase-wavelength column을 제공한다.
- raw/raw와 broad/broad-fixed의 정합성 차이는 작다.
  - `Rrs-I` case-band 평균 MAPE: raw/raw 1.109%, broad/broad-fixed 1.011%
  - `rrs-I` case-band 평균 MAPE: raw/raw 0.329%, broad/broad-fixed 0.314%
- broad는 일부 조건에서 입자산란의 최대 49.84%를 변환한다.
- broad는 `Rrs-Q` 비정합을 해결하지 않는다.
- 따라서 작은 정합성 개선만으로 broad를 물리 기본값으로 채택할 근거가 없다.

### OSOAA 호환·가속용 advanced mode

Broad mode가 필요할 때에는 임의 튜닝하지 않고 다음 OSOAA 규격을 고정한다.

- `mu1 = 0.85`
- `mu2 = 0.92`
- activation threshold `A = 0.10`
- forward fraction `f = A/2`
- `b_eff = b_raw * (1 - A/2)`
- 대응 `omega`, `tau`, vector moments 및 source correction을 양 코드에서 동일 적용
- finite-column 비교에서는 **변환 전 physical depth를 고정**
- auto-depth를 사용할 경우 OCRT와 OSOAA에서 transformed-domain depth/optical depth를 동일화

`0.005°` local cap은 raw 결과와 machine precision 수준으로 동일하므로 별도 production 옵션으로 만들지 않는다.

## 최우선 TODO

### P0-1. Option hygiene audit

다음 기존 옵션을 분류한다.

- advanced physics/numerics
- debug-only test control
- deprecated/remove candidate

대상:

- `--ocrt-mie-truncation`
- `--n-mu-water`
- `--ocrt-mie-moment-mode`
- `--ocrt-mie-moment-nmu`
- `--ocrt-mie-ss-mode`
- `--iop-mie-*`
- `--water-shared-grid`
- `--water-view-as-node-off`

### P0-2. Truncation contract freeze

- `raw/raw`를 공식 validation 기준으로 고정
- `broad/broad`를 advanced compatibility regression으로 유지
- 실행 로그와 산포도에 truncation mode, `A`, `f`, depth policy를 의무 기록
- raw-vs-broad 비교를 최종 Gate에 사용하지 않음

### P0-3. Broad all-view coupled-LUT 구조 수정

- broad constituent phase의 native all-view state 저장
- per-cell replay 금지
- multi-geometry one-run 유지
- 수정 후 Rayleigh 전체 VZA×RAA I/Q/U 산포도 재실행

### P0-4. Raw reference convergence freeze

- OCRT `n_mu_water` convergence
- OSOAA radiance/Mie quadrature convergence
- SOS order convergence
- 4 AHN × 443/555/865 nm
- `Rrs/rrs IQU` 산포도
- 수렴한 raw/raw 결과를 immutable reference로 동결

### P0-5. Rrs-Q reference physics

- Fresnel-transmitted direct-beam Q
- interface Mueller rotation
- Q/U reference plane
- OCRT–OSOAA 동일화 후 재검증
