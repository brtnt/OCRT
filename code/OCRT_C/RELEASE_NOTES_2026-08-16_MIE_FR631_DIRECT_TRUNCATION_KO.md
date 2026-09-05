# OCRT Mie FR631 direct truncation 릴리스 노트

버전: `OCRT-v1.2-2026-08-16-KST-mie-fr631-direct-truncation`

## 실행·병렬 계약

- 하나의 case/행은 native all-view 경로에서 SOS를 한 번 수행하고 요청된 모든 VZA×RAA를 Fourier 합성으로 출력한다.
- production VZA×RAA cell replay는 금지한다. native all-view를 구성할 수 없으면 fail-loud 한다.
- batch 병렬화 단위는 독립 case/행이다. 각 worker는 private cache와 workspace를 사용하며 다른 행의 mutable 상태·중간결과·lock에 의존하지 않는다.

## Mie phase 처리

- 고정 격자 `FR631-2026-08-15`를 적용한다: 0–0.2°/0.005°, 0.2–1°/0.02°, 1–5°/0.05°, 5–20°/0.1°, 20–180°/0.5°.
- P11/P12/P33는 동일한 bracket와 동일한 θ-linear weight를 사용한다.
- production particle kernel에서는 cubic angle interpolation과 particle L=200 reconstruction을 사용하지 않는다.
- 완성된 vector Fourier kernel과 active Mie model은 worker-private RAM에 cache한다.

## Active data

- aerosol 176종과 EAP 17종은 generator 기반 FR631 자료를 설치했다.
- AHN TSM 4종은 legacy 공통 파장과 0.5° node를 exact-preserving remap하여 FR631으로 설치했다.
- TSM 330 nm는 canonical microphysical generator 부재로 350 nm endpoint hold를 사용한다. 이 제한은 자료 header와 validation 문서에 명시한다.
- source-informed TSM 재생성 candidate는 기존 phase와 불일치하여 active 설치하지 않았다.
- active legacy detritus는 variable-grid θ-linear fallback을 사용한다. rejected 330–1100 nm candidate는 계속 제외한다.

## Forward truncation

- OSOAA-compatible real-angle hydrosol contract: `mu1=0.85`, `mu2=0.92`, `A_TRONCA` threshold 0.1.
- P11의 forward cap을 log10-linear continuation으로 대체하고 P12/P33에는 동일 local P11 ratio를 적용한다.
- residual phase는 다시 L=200으로 전개하지 않고 direct θ-linear LUT로 유지한다.
- transport scattering은 `b_eff=b_raw*(1-A_TRONCA/2)`로 축소하며 backscattering closure는 보존한다.
- `A_TRONCA<=0` 또는 threshold 미만인 경우 exact no-op이다.

## 검증 결과

- C/Python truncation phase parity: double precision 범위.
- TSM `A_TRONCA=0`: raw/truncated phase와 full RT가 exact no-op.
- active C/Python Mie 자료 198개 byte-identical; 이 중 FR631 197개가 physical/grid validator를 통과한다.
- Python 전체 회귀시험: 60 passed, exact VZA=0 관련 알려진 경고 1건.
- pressure=0 native all-view: water solve 1회, cell replay 0, diagnostic replay와 동일 출력.
- OSOAA reference executable smoke와 기존 validation-summary integrity를 별도로 확인한다.

## 제한 사항

- EAP 종별 phytoplankton scattering은 constituent production closure에서 계속 비활성이다. FR631/direct-kernel 구현은 준비됐지만 종별 closure의 과학적 활성화는 별도 승인 항목이다.
- historical 2026-07-26 TSM OSOAA campaign은 현재 AHN TSM scalar IOP/phase와 matched input이 아니므로 현재 AHN 출력의 절대 baseline으로 사용하지 않는다.
- exact VZA=0 polarized output extraction singularity는 남아 있다. 진단에는 VZA=0.001°를 사용한다.
