# OCRT 속도 개선 구현 순서 — 2026-08-25 갱신

## 1. S1/S2 표면 Fourier 및 보간 패치 — 완료

- air-side Cox–Munk all-m Fourier kernel cache
- true `m_max` build cap
- air–water interpolation sort/dedup map memoization
- T_wa Q-fold, EAP-disabled, FR631 reader와 병합
- baseline–optimized Rrs I/Q/U 산포도 및 byte regression 완료

## 2. Native Windows fair gate와 warm OSOAA 격차 — 다음 작업

- fresh/warm 각각 3회 이상 median
- 동일 PSSA OFF, 48/48 Gauss, 26/640 layers, order 100, m=32
- OCRT/OSOAA I/Q/U 1,296행 산포도와 source CSV 봉인
- OSOAA `SURF_MATR` warm보다 OCRT가 느린 경우 persistent surface-operator cache 또는 전용 T_aw/T_wa batch build 적용

## 3. Output-only exact-view projection

- near-nadir exact row 때문에 발생하는 coupled solver 중복 진입 제거
- integration node는 변경하지 않고 output projection만 추가

## 4. Contracted interface operator cache

- raw composite-grid row 대신 final quadrature-space `T_aw/T_wa` 저장
- direct-index cache, memory audit, Q-fold one-time storage

## 5. Subsystem R/T operator 및 Jacobian 재사용

- upper Rayleigh / aerosol-mixed / lower Rayleigh / interface / ocean 5-block
- dependency-based invalidation
- ocean perturbation 시 atmosphere 재사용, aerosol perturbation 시 ocean/interface 재사용
- NN operational inference 채택 시 실시간 우선순위는 낮추되 truth/Jacobian 생성 경로는 유지
