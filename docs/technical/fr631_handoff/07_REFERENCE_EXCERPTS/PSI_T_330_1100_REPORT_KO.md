# OCRT 순수해수 흡광 온도계수 ψT 분광확장 — 330–1100 nm

## 판정

- 기존 `Rottgers2014-approx` 400–900 nm 근사표를 source-based 330–1100 nm, 1 nm 표로 교체한다.
- 330–440 nm의 모든 정수 파장에 `ψT = 0`을 데이터로 명시한다.
- 440 nm는 명시적 zero anchor이다.
- 442–1100 nm의 WOPP 원 2 nm node는 그대로 보존하고, 홀수 파장만 PCHIP 내부보간한다.
- 330–1100 nm 내부에서 runtime 외삽이나 파장별 예외분기는 필요하지 않다.
- 가시광·근적외의 음수 ψT는 물리적인 온도응답이므로 유지한다.

## 적용식

`a_w(λ,T) = a_w(λ,20 °C) + ψT(λ) × (T - 20 °C)`

이번 자료는 흡광 `a_w`만 수정한다. 순수해수 산란, 후방산란, depolarization 및 굴절률은 변경하지 않는다.

## 출처

- 주 논문: Röttgers, McKee & Utschig (2014), *Optics Express* 22, 25093–25108, DOI `10.1364/OE.22.025093`.
- 수치자료: R. Röttgers/HZG WOPP `purewater_abs_coefficients_v3.dat`, Nov. 2016, 2 nm grid.
- 보존한 upstream Git blob: `6fd39b0ce2a087a68c97dc3950e88376ce9d14b4`.

## 검증

- 행 수: 771개, 330–1100 nm inclusive, 1 nm.
- 330–440 nm exact zero: 111/111.
- 330–440 nm non-zero count: 0.
- 442–1100 nm source node 최대 절대오차: `8.674e-19`.
- 0, 5, 10, 20, 25, 30 °C 전 범위에서 `a_w(T)>0`; 최소값 `4.418355e-03 m^-1`.
- 20 °C에서는 보정항이 정확히 0이므로 기준 `a_w(20 °C)`가 전 파장에서 보존된다.
- C 실제 loader: 771행 로드 PASS, compiler/linker warning 0.
- Python 실제 loader: 771행 로드 PASS.
- 13개 파장 × 6개 수온 = 78개 비교에서 C/Python 최대 절대·상대차 모두 0.
- C/Python 설치용 데이터 파일은 byte-identical.

## 정확한 경계 상태 플래그

현재 C의 generic interpolation helper는 정확히 최솟값·최댓값을 조회할 때도 `extrapolated=1`로 표시한다. 수치는 정확하지만 상태 플래그만 잘못된 것이다. 전체 통합 때 비교식을 strict `<`/`>`로 정리하며, 파장별 물리 예외분기는 넣지 않는다.

## 설치 경로

- C: `inputs/water_iop/psi_T_rottgers2014_OCRT.txt`
- Python: `data/water_iop/psi_T_rottgers2014_OCRT.txt`

본 묶음은 C/Python paired install-ready data patch이다. 현재 full active OCRT package 자체는 아직 교체하지 않았다.
