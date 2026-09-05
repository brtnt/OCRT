# OCRT C/Python 330–1100 nm 분광자료 통합 검증보고서

## 1. 요약 판정

| Gate | 결과 |
|---|---|
| Handoff paired payload | 197/197 C-Python byte-identical |
| Accepted runtime installation | 196/196 handoff와 일치, C-Python 동일 |
| Conditional detritus | production 불합격, candidate archive 보존 |
| Pure-water table | 첨부/C/Python SHA 완전 동일 |
| Full-range active Mie | 180/180 coverage·finite·basic Mueller bounds PASS |
| Gas tables | 6/6, 40×771, nonnegative finite PASS |
| C/Python scalar IOP | max rel `5.49e-15` |
| C endpoint native LUT | physical outputs exact; diagnostic roundoff <=`1e-13` |
| Python endpoint native LUT | exact equality |
| Native LUT solve-once | C/Python PASS |
| Performance regression | C median 0.72 s -> 0.72 s, 0% |
| Selected C regressions | PASS |
| Selected Python tests | 25 pass; final quick set 9 pass |
| Full long regression suite | 이번 통합과정에서 전체 재완주하지 않음 |

## 2. Data Gate

### 2.1 Handoff mapping

- mapping rows: 197
- accepted: 196
- conditional/rejected: 1
- accepted handoff-C/handoff-Python/installed-C/installed-Python exact equality: PASS
- detritus candidate archive exact equality: PASS

### 2.2 Pure water

세 파일의 SHA-256이 같다.

```text
uploaded water_coef_z09_1nm.txt
C inputs/water_iop/water_coef_z09_1nm.txt
Python data/water_iop/water_coef_z09_1nm.txt

f173e1a4514b653973def109a5415678afbe4f39564a7b19f94a71ef55ea0e46
```

### 2.3 Active Mie

공식 full-range set:

```text
176 atmospheric aerosol + 4 AHN mineral = 180
```

결과:

- range 330–1100: 180/180
- finite: 180/180
- minimum stored P11: `0.00405305235`
- max `|P12|/P11`: `0.9552198543`
- max `|P33|/P11`: `1.0`
- bounds failure: 0

단순 0.5° trapezoidal normalization error의 최대값은 약 0.03076이었다. 이는 매우 강한 forward peak를 coarse output angle로 다시 적분한 진단값이다. 각 파일의 생성/loader normalization과 별개로 기록하며, 이를 자동 rescale에 사용하지 않았다.

### 2.4 Gas

6개 파일 모두:

- shape: 41×771 (wavelength row + 40 layers)
- wavelength: 330–1100 nm
- finite/nonnegative: PASS

### 2.5 PLOPS와 psi_T

- PLOPS: 771 rows, 330–1100, >=800 exact zero
- psi_T: 771 rows, 330–1100

## 3. Loader/Range Gate

### 3.1 Global range

- 330.0 nm: in-range
- 1100.0 nm: in-range
- 329.999 nm: fail-loud
- 1100.001 nm: fail-loud

### 3.2 Table query

- pure water, psi_T, PLOPS, AHN, aerosol, gas: range check 후에만 보간
- silent endpoint clamp를 범위확장으로 사용하지 않음
- exact endpoint는 extrapolation flag가 아님

### 3.3 Detritus

- active range 350–850
- C/Python 330/1100 Chl request: 명확한 오류
- candidate는 active filename으로 사용하지 않음

## 4. C/Python IOP Parity

파장:

```text
330, 340, 349, 350, 443, 555, 750, 800, 865, 940, 1100 nm
```

비교량:

- `a_w(T=15)`, `b_w`, `bb_w`
- `a_phyto`
- active-range 내 `b_det`, `bb_det`
- red-clay `a_min`, `b_min`, `bb_min`

결과:

```text
max absolute difference = 2.22e-16
max relative difference = 5.49e-15
```

350–850 밖 detritus 값은 양쪽 모두 unavailable로 처리하였다.

## 5. RT Endpoint Gate

### 5.1 C native coupled LUT

조건:

- red-clay TSM=1
- Chl=0
- VZA 0/30/60
- RAA 0/90/180/270
- wavelength 330 and 1100

Native vs authoritative cell replay:

| wavelength | rows | max physical difference | max all-column difference |
|---:|---:|---:|---:|
| 330 | 12 | 0 | `1.00e-13` (`T_total_up_view`) |
| 1100 | 12 | 0 | `1.90e-14` (`T_diff_up_view`) |

TOA I/Q/U, Rrs I/Q/U, rrs I/Q/U, Ed/Eu/Lu, IOP/Kd는 동일했다.

Native diagnostic:

```text
water_cold=1
water_views=VZA count
cells=VZA*RAA
```

### 5.2 Python native coupled LUT

330/1100 nm, C50 aerosol + red-clay TSM, 2×2 geometry:

- native rows == legacy rows exact
- atmosphere pass1=1
- water solve=1
- atmosphere pass2=1
- cell solver calls=0

### 5.3 Aerosol endpoint smoke

C50 AOD865=0.2 + red-clay TSM=1 coupled run:

- 330 nm finite TOA/Rrs/rrs I/Q/U
- 1100 nm finite TOA/Rrs/rrs I/Q/U
- explicit band AOD computed at both endpoints

## 6. Detritus Candidate Gate

File-level tests passed, but solver-level test failed.

443 nm reduced test:

```text
L=200 reconstructed P11 minimum: -8.070e-02
rt_solver_sos_pol FAIL m=0 rc=-2
rt_water_rt_sos_pure FAIL rc=-3
```

결론: reject for Stage-2 production.

Validated legacy active phase at 443 nm converged and generated finite results. PLOPS absorption changes mean the complete new 443-nm output is not expected to be byte-identical to the 2026-08-09 model; detritus phase itself was restored.

## 7. Regression Gate

### 7.1 C passed

- release build
- pure-water Z09 table
- 330–1100 spectral IOP test
- phase wavelength PCHIP
- Stage-2 EAP disabled contract
- DTPSIGN diagnostic
- RAA/direct-glint
- Stage-2 water RAA
- immutable cache/native LUT solve-once
- detritus candidate gate

### 7.2 Python passed

Selected suite:

```text
25 passed
```

Final targeted suite after last cleanup:

```text
9 passed
```

Coverage includes spectral contract, data hygiene, aerosol object, native coupled LUT, Rrs/rrs IQU, water internal-reflection sign, DTPSIGN, pure-water and EAP-disabled contract.

### 7.3 Not claimed

The entire historical long-running regression suite was not completed again in this integration pass. In particular, very long direct value-kernel/m-sign campaigns were not repeated end-to-end. Previously frozen fixes were protected by selected source/runtime regressions and unchanged code paths.

## 8. Performance Gate

Fixed-IOP 3 VZA×40 RAA, single-thread, five interleaved runs:

| version | times [s] | median [s] |
|---|---|---:|
| baseline | 0.75, 0.72, 0.72, 0.72, 0.73 | 0.72 |
| integrated | 0.72, 0.72, 0.72, 0.73, 0.72 | 0.72 |

Median change: 0%.  
CSV output: byte-identical.

## 9. Final Acceptance

### Accepted

- pure-water existing canonical data
- psi_T 330–1100
- PLOPS default absorption 330–1100
- AHN scalar/vector
- 176 aerosol Mie
- six-gas tables
- strict range/data contract
- C/Python paired code integration

### Rejected/limited

- 330–1100 detritus candidate: rejected, archived
- Chl>0 outside 350–850: unavailable/fail-loud
- EAP production scattering: disabled
- Raman and O4: excluded
- OSOAA spectral extension: not performed
