# OCRT C/Python DTPSIGN 통합 및 OSOAA 기준 검증 보고서

작성일: 2026-08-02

## 1. 적용 기준

입력 작업지시서 `WORKORDER_DTPSIGN_2026-08-01_KO.md`의 필수 수정은
대기 BOA Fourier 계수의 `phi+pi` 기저를 수중 diffuse-top의 `phi` 기저로
넘길 때 누락된 `(-1)^m` 변환을 복원하는 것이다. 변환은 diffuse-top source
injector 내부가 아니라, 결합장에서 복사한 hand-off 배열에 한 번만 적용한다.

## 2. 최종 버전

- C: `OCRT-v1.2-2026-08-02-KST-dtpsign-phase-diagnostic`
- Python: `pyOCRT-v1.2-2026-08-02-dtpsign-phase-diagnostic`
- OSOAA: `OSOAA_OCRT_reference_tool_2026-07-28-pure-water-z09` — 소스 변경 없음

## 3. 코드 변경

### 3.1 C

`src/rt_solver.c`에서 `coupled.I/Q/U_inwater_per_m`을 `ocrt_dt_I/Q/U`로 복사한
직후 홀수 Fourier mode 블록만 부호 반전한다.

```c
for (int m = 1; m <= m_max; m += 2)
    for (int j = 0; j < n_mu_water; ++j) {
        const size_t z = (size_t)m * n_mu_water + j;
        dtI[z] = -dtI[z]; dtQ[z] = -dtQ[z]; dtU[z] = -dtU[z];
    }
```

source injector는 수정하지 않았다. 따라서 source/operator closure 항등식이 유지된다.

### 3.2 Python

공통 helper `apply_diffuse_top_basis_sign_inplace()`를 추가하고 다음에 적용했다.

- 단일 coupled solve
- native coupled angular LUT solve
- batch R1 atmosphere-to-water hand-off

NumPy/CuPy 배열 모두 같은 in-place 연산을 사용한다.

### 3.3 OSOAA

DTPSIGN 결함은 OCRT의 대기→수중 hand-off에만 존재하므로 OSOAA 코어 수정은
필요하지 않았다. 현행 OSOAA manifest와 Z09 smoke만 재검증했다.

## 4. 의도된 수치 변화

443 nm, SZA=40°, VZA=30°, RAA=0°, 풍속 3 m/s의 coupled pseudo-pure case:

| 항목 | 수정 전 | 수정 후 | 변화 |
|---|---:|---:|---:|
| `Rrs0plus_I` | 1.693286e-02 | 1.733734e-02 | +2.3887% |

Fourier 분해에서 실질적인 변화는 m=1에 집중됐다.

- Rrs I/Q/U의 m=1 진폭 변화: 약 +22.944%
- m=0과 m=2: roundoff 수준
- m>=3: 절대 크기가 약 1e-15 이하인 수치 잡음 영역

이는 작업지시서의 “홀수 mode, 실제로 m=1 중심” 진단과 일치한다.

## 5. 불변 경로

수정 전·후 SHA-256 비교 결과:

- 대기 단독 single stdout/stderr: byte-identical
- 대기 없는 ocean single stdout/stderr: byte-identical
- 대기 단독 LUT CSV: byte-identical

변경 대상은 atmosphere–ocean coupled 실행의 odd azimuthal modes뿐이다.

## 6. Native LUT 구조와 성능

144-cell coupled LUT 진단:

```text
calls=3 water_cold=1 water_views=6 exact_near_nadir=1 cells=144
```

따라서:

- water cold SOS: 1회
- VZA별 target extraction: 6회
- RAA: Fourier reconstruction

수정 전 3회: 2.89, 2.90, 2.88 s, 중앙값 2.89 s
수정 후 3회: 2.91, 2.89, 2.90 s, 중앙값 2.90 s

중앙값 변화 +0.346%는 측정 잡음 범위다. LUT solve-once 구조나 처리속도 하향은
관측되지 않았다.

## 7. C 검증

통과 항목:

- DTPSIGN placement + source/operator closure
- Stage-2 interface FIX1/FIX2/FIX3
- diffuse-top U closure, max abs 2.776e-17 / 5.551e-17
- water RAA convention
- EAP production-disabled contract 및 fixed-bulk truncation 보존
- pure-water Z09 2250-row table
- coupling clamp n_mu_water=48/64/96
- single-target/full-grid parity at n_mu_water=64
- water internal-reflection odd-m sign
- Rrs/rrs IQU and IOP/Kd full-grid
- public RAA/direct-glint branch
- external bottom mode bounds/shape guard
- immutable-cache reuse and water solve-once
- polarized value-kernel native LUT
- 865-nm Chl absorption extension
- P11/P12/P33 wavelength PCHIP

최종 C binary SHA-256:

```text
ee4c6edd3af53157c902f54a1e86e18dfdb4a03864d38ad479d66769dee8b283
```

## 8. Python 검증

- `compileall`: PASS
- CPU quick/paired suite: `28 passed`
- odd-mode-only sign flip + involution: PASS
- batch mode-axis sign flip: PASS
- single-component TSM diagnostic: PASS
- 기존 native LUT 및 coupling tests: PASS

GPU/CuPy가 없는 환경이므로 작업지시서의 7-significant-digit GPU gate는 실행하지
못했다. 두 개의 장시간 4096-azimuth direct value-kernel 시험은 제외했으며, 나머지
value-kernel/native-LUT 시험은 통과했다.

## 9. OSOAA 검증

- `MANIFEST_OSOAA_CURRENT.sha256`: PASS
- active Z09 smoke: 11/11 PASS, relative error < 1e-9

OSOAA 코어는 변경하지 않았다. 작업지시서가 제공한 full-azimuth OSOAA 정합 수치
(I 0.106%, Q 0.120%, U 0.056%, DoLP 0.110 pp)는 패치 채택 근거로 보존했지만,
해당 대규모 캠페인을 이 환경에서 독립 재실행했다고 주장하지 않는다.

## 10. 추가 코드개선 검토

### 10.1 순수수 전용 경로

보고된 0.3–0.4 percentage-point 차이는 현재 canonical tree에서 재현되지 않았다.
정확한 zero-constituent 경로와 1e-6 pseudo-pure 경로의 Rrs 차이는:

- 443 nm: 약 0.0013%
- 555 nm: 약 0.0021%

따라서 기본 경로는 변경하지 않았다.

### 10.2 단일성분 phase diagnostic

채택했다. `OCRT_DEBUG=1 OCRT_DUMP_IOP=1`에서 TSM-only 같은 단일성분도
`OCRT_PHASE_COMPONENT` 및 `OCRT_PHASE_MIX`를 출력한다. C의 일반 출력은
대표 pure/TSM/Chl/mixed case에서 수정 전과 byte-identical이었다. Python도 두
환경변수가 모두 켜진 경우에만 추가 phase evaluation과 stderr 출력을 수행한다.

### 10.3 constituent truncation

채택하지 않았다. 시험 구현은 기존 fixed-bulk parity를 재현하지 못하고 representative
reflectance를 약 40% 변경했다. 현재 정책을 유지한다.

- constituent `--ocrt-mie-truncation`: fail-loud
- fixed-bulk `--iop-mie-truncation`: 기존 검증된 OSOAA parity route로 유지

## 11. Active golden data

DTPSIGN이 coupled odd modes를 의도적으로 변경하므로 coupling-clamp 48/64/96
active references를 2026-08-02 기준으로 갱신했다. 2026-07-25 데이터는 historical
pre-DTPSIGN anchor로 보존했다. single/full-grid test도 새 reference를 사용한다.
