# OCRT PSSA 수정 상태와 고태양천정각 층수 경고

작성일: 2026-07-19  
대상: `OCRT-v1.2-2026-07-19-KST-pssa-numeric-fixes`

## 1. 구현 범위

```text
구면경로 적용:
  - TOA -> 대기층의 하향 태양 직달빔
  - TOA -> 평면 해수면 -> 대기층의 반사 태양 직달빔
  - 위 두 빔이 만드는 대기 단일산란 source

locally plane-parallel:
  - 단일산란 이후의 고차 SOS
  - 관측방향 diffuse transport

의도적 미구현:
  - 0- 이하 수중 굴절 직달빔 PSSA
```

수중 PSSA는 숨겨진 fallback이 아니라 명시적 모델 범위 결정이다.

## 2. 2026-07-19 수정 완료 항목

이전 감사에서 확인된 다음 네 항목은 수정 및 정량검증을 완료했다.

| 항목 | 수정 상태 |
|---|---|
| reflected lower endpoint attenuation 부호 | 완료 |
| gas-only ocean PSSA helper 조기 PP 반환 | 완료 |
| reflected downward order-1의 고차산란 결합 | 완료 |
| reflected-beam layer-local `beta` phase source | 완료 |

상세 수식, ablation 영향 및 성능은 `PSSA_NUMERIC_FIXES_2026-07-19.md`를 참조한다.

## 3. 고 SZA layer 수렴 문제

OCRT는 total optical depth에 균일한 layer grid를 사용하고 물리고도 `z(tau)`를 US Standard Atmosphere 1962에서 역산한다. 고 SZA에서는 상층의 넓은 물리고도 layer 안에서 shell secant가 빠르게 변하므로 40층이 PSSA correction 자체를 충분히 해상하지 못할 수 있다.

412 nm, SZA 85°, US62 Rayleigh 수렴시험:

| layers | 최종 PSSA `rho_I`의 1600층 대비 오차 | PSSA correction magnitude 오차 |
|---:|---:|---:|
| 40 | 1.423% | 14.354% |
| 100 | 0.445% | 4.631% |
| 200 | 0.178% | 1.875% |
| 400 | 0.065% | 0.698% |
| 800 | 0.019% | 0.202% |

이 수치는 특정 조건의 수렴시험이며 모든 파장·광학조건의 보편적 상한은 아니다.

## 4. 경고 정책

`--pssa`가 활성화된 경우 다음 조건에서 stderr 경고를 한 프로세스·severity당 한 번 출력한다.

```text
75 <= SZA < 80 deg and n_layers < 100
80 <= SZA < 84 deg and n_layers < 200
84 <= SZA          and n_layers < 400
```

85° 부근 권고:

```text
총 TOA radiance의 보수적 수렴:    n_layers >= 200
PSSA correction 자체 sub-percent: n_layers 400-800
```

경고는 계산을 자동 변경하지 않는다. 판정은 PSSA setup에서 한 번 수행하며 SOS hot loop에는 연산을 추가하지 않는다.

## 5. 남은 제한

1. 수중 PSSA는 구현하지 않는다.
2. CDISORT/AccuRT 논문 수치 field를 직접 fixture로 재현하지 않았다.
3. SZA 85° 초과는 미검증이다.
4. 향후 물리고도 기반 또는 adaptive PSSA layer grid를 검토할 수 있으나, 현재 자동 layer 재구성은 하지 않는다.

## 6. 검증

```bash
./scripts/smoke_pssa_layer_warning.sh ./build/ocrt
./scripts/smoke_pssa_numeric_fixes.sh ./build/ocrt
```

기하 불변조건:

```bash
# 패키지 검증 절차가 tests/test_pssa_geometry.c를 빌드·실행함
```
