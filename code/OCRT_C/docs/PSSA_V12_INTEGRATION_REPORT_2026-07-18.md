# OCRT v1.2 PSSA 통합 보고서

## 구현 범위
- `--pssa`: He et al. (2018) 기하에 기반한 **대기 direct-beam 부분 구현**.
- 태양 직달 하향빔과 표면 반사 직달빔의 대기 구면 껍질 경로를 계산.
- SOS 다중산란과 관측방향 투과는 locally plane-parallel 유지.
- **수중 굴절 직달빔의 PSSA는 구현하지 않음.** `--pssa`는 0+에서 종료하고 수중 RT는 평면평행을 유지한다.
- Rayleigh 연직구조는 OCRT v1.2의 US Standard Atmosphere 1962 누적 분자 column 사용.
- 에어로졸은 OCRT 2 km exponential, 흡광가스는 기존 AFGL profile 유지.
- 본 보고서는 PCOART-SA 전체 paper-equivalence를 주장하지 않는다.

## v1.2 리팩터링 대응
원 포팅 패키지의 `rt_solver.c` 표면소스 코드는 v1.2에서 `rt_surface_boundary.c`로 분리되어 있다. 따라서 다음을 해당 모듈에 이식했다.
- 평면 Fresnel 반사 직달빔의 레벨별 `xi_refl` 및 `alpha` 적용.
- Cox-Munk 1차 표면소스의 태양 하향 `xi_dn` 적용.
나머지 TOA sunglint, direct transmittance, LUT grid 및 ocean driver 연결은 `rt_solver.c`에 적용했다.

## 검증 결과
### PSSA off 회귀
6개 대표 atmosphere/surface 조건에서 v1.2 배포 바이너리와 stdout/stderr byte-exact.
기존 기능 스모크:
- water branch 14/14 PASS
- CCRR Chl 11/11 PASS
- OCRT organic Chl 10/10 PASS
- Ahn TSM 6/6 PASS
- US62 atmosphere regression PASS

### US62 PSSA 나디르 보정량, n_layers=400, gas off
| wavelength | SZA 70 | 75 | 80 | 85 |
|---|---:|---:|---:|---:|
| 412 nm | 0.347651% | 0.854337% | 2.556747% | 11.069548% |
| 443 nm | 0.257797% | 0.685831% | 2.260429% | 10.869931% |

원 v1.1 exp-8 anchor보다 작아진 것은 US62가 상층 Rayleigh column을 줄이기 때문이다.

### 층수 민감도, 412 nm, SZA 85
- 40층: 12.570778%
- 100층: 11.501903%
- 200층: 11.199014%
- 400층: 11.069548%

고태양각에서 100~200층 이상 권장. 40층이 수렴값을 과대평가하는 방향은 원 포팅 문서와 동일하다.

### 해양 결합 스모크, 443 nm, SZA 75, wind 3
- TOA rho_I: +0.742382%
- Lu0plus: +1.147995%
- Ed0plus: +1.177866%
- Rrs0plus: -0.029541%
- rrs0minus: -0.006170%

Rrs 상쇄 구조가 0.05% 이내로 유지된다.

### 메모리 안전성
ASan/UBSan atmosphere PSSA 및 ocean PSSA 스모크에서 sanitizer finding 0.

## 검증 범위
논문 및 포팅 패키지와 동일하게 SZA 70~85도. 85도 초과는 미검증 외삽이다.


## 2026-07-18 범위 확정 및 복사량 출력 추가

사용자 결정에 따라 수중 PSSA는 현재와 향후 본 작업단계에서 구현하지 않는다. 코드 주석과 사용자 문서에 다음 경계를 명시했다.

```text
--pssa 적용: TOA → 0+ 대기 direct-beam 경로
수중 경로:    0− 이하 plane-parallel
```

추가 출력:

```text
Ed0plus_direct,  Ed0plus_diffuse
Ed0minus_direct, Ed0minus_diffuse
Eu0minus_direct, Eu0minus_diffuse
```

이들은 이미 계산된 flux로 구성한 angular radiometric direct/diffuse 분해이며 추가 SOS solve를 수행하지 않는다. 상세 정의는 `PSSA_SCOPE_AND_RADIOMETRY_SPLIT_2026-07-18.md`를 참조한다.

### 2026-07-19 대기 PSSA 수치 수정

이전 충실도 감사에서 확인한 다음 항목은 후속 패치에서 수정했다.

- reflected-beam layer-local `beta` phase source
- lower-endpoint attenuation 부호 및 spherical amplitude 직접평가
- reflected downward order-1의 고차 SOS 결합
- gas-only ocean PSSA helper

정량 영향과 검증은 `PSSA_NUMERIC_FIXES_2026-07-19.md`를 참조한다. 수중 PSSA는 의도적으로 미구현이며 CDISORT/AccuRT 수치 field를 직접 재현하지 않았으므로 전체 PCOART-SA paper-equivalence는 여전히 주장하지 않는다.
