# OCRT 350–1100 nm 분광자료 및 산란모델 상세 검토

작성일: 2026-08-09  
대상 기준:

- OCRT C: `OCRT-v1.2-2026-08-02-KST-dtpsign-phase-diagnostic`
- OCRT Python: `pyOCRT-v1.2-2026-08-02-dtpsign-phase-diagnostic`
- 이번 수정본: `2026-08-09 chl-zero-extension-legacy-mie-cleanup`
- OSOAA: 변경하지 않음

## 1. 먼저 용어와 계산구조

해수의 고유광학특성(IOP)은 이 문서에서 다음 세 값으로 구분한다.

- `a(λ)` — 흡수계수, 단위 m⁻¹. 광자가 매질 안에서 소멸되는 정도.
- `b(λ)` — 총 산란계수, 단위 m⁻¹. 모든 방향으로 산란되는 양의 합.
- `bb(λ)` — 후방산란계수, 단위 m⁻¹. 산란 중 90–180° 후방반구로 가는 부분.
- `bb/b` — 후방산란비. 입자 위상함수의 후방 비중을 나타내는 무차원 값.
- `P11, P12, P33` — 편광 산란 위상행렬의 핵심 성분. P11은 세기 각분포, P12와 P33은 Q/U 편광 결합을 결정한다.

OCRT의 해수 IOP는 순수해수, CDOM, phytoplankton absorption, Chl-linked detritus, mineral TSM을 합산해 만든다. 산란량 `b`, 후방산란 `bb`, 그리고 방향분포인 위상행렬은 서로 다른 자료 또는 closure에서 올 수 있으므로 각각 따로 검토해야 한다.

## 2. ψT가 무엇인가

`ψT(λ)`는 **순수해수 흡수계수의 수온 보정계수**다. 현재 코드의 식은 다음과 같다.

\[
a_w(\lambda,T)
=
a_w(\lambda,20^\circ\mathrm{C})
+
\psi_T(\lambda)\,[T-20^\circ\mathrm{C}].
\]

따라서 ψT의 단위는 `m⁻¹ °C⁻¹`이다. 기존 파일 머리글에는 `1/°C`로만 적혀 있었으나 차원이 불완전하므로 이번 패키지에서 머리글을 바로잡았다. 수치 값은 바꾸지 않았다.

현재 파일:

```text
inputs/water_iop/psi_T_rottgers2014_OCRT.txt       # C
data/water_iop/psi_T_rottgers2014_OCRT.txt         # Python
```

현행 표 범위는 400–900 nm이며, 표 자체도 `Rottgers2014-approx`로 표시된 근사표다.

### 2.1 언제 결과에 영향을 주는가

기본 수온이 20°C이면 `T−20=0`이므로 ψT는 어떤 값을 갖더라도 계산에 전혀 기여하지 않는다. 따라서 지금까지 20°C 기준 검증은 ψT 범위 문제와 무관하다.

수온이 20°C가 아니면 영향을 준다. 현재 보간기는 표 밖에서 끝값을 유지하므로:

- 350–399 nm에서는 400 nm ψT를 사용한다.
- 901–1100 nm에서는 900 nm ψT를 사용한다.

이는 실행은 가능하게 하지만 해당 파장의 측정자료에 기반한 보정이라고 볼 수 없다.

예시:

| λ (nm) | 사용 ψT (m⁻¹ °C⁻¹) | 표 범위 여부 |
|---:|---:|---|
| 350 | 1.0000×10⁻⁴ | 400 nm 끝값 유지 |
| 400 | 1.0000×10⁻⁴ | 표 내부 |
| 550 | 1.000558×10⁻⁴ | 표 내부 |
| 900 | 1.017950×10⁻⁴ | 표 내부 |
| 1100 | 1.017950×10⁻⁴ | 900 nm 끝값 유지 |

이번 작업에서는 검증 가능한 ψT 확장자료가 제공되지 않았으므로 수치 모델은 바꾸지 않았다. 20°C 이외의 350–399 및 901–1100 nm는 여전히 후속 과제다.

## 3. “해저표”가 무엇인가

앞선 감사에서 말한 해저표는 OCRT의 해수면(surface) 자료가 아니라 **OSOAA의 얕은 바다 바닥 반사율 스펙트럼**이다.

OSOAA 파일:

```text
fic/OSOAA_SEABED_REFLECTANCES.txt
```

열은 파장과 네 종류의 바닥 반사율이다.

```text
Wa, rsand, rgreen, rbrown, rred
```

예를 들어 모래 바닥, 녹색 식생성 바닥, 갈색/적색 바닥을 유한수심 모델의 하부 경계로 사용할 때 적용한다. 범위는 350–800 nm다.

이 표는 다음 경우에는 사용되지 않는다.

- 무한수심 또는 충분히 깊은 해양
- black bottom
- 일반적인 open-ocean OCRT–OSOAA 검증
- air–water Fresnel surface 자체

사용자가 이번에는 OCRT 완성에만 집중하고 OSOAA를 1100 nm까지 확장하지 않기로 했으므로 OSOAA 해저표와 코드는 변경하지 않았다.

## 4. 첨부 순수해수 IOP 표의 적용 여부

사용자 제공 `water_coef_z09_1nm.txt`는 다음 계약을 갖는다.

- 범위: 200–2449 nm
- 간격: 1 nm
- 행수: 2250
- `a_w`: Pope & Fry 계열 visible + Kou 계열 NIR 조합
- `b_w`: Zhang et al. 2009, S=38.4 g kg⁻¹, T=20°C
- `b_w`는 총산란이며 OCRT가 내부에서 `bb_w=0.5 b_w`를 계산

첨부 파일과 현재 C/Python 파일의 SHA-256은 모두 다음으로 완전히 같았다.

```text
f173e1a4514b653973def109a5415678afbe4f39564a7b19f94a71ef55ea0e46
```

따라서 현행 OCRT는 이미 첨부 자료 전체 범위를 사용하고 있었으며 실제 수치 교체는 필요하지 않았다. 배포본에는 사용자가 첨부한 파일을 다시 복사해 canonical source를 명확히 했지만 byte diff는 없다.

주요 값:

| λ (nm) | a_w (m⁻¹) | b_w (m⁻¹) | bb_w (m⁻¹) |
|---:|---:|---:|---:|
| 443 | 7.067186×10⁻³ | 4.331984×10⁻³ | 2.165992×10⁻³ |
| 850 | 4.202300 | 2.854538×10⁻⁴ | 1.427269×10⁻⁴ |
| 865 | 4.605200 | 2.656188×10⁻⁴ | 1.328094×10⁻⁴ |
| 1100 | 19.08413 | 9.904724×10⁻⁵ | 4.952362×10⁻⁵ |

850 nm 이후에는 순수해수 흡수가 매우 크므로 깊은 물의 water-leaving signal은 빠르게 약해진다. 그러나 이것만으로 입자 산란 외삽을 임의로 해도 된다는 뜻은 아니다. 고탁도 또는 유한수심에서는 입자 산란의 영향이 남을 수 있다.

## 5. OCRT 산란 경로별 350–1100 nm 지원상태

### 5.1 순수해수 분자산란

- 자료: `water_coef_z09_1nm.txt`
- 범위: 200–2449 nm
- 총산란 `b_w`를 읽고 `bb_w=0.5b_w`
- 350–1100 nm 완전 지원

### 5.2 phytoplankton 자체 산란

현재 Stage-2 production constituent model에서는 phytoplankton 자체의 산란을 의도적으로 비활성화했다.

```text
b_phyto = 0
bb_phyto = 0
```

17종 EAP 파일은 generator/API/향후 검증을 위해 유지되지만 production Chl 산란에는 사용하지 않는다. 따라서 EAP 자료가 850 nm에서 끝나는 문제는 현재 기본 산출물의 직접적인 분광 결손이 아니다.

### 5.3 Chl-linked detritus 산란

현재 Chl>0일 때 실제로 활성인 입자 산란은 `Detritus_Stramski2001.mie`를 사용하는 detritus 경로다.

자료 범위:

- bulk metadata: 350–850 nm
- phase P11/P12/P33: 350, 400, 443, 490, 550, 620, 670, 745, 850 nm
- 비대칭인자 g: 약 0.96673–0.96683

`.mie` 파일의 Extinct/Scatter 열은 1.0 placeholder이며 절대 산란량을 정하지 않는다. 실제 closure는 다음과 같다.

\[
b_{b,det}(\lambda)
=
[1-f_{ph}(Chl)]\,
 b_{bp}^{Huot}(550,Chl)
\frac{r(\lambda)}{r(550)},
\]

\[
b_{det}(\lambda)=\frac{b_{b,det}(\lambda)}{r(\lambda)},
\qquad r(\lambda)=\frac{bb}{b}.
\]

두 식을 결합하면:

\[
b_{det}(\lambda)
=
\frac{[1-f_{ph}(Chl)]\,b_{bp}^{Huot}(550,Chl)}{r(550)}.
\]

즉 **현재 코드에서 detritus 총산란 b는 파장에 무관하게 상수**다. `.mie`의 파장 의존성은 주로 다음 두 항에 들어간다.

1. `bb/b`에 따른 후방산란량 `bb_det`
2. P11/P12/P33에 따른 방향·편광 분포

Chl=1 mg m⁻³ 예시:

| λ (nm) | bb/b | b_det (m⁻¹) | bb_det (m⁻¹) |
|---:|---:|---:|---:|
| 350 | 0.00569867 | 0.37194643 | 0.00211960 |
| 550 | 0.00588164 | 0.37194643 | 0.00218765 |
| 745 | 0.00591564 | 0.37194643 | 0.00220030 |
| 850 | 0.00592017 | 0.37194643 | 0.00220198 |
| 1100, 현행 | 0.00592017 | 0.37194643 | 0.00220198 |

850 nm 초과에서 현행 generic interpolation은 850 nm 끝값을 유지한다. 따라서 851–1100 nm에서 `bb/b`와 P11/P12/P33가 모두 850 nm 값으로 고정된다.

### 5.4 mineral TSM

AHN 네 종류 자료는 350–3750 nm를 포함한다.

- Red clay
- Brown earth
- Yellow clay
- Calcareous sand

따라서 mineral TSM 경로는 350–1100 nm에서 별도 외삽 문제가 없다.

## 6. detritus 장파장 외삽 상세 평가

### 6.1 850 nm까지 관측되는 변화

745→850 nm 구간에서:

- `bb/b` 변화: +0.0764%
- 10–170° 구간 P11 중앙 상대변화: 약 0.123%
- P12 중앙 상대변화: 약 0.178%
- P33 중앙 상대변화: 약 0.136%
- 각 성분의 95 percentile 상대변화는 대체로 0.75% 미만

즉 마지막 두 파장 사이에서는 위상형상이 비교적 완만하다. 이 사실은 끝값 유지가 수치적으로 폭발하지 않는 이유지만, 1100 nm까지 동일하다는 물리적 증거는 아니다.

### 6.2 bb/b만 단순 외삽할 경우

마지막 2–4개 파장을 λ, 1/λ, log λ에 대해 직선 외삽한 민감도 시험에서 1100 nm `bb/b`는 현재 끝값 유지보다 약 +0.12%~+0.43% 높게 나왔다.

이는 scalar backscatter ratio만 보면 차이가 작을 수 있음을 뜻한다. 그러나 편광 RT에는 P11뿐 아니라 P12/P33가 필요하므로 bb/b 하나만 연장해서는 충분하지 않다.

### 6.3 P11/P12/P33를 각각 선형 외삽하면 안 되는 이유

745–850 nm의 각 angle별 값을 독립적으로 1100 nm까지 선형 외삽해 시험했다.

- P11 음수는 발생하지 않았다.
- 그러나 일부 각도에서 `|P33| ≤ P11` 물리성 조건을 위반했다.
- 최대 위반량은 외삽 변수에 따라 약 0.016–0.024였다.

즉 각 Mueller 성분을 독립적으로 선형 외삽하면 겉보기에는 매끄러워도 물리적으로 허용되지 않는 산란행렬을 만들 수 있다. SOS 불안정 또는 잘못된 U 편광으로 이어질 수 있으므로 채택하지 않았다.

### 6.4 권고 순서

1. **최선:** detritus의 입경분포와 복소굴절률을 명시해 350–1100 nm의 vector phase를 microphysical model로 재생성한다.
2. **차선:** 1100 nm까지 검증된 측정/모델 P11/P12/P33 자료를 도입한다.
3. **임시 연구안:** Mueller 성분 자체가 아니라 generalized moments 또는 제한된 parameter space에서 외삽하고 다음 제약을 매 파장에서 강제한다.
   - P11 정규화
   - P11≥0
   - |P12|≤P11
   - |P33|≤P11
   - g와 bb/b의 파장연속성
   - reciprocity/positivity 시험
4. **현재 production:** 850 nm endpoint hold를 유지하되, 이를 `SPECTRAL_SUPPORT_350_1100.json`에 provisional policy로 명시한다.

이번 패치에서는 4번을 유지했다. 즉 detritus의 수치는 바꾸지 않았고, 850–1100 nm가 완료된 source-backed model이라고 표기하지 않았다.

## 7. Chl 흡수 장파장 0 외삽 적용

사용자 결정에 따라 phytoplankton absorption은 850 nm 초과에서 정확히 0으로 정의했다.

변경 내용:

- `phyto_absorption_default.csv`를 1100 nm까지 확장
- 855–1100 nm 모든 행에 0 기록
- 코드에 source cutoff 850 nm를 별도 유지
- 850.000001 nm처럼 첫 zero row보다 짧은 비정수 파장도 즉시 0
- 과거 one-time warning 제거: 이제 예외가 아니라 명시적 model policy

이 변경은 이전 코드도 이미 850 nm 초과에서 0을 반환했기 때문에 수치 결과를 바꾸지 않는다.

## 8. legacy narrow Mie 정리

C 패키지에서 다음 파일을 삭제했다.

```text
inputs/Phytoplankton.mie
inputs/TSM_mineral.mie
inputs/Calcareous_sand_from_OSOAA_phase_443_550.mie
```

이 파일들은 각각 443–670 또는 443–550 nm에 불과하고 current production source에서 참조하지 않았다.

삭제하지 않은 파일:

- `Detritus_Stramski2001.mie`: 현재 active Chl-linked scattering
- EAP 17종: legacy가 아니라 current advanced/generator catalog
- AHN 4종: current mineral TSM model, 350–3750 nm

현재 water-related Mie inventory에는 위 세 범주만 남는다.

## 9. 수정 후 검증

### C

- release build: PASS
- pure-water Z09 2250-row/table exact test: PASS
- Chl 850 anchor 및 865 zero-extension: PASS
- legacy narrow Mie absence: PASS
- DTPSIGN phase diagnostic: PASS
- IOP/Kd full-grid: PASS
- water internal-reflection m-sign: PASS
- LUT cache physical-output parity: PASS
- polarized value-kernel native LUT: PASS

443, 850, 850.5, 851, 865, 1100 nm의 축소 coupled 실행에서 수정 전후 stdout은 byte-identical했다.

### Python

- compileall: PASS
- selected paired tests: 16 passed
- 동일 파장의 aw, bw, phyto a/b/bb, detritus a/b/bb는 수정 전후 정확히 동일

## 10. 남은 작업

1. ψT의 source-backed 350–1100 nm 확장 또는 범위 밖 정책 결정
2. detritus vector phase의 850–1100 nm 물리 모델 선정
3. detritus `b_det`가 파장 불변이 되도록 만든 현재 closure가 의도된 것인지 재검토
4. detritus 확장 후 C–Python parity, Mueller positivity, OSOAA/fixed-bulk 비교, LUT 비회귀시험


## 참고문헌

- Röttgers, R., McKee, D., & Utschig, C. (2014). Temperature and salinity correction coefficients for light absorption by water in the visible to infrared spectral region. *Optics Express*, 22, 25093–25108. DOI: 10.1364/OE.22.025093.
- Huot, Y., Morel, A., Twardowski, M. S., Stramski, D., & Reynolds, R. A. (2008). Particle optical backscattering along a chlorophyll gradient in the upper layer of the eastern South Pacific Ocean. *Biogeosciences*, 5, 495–507.
- Stramski, D., Bricaud, A., & Morel, A. (2001). Modeling the inherent optical properties of the ocean based on the detailed composition of the planktonic community. *Applied Optics*, 40, 2929–2945.
- Doxaran, D. et al. (2007). Near-infrared light scattering by particles in coastal waters.
