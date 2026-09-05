# OCRT 의사구면근사(PSSA) 구현 상세

최초 작성: 2026-07-18  
최종 갱신: 2026-07-19  
실행 옵션: `--pssa`

## 1. 범위

현재 OCRT는 He et al. (2018)의 pseudo-spherical direct-beam 처리 중 다음 대기 항을 구현한다.

- 하향 태양 직달빔 TOA→각 대기 level
- 평면 해수면 반사 직달빔 TOA→surface→각 대기 level
- 위 직달빔들이 만드는 single-scattering source
- single scattering 이후 locally plane-parallel SOS

수중 굴절 직달빔 PSSA는 사용자 결정에 따라 구현하지 않는다. `--pssa`는 0+에서 종료하고 0− 이하 수중 RT는 평면평행이다. 따라서 전체 PCOART-SA paper-equivalence는 주장하지 않는다.

## 2. 하향 직달빔

목표 level `k`의 고도를 `z_k`, 지구반지름을 `R`, surface solar zenith를 `theta0`라 하면

\[
D_k=(R+z_k)\sin\theta_0
\]

\[
F(z)=\sqrt{(R+z)^2-D_k^2}
\]

\[
\xi_{dn}(k)=\sum_{i<k}\Delta\tau_i
\frac{F(z_i)-F(z_{i+1})}{z_i-z_{i+1}}
\]

이다. 큰 수의 뺄셈을 피하기 위해 코드에서는 다음 항등식을 사용한다.

\[
F(z)^2=((R+z)\mu_0)^2+(z-z_k)(2R+z+z_k)\sin^2\theta_0
\]

\[
\sec_i(k)=\frac{2R+z_i+z_{i+1}}{F_i+F_{i+1}}
\]

`rt_pssa_apply()`는 최종 Rayleigh+aerosol+gas `h[]`가 완성된 뒤 호출되고,

\[
ch[k]=\frac12\exp[-\xi_{dn}(k)]
\]

로 direct source attenuation을 설정한다.

## 3. 반사 직달빔

He et al. Fig. 1(b), Eq. (8):

\[
\frac{\sin\alpha_k}{\sin\beta_k}=\frac{R+z_k}{R},
\qquad 2\alpha_k=\beta_k+\theta_0
\]

`alpha_k`는 이분법으로 풀고,

\[
\beta_k=2\alpha_k-\theta_0
\]

를 저장한다. `alpha_k`는 surface Fresnel Mueller 계수에, `beta_k`는 reflected-beam phase source의 입사방향에 사용한다.

왕복 경사광학두께는

\[
\xi_{refl}(k)=
\sum_{i<k}\Delta\tau_i\sec_i(\alpha_k)
+2\sum_{i\ge k}\Delta\tau_i\sec_i(\alpha_k)
\]

이며 amplitude를 직접

\[
B_{refl}(k)=\frac12\exp[-\xi_{refl}(k)]
\]

로 평가한다. 이는 plane-parallel baseline×correction 계산에서 생길 수 있는 underflow/overflow를 피한다.

## 4. Layer-local beta phase source

`rt_angular_basis_eval_scalar()`와 `rt_angular_basis_eval_spin2()`가 임의 signed direction cosine에서 기존 quadrature-table builder와 같은 generalized angular basis를 계산한다.

PSSA flat reflected-beam source에서는 direction-reversal parity를 유지하기 위해 입사 basis를 `mu=-cos(beta_k)`에서 계산하고 기존 reflected-source output hemisphere mapping을 그대로 사용한다.

- aerosol active: P11/P12/P33의 beta/gamma/alpha/zeta moments를 endpoint별 kernel로 사전 contraction
- pure Rayleigh: l=2 basis만 사전계산
- SOS order loop: 준비된 source field와 기존 operator만 사용

phase/moment/kernel은 산란차수마다 다시 계산하지 않는다.

## 5. First-order up/down coupling

반사 직달빔 single scattering은 upward와 downward field를 모두 생성한다. 두 방향을 `order1_*` 및 누적 field에 저장하여 order 2 이상 source에 모두 참여시킨다.

이 처리에 추가 SOS solve나 loop가 생기지 않는다. 기존 source operator가 이미 양 방향 배열을 순회한다.

## 6. Gas-only ocean helper

해양 결합용 `T_dir_dn` helper는 Rayleigh와 aerosol이 모두 0이어도 US62 altitude grid를 만들고 AFGL gas absorption을 layer별로 적분한다. 그 뒤 동일한 `rt_pssa_apply()`를 사용한다.

```text
scattering optical depth = 0
```

은

```text
direct-beam extinction = 0
```

을 의미하지 않는다.

## 7. 자료구조

`rt_atm_t`:

```text
pssa_active
pssa_xi_dn[n_layers+1]
pssa_xi_refl[n_layers+1]
pssa_alpha[n_layers+1]
pssa_beta[n_layers+1]
```

모두 atmosphere object가 소유하며 `rt_atm_free()`에서 해제한다.

## 8. 성능

- shell geometry: case/wavelength별 setup 1회
- beta angular basis: Fourier mode별 setup 1회
- SOS order별 phase 재생성: 없음
- order별 동적할당 증가: 없음
- pure Rayleigh: l=2 전용 fast path

교대 benchmark 결과는 `PSSA_NUMERIC_FIXES_2026-07-19.md` 참조.

## 9. 고 SZA 층수 경고

```text
75-80 deg: n_layers < 100
80-84 deg: n_layers < 200
>=84 deg:  n_layers < 400
```

경고만 출력하고 자동으로 layer 수를 변경하지 않는다. SZA 85° 부근에서 총 TOA 수렴은 200층 이상, correction magnitude 자체의 sub-percent 수렴은 400-800층을 권장한다.

## 10. 검증 상태

완료:

- 하향 shell geometry 독립 US62 적분
- Eq. (8), surface, attenuation machine-precision invariant
- arbitrary-direction basis 대 table builder raw-bit identity
- Rayleigh/aerosol/PP/gas-only differential matrices
- black/Cox-Munk unaffected control byte identity
- sanitizer 및 기존 해수 회귀

미완료/범위 밖:

- 수중 PSSA
- 논문 CDISORT/AccuRT 수치 field 직접 재현
- SZA>85° 검증

상세: `PSSA_NUMERIC_FIXES_2026-07-19.md`.
