# OCRT 출력 변수 정의서 — TOA, 0+, 0-

## 좌표면

- `TOA`: 대기 상단, 센서 방향 상향 Stokes 및 반사도.
- `0+`: 해수면 바로 위 공기측.
- `0-`: 해수면 바로 아래 수중측.

## 정규화

- TOA `rho_*`: 무차원 방향 반사도(BRF 계열), OCRT Convention B에서 `rho = I_TOA/mu_s`.
- `Rrs0plus_* = {Lu,Qu,Uu}(0+)/Ed(0+)`, 단위 `sr^-1`.
- `rrs0minus_* = {Lu,Qu,Uu}(0-)/Ed(0-)`, 단위 `sr^-1`.
- `BRF0plus_* = pi*Rrs0plus_*`, `BRF0minus_* = pi*rrs0minus_*`, 무차원.
- 엄밀한 BRDF는 입사 방향별 복사휘도 경계조건을 기준으로 정의해야 한다. OCRT의 Rrs/rrs는 총 하향조도 정규화 방향반사량이며, BRDF 계산에 필요한 기본 변수 `L`과 `Ed`를 함께 출력한다.

## TOA 출력

| 변수 | 의미 |
|---|---|
| `TOA_rho_I/Q/U` | 직달 선글린트를 분리한 최종 TOA Stokes 반사도 |
| `TOA_rho_atm_I/Q/U` | 대기 경로복사 성분 |
| `TOA_rho_water_direct_I/Q/U` | 직달 태양광에 의해 생성된 수광신호의 TOA 기여 |
| `TOA_rho_water_sky_I/Q/U` | 대기 확산광에 의해 생성된 수광신호의 TOA 기여 |
| `TOA_rho_water_total_I/Q/U` | direct + sky 수광신호 |
| `TOA_rho_glint_direct_I/Q/U` | 직달 태양 Cox–Munk 선글린트 추가분 |
| `TOA_rho_glint_on_I/Q/U` | 선글린트 포함 최종 TOA 반사도 |
| `I_TOA/Q_TOA/U_TOA` | 내부 radiance 상태. 공개 stdout에서는 rho를 기본 제공 |

해양 결합 실행에서는 다음 합이 성립한다.

`TOA_rho = TOA_rho_atm + TOA_rho_water_total`

`TOA_rho_glint_on = TOA_rho + TOA_rho_glint_direct`

비해양 surface 실행은 별도 black-surface companion solve를 수행하지 않으므로 `TOA_rho_atm`은 현재 solved atmosphere+surface path 전체를 뜻한다. 표면·대기 성분을 완전히 분리하려면 추가 companion solve 옵션이 후속으로 필요하다.

## 대기투과율

| 변수 | 의미 |
|---|---|
| `T_dir_dn` | 태양 직달 하향 투과율. 비산란 직달광의 정확한 Beer–Lambert/PSSA 경로 |
| `T_diff_dn_hemi` | 대기 SOS가 계산한 BOA 하향 확산 radiance의 반구 적분값 |
| `T_total_dn_hemi` | `T_dir_dn + T_diff_dn_hemi`; 총 하향조도 투과율 |
| `T_diff_dn_dir` | 대기 SOS 하향장으로부터 재구성한 지정 방향 확산 성분 |
| `T_dir_up_view` | 지정 관측방향의 비산란 직접 상향 Beer–Lambert 성분 |
| `TOA_water_signal_I` | `I_TOA_total - I_atm_path`; 대기 SOS bottom-source pass가 전달한 수광신호의 TOA radiance |
| `T_total_up_view` | `TOA_water_signal_I / Lu0plus`; 실제 해수 BRDF를 보존한 RT 유효 상향투과율 |
| `T_diff_up_view` | `T_total_up_view - T_dir_up_view`; RT 유효 확산·각도재분배 성분 |
| `T_up_rt_valid` | 위 상향값이 엄밀 bottom-source SOS로 계산되었으면 1. 정의 불가 또는 fallback이면 0 |
| `diag_sunglint_up_flux_ratio` | `rho_glint*mu_sun`; 투과율이 아닌 legacy 진단비 |
| `diag_TOA_up_flux_ratio` | `rho_TOA*mu_sun`; 투과율이 아닌 legacy 진단비 |

상향 대기투과율은 대기의 고유 스칼라 상수가 아니다. 실제 `Lu(0+)`의 각도분포와
해수 BRDF에 의존하는 source-dependent 유효값이다. 따라서 하향 반구투과율을
상반성으로 복사하거나 `exp(-tau/mu_view)`에 더하지 않는다. 엄밀한 bottom-source
RT가 완료되지 않은 실행은 `T_diff_up_view`와 `T_total_up_view`를 `NaN`,
`T_up_rt_valid=0`으로 출력한다.

`T_diff_up_view`는 방향별 에너지의 각도 재분배를 포함하므로 음수가 될 수도 있다.
이는 총 투과율이 비물리적이라는 뜻이 아니라, 동일 view의 비산란 직접 성분에
대한 RT 산란 보정이 음수임을 뜻한다.

## 0+ 공기측

| 변수 | 의미 |
|---|---|
| `Lu0plus/Qu0plus/Uu0plus` | 지정 관측방향 상향 Stokes radiance |
| `Ed0plus` | 해수면 직상 총 하향조도 |
| `Ed0plus_direct` | 0+에 도달한 비산란 collimated solar-beam 하향조도 |
| `Ed0plus_diffuse` | 대기 확산 하향조도. `Ed0plus-Ed0plus_direct` |
| `Rrs0plus_I/Q/U` | 0+ Stokes 원격반사도 |
| `BRF0plus_I/Q/U` | `pi*Rrs0plus`, 무차원 방향반사계수 |

현재 `Eu0plus`는 출력하지 않는다. 이유는 air-side 전체 상향 반구를 별도로 재구성·적분하는 production 경로가 아직 없기 때문이다. 단일 view `Lu0plus`를 이용해 임의로 `Eu0plus`를 추정하지 않는다.

## 0- 수중측

| 변수 | 의미 |
|---|---|
| `Lu0minus/Qu0minus/Uu0minus` | 지정 관측방향 수중 상향 Stokes radiance |
| `Ed0minus` | 수면 직하 총 하향조도 |
| `Ed0minus_direct` | Fresnel/Snell 전달 후 비산란 collimated solar-beam 하향조도 |
| `Ed0minus_diffuse` | 수중 하향 확산조도. `Ed0minus-Ed0minus_direct` |
| `Eu0minus` | 수면 직하 총 상향조도 |
| `Eu0minus_direct` | 현재 deep-water/black-bottom 모델에서는 0 |
| `Eu0minus_diffuse` | 수면 직하 상향 확산조도. 현재 `Eu0minus`와 동일 |
| `rrs0minus_I/Q/U` | 수중 원격반사도 |
| `BRF0minus_I/Q/U` | `pi*rrs0minus` |
| `Kd0minus`, `Ku0minus` | 수면 직하 하향/상향 diffuse attenuation 진단값 |

### direct/diffuse의 의미

이 여섯 변수는 각도공간의 표준 radiometric 분해다. `direct`는 해당 경계면의 비산란 collimated solar beam이고, `diffuse`는 나머지 방향분포의 반구 적분이다. 태양 기원과 skylight 기원을 모든 다중산란 차수까지 분리하는 source-provenance 출력은 아니다.

다음 합이 내부 double 정밀도에서 성립한다.

```text
Ed0plus  = Ed0plus_direct  + Ed0plus_diffuse
Ed0minus = Ed0minus_direct + Ed0minus_diffuse
Eu0minus = Eu0minus_direct + Eu0minus_diffuse
```

상세 정의와 성능 정책은 `PSSA_SCOPE_AND_RADIOMETRY_SPLIT_2026-07-18.md`를 참조한다.

## IOP 출력

`a_w,b_w,bb_w,a_dom,a_pig,b_pig,bb_pig,a_min,b_min,bb_min,a_total,b_total,bb_total,omega_total`을 출력한다. 단위는 계수는 `m^-1`, omega는 무차원이다.

## 후속 추가가 필요한 변수

1. air-side 전체 상향 반구 적분 `Eu0plus`.
2. non-ocean surface 실행의 `rho_atm_path`, `rho_surface_diffuse`, `rho_surface_direct` 완전 분해.
3. 전체 방향·Stokes를 연결하는 상향 Mueller transmission operator의 선택적 공개 출력.
4. LUT/배치 CSV에 단일실행과 동일한 전체 변수 집합을 선택적으로 출력하는 `--output-radiometry` 모드.
