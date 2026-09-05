# item#4 aDOM(CDOM/yellow substance) 검증 (2026-06-29, sza40)

## 약어
aDOM/CDOM = yellow substance(용존유기물) 흡수, a_CDOM(440) = 440nm 흡수계수[m⁻¹], slope = 지수 spectral slope[nm⁻¹], a_total = 총흡수(a_w + a_CDOM), rrs(0−) = 수면 바로 아래 remote-sensing reflectance, TOA ρ_I = 대기상한 반사도 I, MAPE = mean absolute percentage error, YS = OSOAA yellow substance.

## 셋업 & parity
- CDOM = **흡수만**(b=bb=0). 산란·위상·delta-M 절단 없음 → in-water 흡수 응답만 격리 검증.
- 모델 동일: OCRT `a_CDOM(λ)=a_CDOM(440)·exp(−S·(λ−440))` == OSOAA `YS_A440·exp(−YS_SWA·(λ−440))` (OSOAA_PROFILE.F line 864).
- parity 인자: OCRT `--cdom-a440 X --cdom-slope 0.014` == OSOAA `-YS.Abs440 X -YS.Swa 0.014`. DET=SED=PHYTO.Chl=0.
- IOP 검증: OCRT a_total − a_w = a_CDOM이 formula와 정확 일치(412/adom0.1: 0.14799 vs 0.147994). 물 IOP는 #2와 동일(OCRT 기본 = OSOAA).
- aDOM(440) = {0.01, 0.1, 1.0} m⁻¹, slope 0.014, sza40, vza{0,30,60}, raa90, 6밴드, wind3.

## 결과 (sza40)

**TOA ρ_I**: 전 aDOM MAPE **0.49~0.57%** (일관 양호; 1:1선 일치).

**rrs(0−) nadir bias** (흡수 클수록 음의 bias 증가):
| band | adom=0.01 | adom=0.1 | adom=1.0 |
|---|---|---|---|
| 412 | +0.24% | −1.74% | −1.83% |
| 443 | +0.13% | −0.97% | −1.78% |
| 490 | +0.03% | −0.42% | −1.71% |
| 555 | −0.34% | −0.52% | −1.54% |
| 660 | −1.65% | −1.73% | −1.91% |
| 865 | −1.52% | −1.53% | −1.55% |

평균 rrs bias: adom0.01 −0.52%, adom0.1 −1.15%, adom1.0 −1.72%.

## 핵심 발견: 잔차는 흡수 구동(constituent 무관)

rrs bias vs 총흡수 a_total을 그리면(그림 Panel C) 모든 aDOM/band 점이 **한 추세로 모인다**: 저흡수(a_total~0.01) ~0% → 고흡수(a_total>0.3) ~−1.8%.
- 412가 pure water서 +4.27%(분자위상 bias, Rayleigh)였으나 CDOM 강흡수 추가 후 −1.83%로 전환 = **고흡수 regime 진입**.
- #2의 660/865 −1.5% 잔차와 **동일 정체**: in-water RT **고흡수 regime 잔차**. band·constituent가 아니라 **a_total로 결정**.
- 즉 skylight 버그(해결)·Rayleigh 분자위상(OSOAA 의심, 보류)과 **별개인 제3의 잔차**.

## 결론
- aDOM 검증: TOA 양호(0.52%), Rrs는 저흡수서 sub-1%, **고흡수서 ~−1.5~−1.9%**(흡수 구동 잔차).
- 이 고흡수 잔차는 chl(red band)·복합에도 나타날 cross-cutting 이슈 → 별도 규명 필요(in-water RT 고흡수 처리: depth/layering/single-scatter 균형 의심).
- 남은 것: sza0/sza80 geometry(이번은 sza40만), 고흡수 잔차 규명.
