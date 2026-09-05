# FIX-SKY-EDLU: rrs(0−) skylight Ed/Lu sub-cone leak 일관성 수정 (2026-06-28)

## 1. 버그 (확정)

rrs(0−)가 atmosphere 추가 시 spurious하게 하락하며, **하락폭이 wind(해수면 거칠기)에 비례**한다.
물리적으로 rrs(0−)는 atm/wind에 거의 무관해야 한다(OSOAA flat; in-water IOP만의 함수에 가까움).

측정 (555 nm, sza=40, nadir, no-atm vs Rayleigh):
- wind 0: −0.23%
- wind 3: **−0.93%**
- wind 12: **−1.52%**
sza=80(grazing)에서는 더 컸다.

## 2. 원인 (코드 + 덤프로 규명)

rrs(0−) = Lu_diff(0−) / Ed_diff(0−). 둘이 **서로 다른 투과 skylight 총량**을 썼다:

- **Ed_diff(분모, rt_solver.c)**: air-side `surface_T_aw_coxmunk_direct` 적분.
  거친 해수면이 flat 임계각 너머(μ_w < μ_crit)로 투과시키는 **sub-cone leak을 포함**한다(정확한 투과 irradiance).
- **Lu_diff(분자)**: equivalent-beam superposition이 굴절 cone(μ_w > μ_crit)만 돌았다.
  air-refraction 등가빔은 sub-cone 방향에 대응하는 **실제 air 입사각이 없어**(TIR) `if (mu_w <= mu_crit) continue`로 건너뛰었다.

→ 분모는 sub-cone leak 포함, 분자는 누락 → rrs 하락. 거칠수록(wind↑) sub-cone leak↑ → 하락폭↑.
거친 해수면 결합(`rt_air_water_couple_atm_to_water`)은 이미 Cox-Munk BTDF를 쓰지만(line 176-208),
Lu 등가빔이 sub-cone를 **표현할 수 없다는 한계는 line 174 주석에 명시**돼 있었다.

정량 확인 (555, sza40, wind3): air-side Ed_diff = 4.16e-2 vs cone-only field Ed = 3.58e-2 (16% 차 = sub-cone leak).
sza80에서는 **sub-cone가 diffuse field의 25%**를 차지(grazing → near-horizon skylight 투과 큼) → 버그가 컸던 이유.

## 3. 수정 (direction B — 엄밀)

분자 Lu도 sub-cone를 포함하게 하여 분자/분모가 **동일한 전체 투과 field**의 물 응답이 되도록 함.
임의 튜닝/rescale 아님(Rule 3 준수). air-side Ed는 정확하므로 유지하고 Lu를 보완.

수정 위치:
1. `rt_water_rt.h`: `rt_water_rt_options_t`에 `double mu_sun_water_override` 추가(default 0).
2. `rt_water_rt.c` (wrapper `rt_water_rt_sos_pure`):
   - `mu_sun_water = (override>0) ? override : snell_down(...)` — override 시 in-water 빔 cosine 직접 사용(굴절 우회).
   - override 시 `beam_q = 0`(sub-cone 빔은 대응 air 입사각 없음 → unpolarized 처리; 아래 한계 참조).
3. `rt_solver.c`:
   - wind>0 시 in-water diffuse field grid를 **전반구 [0,1]**로(BTDF 결합이 sub-cone field도 채움).
     wind=0(flat)은 cone-map [μ_crit,1] 유지(TIR 아래 0, kink 회피 — BUG-WRT-003).
   - Lu 등가빔 루프에 **sub-cone 분기**: μ_w ≤ μ_crit이고 wind>0이면 건너뛰지 않고
     `mu_sun_water_override = μ_w`로 in-water 빔을 직접 구동.

## 4. 검증

**정규화 일관성**: full-grid field 적분 / air-side Ed_diff = **1.0024** (≈1.0) — 분자/분모가 동일 field에서 나옴(일관).

**wind 무관성 복원** (555, no-atm vs Rayleigh):
| wind | sza40 전 | sza40 후 | sza80 후 |
|---|---|---|---|
| 0 | −0.23% | −0.23% | +3.71% |
| 3 | −0.93% | **−0.14%** | +3.59% |
| 12 | −1.52% | **−0.05%** | +3.71% |
→ wind 의존성 제거됨(양쪽 sza). sza80의 +3.7%는 wind-무관 = 버그 아님(grazing서 직달빔 거의 없어 skylight가 in-water 장 지배 → no-atm vs with-atm rrs가 물리적으로 다름).

**OCRT-B vs OSOAA-with-atm rrs(0−)** (wind3, nadir):
| sza | 490 | 555 | 660 |
|---|---|---|---|
| 40 | **+0.05%** | **−0.33%** | −1.57% |
| 80 | +4.47% | +2.32% | −1.03% |
→ sza40 완벽 일치(490/555 ~0.3%). sza80 +2~4% 초과.

**TOA**: sza40 TOA_I MAPE 0.84% → **0.63%** (개선; 더 정확한 water-leaving Lu가 TOA에 기여).

## 5. 남은 한계 (별개 이슈, 본 수정과 무관)

- **sza80 +2% 초과**: 정규화는 일관(field/air=1.0024)하므로 정규화 오류 아님. 원인 = (1) OSOAA 자체 P7 grazing 부정확(TOA도 sza80서 3.10%), (2) sub-cone 빔의 unpolarized(beam_q=0) 근사 — Rayleigh skylight 강편광인 grazing에서 1차 효과가 됨. 완전 해결은 sub-cone 빔에 coupled.Q/U_inwater 전달 필요(향후).
- **412/443 +4.27%/+1.78%**: 본 수정이 skylight 버그(음의 기여)를 제거하니 **가려져 있던 high-omega 분자 위상/depolarization bias가 노출**됨(이전엔 +3.01%/−0.02%로 부분상쇄). coupled Rayleigh +2% mid-VZA@412와 동일 근원 추정. 별도 처리 대상.
- **660/865 −1.5%**: 고흡수·미소-rrs 밴드의 별개 잔차.

## 6. 결론

skylight Ed/Lu sub-cone leak 버그는 **확정·수정**됨(direction B). wind 무관성 복원, sza40 OSOAA 일치, TOA 개선,
정규화 일관. sza80 grazing 잔차와 노출된 분자위상 bias는 별개 후속 이슈.
