# ⛔ OCRT↔OSOAA 비교 실수 registry + 필수 절차 (세션 시작 시 **맨 먼저** 읽을 것)
v1.08 · 2026-06-28 · 이 문서를 안 읽고 비교하면 몇달째 같은 실수를 반복한다. 비교 전 §2 체크리스트 전부 통과 필수.

> **왜 이 문서가 존재하나:** OCRT↔OSOAA 비교가 몇달간 같은 설정/규약 실수로 가짜 물리버그(−20% 등)를 양산하며 진도를 못 냈다. 근본원인 2개와 그 해결을 아래에 박는다. **실수할 때마다 §3 registry에 append.**

---

## 0. 근본원인 2개 (반드시 내면화)
1. **기억/prose 참조값과 비교 = 금지.** 과거 세션이 적어둔 OSOAA 숫자(예: 구 nadir 2.5061e-2)는 그 자체가 **틀린 설정**으로 나온 것. 그걸 기준으로 OCRT를 맞추려다 −20% 같은 가짜 gap을 surface/MS로 오진단하고 몇 세션 낭비. **해결: OSOAA를 매 세션 sandbox서 fresh 실행(§4 recipe). 기억/prose/과거 CSV 숫자와 절대 비교 금지.**
2. **규약을 코드에 안 박음.** 매 세션 REFL/glint/U부호/Snell을 재유도하다 또 틀림. **해결: 규약을 canonical 스크립트(§5)에 박고, 새 비교는 그걸 호출.**

## 1. 절대 규칙 (위반 = 비교 무효, 결과 폐기)
- **R1.** OSOAA는 sandbox서 직접 실행. 기억/prose/과거 참조값 비교 금지.
- **R2.** **큰 gap(>~5%)은 거의 항상 설정/값 실수**다. 물리버그로 결론짓기 전 §2 체크리스트 전부 통과 + §3 registry 대조. 수렴된 벡터 RT 둘은 보통 ~1% 일치한다.
- **R3.** 임의 튜닝 금지. gap은 **규약/설정 교정으로만** 줄인다(Rule 3).
- **R4.** 큰 gap을 풀 땐 **분해/이론 전에** 먼저 OSOAA 설정 변수(아래 registry)부터 1개씩 토글해 원인 격리.

## 2. PRE-FLIGHT 체크리스트 (결과 신뢰 전 전부 ✓)
```
[ ] OSOAA fresh run (이번 세션 sandbox 빌드), NOT 기억/prose 숫자
[ ] 비교 물리량: OCRT rho_I ↔ OSOAA REFL 컬럼
        REFL = π·L/Ed (표준 반사도).  I 컬럼 = π·L/Esun = REFL×μ_sun (다른 값! 혼동시 ~13%/cos오차)
[ ] OSOAA direct sunglint 억제: env OSOAA_NO_DIRECT_GLINT=1  (OCRT는 glint-free 출력)
        안 하면 specular 근처(특히 sza0/vza0)서 −50~−98% 가짜 gap
[ ] matched 물리입력 전부:
        - Rayleigh: OSOAA -AP.MOT = OCRT bodhaine tau_R (--dump-atm으로 추출)
        - aerosol: -AER.AOTref = OCRT --aod(밴드값), -AP.HA = OCRT --aer-h-km, ExtData ssa/phase 일치
        - truncation: OCRT --trunc-aer-loglin ↔ OSOAA -AER.Tronca 1 (방법差 작음, gap시 의심)
[ ] VZA frame: TOA(Level1)=air-side, Snell 불필요. 0−(Level4)=water-side, air=arcsin(1.34·sin(water))
[ ] RAA 매핑: OCRT raa=90 = OSOAA Phi=90 (convention-invariant, 여기부터). 그 외 raa는 SCA(Θ) 일치로 검증
[ ] U 부호: OCRT_U = −OSOAA_U (handedness 규약). 적용 후 비교
[ ] Q at nadir(vza=0): meridian 평면 특이점 → 부호 모호. |Q| 비교, 부호비교는 vza>0만
[ ] 보간: OSOAA Gauss격자 → target VZA는 PCHIP (nearest 절대 금지)
[ ] 기하 축퇴: sza=0 또는 vza=0이면 raa 축퇴 — 대표 1개만, lookup시 fallback
```

## 3. 실수 registry (증상→원인→교정) — **실수할 때마다 여기 append**
| # | 증상 | 진짜 원인 | 교정 | 발견 |
|---|---|---|---|---|
| P1 | aerosol nadir −20% (OCRT 2.003e-2 vs "2.506e-2") | 구 참조 2.506e-2가 **틀린 값**(잘못된 설정의 산물) | OSOAA fresh REFL+glint-free → MAPE 0.56% | 2026-06-28 |
| P2 | sza0/vza0서 −54~−98% | OSOAA **direct sunglint** 포함(OCRT는 glint-free) | OSOAA_NO_DIRECT_GLINT=1 → 0.4% | 2026-06-28 |
| P3 | U 부호 전부 반대(24/24) | **handedness 부호 규약**差 | OCRT_U = −OSOAA_U | 2026-06-28 |
| P4 | Q at vza=0 d=−198% | nadir **meridian 특이점**(Q 부호 모호) | |Q| 비교, vza>0만 부호 | 2026-06-28 |
| P5 | (잠재) ~13% 일정 offset | OSOAA **I 컬럼**(=REFL×μ_sun) 사용 | REFL 컬럼 사용 | 규약 |
| P6 | "OSOAA sandbox서 못 돌림/파서 못 만듦"으로 비교 회피 | full source는 **GitHub**서 받아 patch+빌드 가능; 파서는 패키지 parse_osoaa.py에 있었음 | §4 recipe로 빌드, §6 포맷 | 2026-06-28 |
| P7 | sza=80서 gap(Rayleigh −2~−7%@λ↑; aerosol thin-aot 6%→thick-aot 0.5%) | **미해결, surface 확정** — Rayleigh·aerosol 공통, **aot↑에 감소** = surface BRDF 항(atmosphere RT 아님). grazing Cox-Munk 확산반사/skylight 결합 의심 | 조사 대기 (sza≤40은 영향 없음) | 2026-06-28 |

| P8 | ray_aer서 −85% (OCRT≪OSOAA) | 비교 스크립트가 **OCRT 잘못된 atm subset 로드**(aerosol-only를 ray_aer와 비교; sed가 `r['atm']==` 안 바꿈) | 로드 atm 필터 확인 → ray_aer | 2026-06-28 |
| P9 | (과거) "IOP parity b_w 1.67e-3" | **실제 OSOAA b_w(555)=1.859e-3 = Morel 계열**(n=4.29), memory의 z09 1.67e-3가 아님 → 과거 OCRT(1.67e-3)와 OSOAA(1.859e-3)는 **~11% b_w 불일치**였을 가능성. a_w는 Pope&Fry1997로 정확 일치 | water IOP LUT를 **OSOAA 실측 6밴드값으로 재구성**(2026-06-28), OCRT 출력 a_w/b_w/ω가 OSOAA와 exact match 확인 | 2026-06-28 |
(| P10 | Rrs(0−) 0+ 비교가 ~4배 틀림 | OSOAA Adv level 26(0+)는 **surface 반사 skylight 포함**, OCRT Rrs0plus는 water-leaving만 | **0−(level 27)에서 비교**(surface 반사 없음): OCRT rrs0minus ↔ OSOAA I(27)/Ed(27) | 2026-06-28 |
| P11 | rrs(0−)가 atm 추가 시 **wind 비례 하락**(wind3 −0.93%, wind12 −1.52%@555/sza40) | **FIX-SKY-EDLU 해결됨**. skylight Ed/Lu sub-cone leak 불일치: Ed_diff(분모)는 air-side 적분으로 거친면이 μ_w<μ_crit로 투과시키는 sub-cone leak 포함, Lu(분자) equivalent-beam은 cone(μ_w>μ_crit)만 — air-refraction 등가빔이 sub-cone 표현 불가 | **direction B**: wind>0 시 in-water field grid 전반구[0,1] + Lu 루프 sub-cone 분기(`mu_sun_water_override`로 in-water 빔 직접 구동). 검증: wind 무관성 복원, field/air-side Ed=1.0024 일관, sza40 OSOAA 일치(490/555 ~0.3%), TOA 0.84→0.63%. **sza80 +2% 잔차는 별개**(OSOAA P7 + sub-cone unpolarized 근사). CHANGELOG_FIX-SKY-EDLU_2026-06-28.md | 2026-06-28 |

P7은 아직 open. P8 = 단일물량의 가장 흔한 실수: 큰 gap이면 **먼저 어느 subset/값을 로드했는지** 확인.

## 4. OSOAA sandbox 빌드 recipe (매 세션 — 파일시스템 초기화됨)
```bash
cd /home/user && git clone --depth 1 https://github.com/CNES/RadiativeTransferCode-OSOAA.git osoaa_build
# 수정 4파일 덮어쓰기 (Jae 패키지 osoaa_pkg/ 또는 uploads zip)
cp osoaa_pkg/OSOAA.h osoaa_build/inc/ ; cp osoaa_pkg/OSOAA_{MAIN,PROFILE,TRPHI}.F osoaa_build/src/
apt-get install -y gfortran
export OSOAA_ROOT=/home/user/osoaa_build && cd $OSOAA_ROOT/gen && make -f Makefile_OSOAA.gfortran
mkdir -p $OSOAA_ROOT/DATABASE/SURF_MATR   # surface matrix 첫 run 자동 캐시
```
- 버전: master = OSOAA V2.0(2025-01-30) = patch base 일치. patch는 full-file 교체라 클린 적용.
- gfortran 13.x. 빌드 ~1분. 단일 OSOAA run ~1s(matrix 캐시 후), 첫 sza/wind은 matrix 생성으로 느림.

## 5. canonical 비교 절차 (executable — 규약 내장, 새 비교는 이걸 호출)
1. **ExtData 생성**(aerosol): `mie_to_osoaa_extdata.py <mie> <band_nm> <out>` — EXTINCTION/SCATTERING_COEF, NB_LINES, `ANGLE F11 -F12/F11 F22/F11 F33/F11`(오름차순). 구면: F22/F11=1, F33/F11=P33/P11. 밴드 비정확매치는 wl-PCHIP 보간(검증 요).
2. **OSOAA 실행 필수 플래그**(빠지면 에러): `-AER.DirMie <dir>`(필수!), `-AER.Model 4`(external ExtData), `-AER.ExtData`, `-AER.Tronca 1`, `-AER.Waref <µm>`, `-AER.AOTref <aod>`; Rayleigh는 `-AP.MOT <tauR>`; `-OSOAA.View.Level 1`(TOA)/`4`(0−); `-YS.Abs440 1000`(black water); `-SEA.Wind 3 -SEA.Ind 1.34 -SEA.Dir $SURF`; env `OSOAA_NO_DIRECT_GLINT=1`.
3. **추출**: ρ_I = vsVZA **REFL**; ρ_Q = Adv_UP(level0=TOA) Q/μsun; ρ_U = Adv_UP U/μsun; PCHIP로 target VZA 보간.
4. **비교 스크립트**: `osoaa_ray_compare.py`(Rayleigh, 규약 내장), aerosol은 동형 확장. U 부호반전·축퇴 fallback 포함.

## 6. OSOAA 출력 포맷 (파싱 — parse_osoaa.py 기반)
- **Standard/RESLUM_vsVZA.txt**: `VZA SCA_ANG I REFL POL_RATE LPOL REFL_POL` (7컬럼). **REFL=π·L/Ed ← 이걸 OCRT rho_I와 비교.** VZA<0=azimuth(Phi+180), VZA>0=azimuth Phi.
- **Advanced/RESLUM_Adv_UP.txt**: `LEVEL Z VZA SCA I Q U ...`. **level 0=TOA, 26=0+, 27=0−.** I/Q/U = π·L/Esun (ρ는 /μsun). 'D'→'E'.
- **Advanced/Flux.txt**: `Level Z Direct_Down Diffuse_Down Total_Down(=Ed) ...`. Esun_TOA=π.
- **RAA 규약(출력에 명시)**: Phi=180 ↔ 위성·태양 same half-plane, Phi=0 ↔ opposite. nadir SCA≈158°(0−)/150°(TOA)로 water/air-side 확인.

## 7. 검증된 사실 (재유도 금지 — frozen)
- OCRT bodhaine tau_R(P=1013.25): 412=0.318540221, 443=0.236054530, 490=0.155974381, 555=0.093751620, 660=0.046362496, 865=0.015540855. (`--dump-atm --wavelength W --pressure 1013.25 --sza S --vza 0 --raa 90`)
- M80C 밴드별 AOD ratio(aot865→밴드): 412=1.1309,443=1.1152,490=1.0937,555=1.0691,660=1.0409,865=1.0000.
- M80C@412 ExtData: ext=5.0667e-2, sca=5.0480e-2, ssa=0.99631, 361각.
- 결과: Rayleigh sza≤40 MAPE ~0.4%; aerosol M80C@412 sza30 MAPE 0.56%. **−20%는 폐기된 아티팩트.**

## 8. water IOP LUT (재구성, inputs/water_iop/ — #2~#6 in-water RT 필수)
- **OCRT가 요구**: `water_coef_z09_1nm.txt`(3컬럼 λ_nm a_w b_w), `psi_T_rottgers2014_OCRT.txt`(2컬럼 λ_nm psi_T). 없으면 `--surface ocean` 실행 불가.
- **로더 포맷 함정**: 헤더는 `#` 주석이 아니라 **`/end_header` 줄로 종료** 표시 필수. 데이터는 그 뒤. (이거 빠지면 "too few data rows (0)".)
- **재구성 방법(2026-06-28)**: OSOAA pure-water PROFILE_SEA.txt에서 6밴드 a_w/b_w 추출(c=dTAU/d깊이, b_w=MOL_PC×c, a_w=c−b_w) → 1nm log-PCHIP. **6밴드서 OSOAA와 exact match**(검증됨). a_w=Pope&Fry1997, b_w=Morel-like(n=4.29).
- **검증**: `--surface ocean --simple-chl 0 --pressure 0`서 출력 a_w/b_w/omega_total이 OSOAA와 일치. 555: a=5.960e-2 b=1.859e-3 ω=0.03025.
- OCRT 출력은 named-field(comma 아님): rho_I/Q/U(첫 3토큰=TOA), rho_I_glint, Lu0±, Ed0±, Rrs0plus/rrs0minus, Kd/Ku, a_w/b_w/bb_w, a_total/b_total/omega_total.
