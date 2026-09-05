# OCRT ↔ OSOAA 비교 검증 가이드라인 (통합본)
버전 v1.08 · 2026-06-27 · 출처: README(harness), handoff rev3, OSOAA MODS MANIFEST, userMemories, 2026-06-27 세션 findings

> 이 문서는 흩어져 있던 비교 규약·matched config·검증 방법론을 **하나로 통합**한 것이다. 매 비교 세션 시작 시 먼저 읽는다. 비교 파라미터는 여기 frozen 값에서만 쓰고, 벗어나면 명시·정당화한다(묵시적 변경 금지).

---

## 0. 두 비교 케이스 — 규약이 다르다 (최대 함정)

| 구분 | Case A: aerosol/atmosphere TOA | Case B: in-water rrs(0−) |
|---|---|---|
| 비교량 | TOA reflectance rho_I | rrs(0−) |
| OSOAA View.Level | **1 (TOA)** | **4 (0− just below surface)** |
| VZA 면 | **air-side** | **water-side** |
| **Snell 변환** | **불필요** (양쪽 다 air-side) | **필수** (water→air) |
| 대기 | aerosol-only(`--pressure 0`) 또는 coupled | pure-water(무대기 또는 Rayleigh) |
| 현 용도 | aerosol phase −20% 디버깅 | frozen baseline (regression guard) |

**Snell 함정**: Case B에서 OCRT `--vza`(air-side)와 OSOAA Level4 VZA 컬럼(water-side)을 직접 비교 금지. 변환: `air_VZA = arcsin(n_water · sin(water_VZA))`, n_water=1.34. nadir(0°)만 우연히 일치 → nadir 외에서 어긋나면 VZA 해석오류 신호. Case A(TOA, Level1)는 둘 다 air-side라 변환 불필요.

---

## 1. 공통 규약 체크리스트 (모든 비교에 적용)

- **(A) VZA 면**: Case A air-side 직접 / Case B water→air Snell(n=1.34). ← §0 함정.
- **(B) glint**: wind>0이면 **sunglint-free RAA**에서만 비교, 기본 **RAA=90**(nadir vza=0에선 glint가 기하적으로 부재 — glint는 vza=sza). glint가 RAA=0인지 180인지는 **convention 의존**이고 코드(V3: specular at raa=180)와 memory 체크리스트(RAA=0 glint)가 **충돌**(§1.1 ⚠) → surface 비교 시 경험적 확인 필수. 오염 max ~6%.
- **(C) wind**: OSOAA는 **wind 3/5/10만 유효**(wind=0 미구현). wind=0 비교 금지.
- **(D) 산란차수 IGmax**: **무제한(full, 권장 IGmax=100)**. IGmax=1은 단일산란만이라 부적합(−22%).
- **(E) 휘도 규약**: **OCRT rho_I = OSOAA REFL 컬럼** (= π·L/Ed = 표준 reflectance ρ). OSOAA의 `I` 컬럼(= π·L/Esun = ρ·μ_sun)은 **아님**. 혼동 시 μ_sun=cos(sza) factor만큼 어긋남(sza30 → −13.4%). ※ 이 차이는 **RAA에 무관**(μ_sun factor는 전 RAA에 존재) — azimuth convention-invariance(§1.1)와는 별개 사안.
- **(F) 보간**: OSOAA Gauss VZA 격자 → OCRT vza는 **linear 보간**(nearest 금지).
- **(G) 참조도 틀릴 수 있음**: OSOAA를 무조건 ground truth로 두지 말 것. 불일치 시 양쪽 다 물리적으로 옳은지 확인. "OCRT가 틀렸다" 단정 금지.

### 1.1 RAA/azimuth 규약 — OCRT ↔ OSOAA 180° 차이 (★재발 함정 · 과거 거짓 "OCRT 버그" 판정 유발)

**확정 사실 (코드 검증):**
- OCRT 공개 CLI `--raa`는 V3 convention. **전 경로에 일관된 +π(180°) azimuth shift**: 물(rt_water_rt.c:139 `phi_arg=raa+π`, :3001 `phi_v=π−RAA`), 대기 TOA(rt_solver.c:863 `V3 raa=(raa_AF+180) mod 360`; :859 specular at raa=180). 6SV 검증된 **의도적 규약 — 코드 정상, 버그 아님**.
- 즉 물리 Δφ = OCRT_raa + 180° (V3↔AF). OSOAA `View.Phi`는 자체 relative azimuth.

**핵심 귀결 — RAA=90는 invariant:**
φ↔180−φ(또는 +180) 관계의 고정점이 90° → **RAA=90에서 OCRT·OSOAA azimuth 규약 차이가 소거된다.** 이게 RAA=90을 표준 비교 방위로 쓰는 진짜 이유(glint-free + **convention-invariant**). 현 frozen config(Case A·B 둘 다 RAA=90)에선 이 차이가 안 문다 — 단 **non-90 RAA 비교 시 반드시 매핑·검증**.

**non-90 RAA 매핑 (과거 기록 — 출력 레벨별로 다름, 경험적 확인 전제):**
- in-water rrs(0−) air-side: OCRT --raa = OSOAA Phi (flip 없음) [v1.06].
- above-water Rrs(0+): OCRT --raa = 180 − OSOAA Phi [v1.06].
- in-water Stokes Q: OCRT --raa = OSOAA Phi − 180 [v1.07-2].
- 대기 coupled: v1.07-5는 "View.Phi = RAA_std identity" 보고(empirical: View.Phi=45 → Θ=115.86° 일치). ← 위 in-water 매핑들과 표면상 충돌 = **레벨 의존 + 버전별 정제** 결과이니 단일 규칙으로 단정 금지.

**경험적 확인 방법 (non-90 비교 전 필수):** 테스트 geometry에서 **산란각 Θ를 양 코드로 계산해 일치** 확인(예 OSOAA View.Phi=45 → Θ=115.86°). glint는 specular 기하(vza=sza, 해당 azimuth)로 확인. 라벨(forward/back/glint at 0 vs 180) **암기 금지 — 매번 Θ로 검증**.

**⚠ 미해결 충돌:** 코드(V3: glint at raa=180)와 §1B 출처 memory("RAA=0 glint")가 충돌. 서로 다른 convention 참조로 보이나 **확정 필요**. glint azimuth는 경험적으로 잡을 것.

**역사적 교훈:** 이 규약을 놓쳐 과거 OSOAA가 "OCRT outlier(버그)"로 오판(v1.07-2). 실제 원인은 OCRT 측 azimuth 라벨링 오류 — OCRT 코드·OSOAA 판정 모두 정상이었음. → 모든 Stokes/non-90 비교에서 이 매핑 반드시 적용·Θ검증.

---

## 2. Matched configuration (frozen — 재유도 금지)

### 2.1 Case A: aerosol-only TOA (현 active 디버깅)

약어: AOD=aerosol optical depth, HA=aerosol scale height(km), HR=Rayleigh scale height(km), MOT=molecular optical thickness, Tronca=forward-peak 절단 플래그.

**OCRT 명령** (value/moment 각각):
```
./build/v2_solver_vk --surface coxmunk --wind-speed 3 --sza 30 --vza 0 --raa 90 \
  --wavelength 412 --pressure 0 --mie inputs/M80C.mie --aod 0.2 \
  --aer-l-max 80 --m-max 16 --aer-h-km 2.0 --n-layers 400 --trunc-aer-loglin \
  --aer-phase-kernel value --vector --lut --lut-vza-step 10 --lut-vza-max 80 \
  --lut-raa-step 30 --output-full-grid /tmp/vk_value.csv
```
moment 비교는 `--aer-phase-kernel moment`로 동일 실행.

**필수 3플래그 (coupled/aerosol 정합, 빠지면 틀린 LUT)**:
- `--aer-h-km 2.0` : OCRT 기본 수직프로파일은 6SV an23(V=23km)라 OSOAA `AP.HA 2.0`과 불일치 → 2km exp로 맞춤.
- `--n-layers 400` : 기본 40은 coupled under-converged.
- `--trunc-aer-loglin` : **핵심 함정** — loglin forward-peak 절단이 aerosol-only(Rayleigh OFF)에서만 자동 활성(main.c:1349 `cs.aerosol_on && !cs.rayleigh_on`). **coupled에선 자동 안 켜져** Tronca=0이 되어 OSOAA Tronca=1과 불일치 → coupled 비교 시 반드시 명시.

**OSOAA 명령**:
```
export OSOAA_ROOT=$(pwd)
python <harness>/mie_to_osoaa_extdata.py M80C.mie 412 /tmp/M80C_412.extdata
exe/OSOAA_MAIN.exe -OSOAA.Wa 0.412 -OSOAA.View.Level 1 -OSOAA.View.Phi 90 -SOS.IGmax 100 \
  -AER.Model 4 -AER.ExtData /tmp/M80C_412.extdata -AER.Tronca 1 -AER.Waref 0.412 -AER.AOTref 0.2 \
  -AP.HA 2.0 -AP.HR 8.0 -AP.Pressure 0.0 \
  -SEA.Dir $OSOAA_ROOT/DATABASE/SURF_MATR -SEA.Ind 1.34 -SEA.Wind 3 -SEA.Depth 100 \
  -YS.Abs440 1000 -PHYTO.Chl 0 -SED.Csed 0 -DET.Abs440 0 \
  -OSOAA.Log run.log -OSOAA.ResRoot /tmp/OSOAA_AERONLY_M80C_412
```
규약 매핑: `OCRT --aod` ↔ `AER.AOTref`(같은 파장), `OCRT --wavelength nm` ↔ `OSOAA.Wa`/`AER.Waref`(µm), 흑수면 = `YS.Abs440 1000` + Chl/SED/DET 0.

**coupled(Rayleigh ON)일 때만 변경**: OCRT `--pressure 1013.25`, OSOAA `-AP.Pressure 1013.25`. (aerosol-only가 기본 — coupled@412는 Rayleigh가 nadir 지배해 aerosol backscatter 결함을 가린다.)

### 2.2 Case B: pure-water 0− (frozen baseline, regression guard)

**고정 파라미터** (매 세션 재유도 금지):
- 파장 555 nm, wind 3, SZA 30°, RAA 90(glint-free), **depth 200**(5000 금지: layer 이산화 거칠어 rrs +0.8% 과대 / 1000도 baseline과 다름), **chl 0.0**(0.00001 금지: +0.13%).
- VZA: air-side Snell `air = arcsin(1.34·sin(water_VZA))`.
- **OSOAA rrs(0−) 추출 = Adv_UP level27 의 I / Flux level27 의 Ed** (출력 profile의 0− 레벨 행. `-OSOAA.View.Level 4` 인자와 별개의 출력 레벨 인덱스).
- IGmax 무제한.
- **IOP parity 확인**(필수): ω 0.02728, b_w 1.67e-3, a_w 5.96e-2.

baseline 산출물(FROZEN 2026-06-25): `OCRT_purewater_atm_BASELINE_golden_v1.08_2026-06-25.csv`, `ocrt_purewater_atm_REGRESSION_GUARD.py`, `OCRT_BASELINE_RULES_purewater_atm_2026-06-25.md`. Gate: 정확도가 baseline보다 나빠지면 FAIL, 좋아지면 golden 재동결.

---

## 3. OSOAA 측 셋업 (patches + 변환)

### 3.1 적용 patches (pristine OSOAA V2.0 대비, 재컴파일 필수)
| Patch | 파일 | 내용 | 거동 |
|---|---|---|---|
| P1 | inc/OSOAA.h | `CTE_MAXNB_ANG_EXT` 200→**400** | 배열 상한만. .mie 361각도 수용(<361이면 ERROR_950). 물리 불변 |
| P2 | src/OSOAA_TRPHI.F | env `OSOAA_NO_DIRECT_GLINT` → direct glint 제외 | **default-OFF**(안 주면 정상 glint). glint 분리 진단용. RAA=90선 무관 |
| P3a | src/OSOAA_PROFILE.F | Rayleigh bw를 **Z09 table**로(항상-on) | **물리 변경**(pristine과 Rayleigh 결과 다름). OCRT parity 의도 |
| P3b | src/OSOAA_PROFILE.F | user-depth J-1 레벨 선형 보간 | depth 격자 정합 |
| P4 | src/OSOAA_MAIN.F | seawater index OCRT parity 기본값(1.34) | `-SEA.Ind` 명시하면 무관 |

⚠ P3a 함의: Rayleigh bw는 이미 일치 → 잔여 Rayleigh @412 차이는 **bw 아니라 depolarization factor + Rayleigh phase matrix**에서 옴(어느 쪽이 옳은지 미해결, memory #19). aerosol-only(`AP.Pressure 0`)엔 무관.

### 3.2 `.mie` → OSOAA AER.ExtData 변환 (검증된 원본 `mie_to_osoaa_extdata.py`)
- **사용**: `python mie_to_osoaa_extdata.py <file.mie> <wavelength_nm> <out.extdata>`
- **컬럼**: `ANGLE  F11  −F12/F11  F22/F11  F33/F11` (구형 가정 F22=F11=1).
- **핵심 규약**: col3 = **−F12/F11** (OSOAA가 −F12/F11을 읽음; F12<0 for Rayleigh). 각도 **0→180 오름차순**(.mie는 180→0 내림차순이라 역순으로 씀).
- ExtData 컬럼 레이아웃의 ground truth — 추측 재작성 금지.

### 3.3 OSOAA 운영 주의
- `OSOAA_ROOT` env 필수. `-OSOAA.Log`는 **relative path만**(absolute → ERROR_900).
- 전체 배포본(OSOAA_SOS.F, OSOAA_SURFACE.F, Makefile, DATABASE/SURF_MATR 등)은 별도. 현 패키지엔 수정 4파일 + harness만 → **sandbox 빌드 불가**, Jae 로컬 실행.
- surface는 precomputed Cox-Munk BRDF matrix(`DATABASE/SURF_MATR`)에서 옴.

---

## 4. OCRT 측 셋업

- **빌드**: `gcc -std=c11 -O2 -fopenmp -Isrc $(find src -name '*.c') -o build/v2_solver_vk -lm` (함정: `src/*.c shared/*.c` 금지 — shared는 src/shared/에 있음).
- **대기 제어**(v1.08): 토글 없음, 물리량만 — Rayleigh는 `--pressure`(기본 1013.25; τ_R∝P/1013.25; `--pressure 0`이면 없음), aerosol은 `--aod`. **기본 run에 대기 ON** → in-water 통제엔 `--pressure 0` 명시 필수.
- **phase kernel**: `--aer-phase-kernel moment`(기본, Legendre) / `value`(angle-space raw LUT). value는 2026-06-27 세션서 2버그 수정 완료.
- **진단 env-gate**(전부 default-off): `OCRT_DEBUG=1` 필수 전제. `OCRT_DUMP_PFM=1`(phase_fourier 덤프), `OCRT_DUMP_SKY=1`(skylight). 신규 env var는 반드시 default-off·env-gate, production path 우회 금지.

---

## 5. 2-layer 검증 방법론 (재발 방지의 핵심)

과거 "풀 RT 복사휘도 비교"만으로는 **coupled@412서 Rayleigh가 aerosol backscatter 결함을 가려** 통과해버렸다. 두 층으로 막는다.

**층 A — phase-consistency (풀 RT 불필요, 빠름, 신규 aerosol 작업 전 필수 gate)**
- `mie_phase_consistency.py <file.mie> <wl> [L...]`: raw P11을 OCRT 규약 Σ(2l+1)βₗPₗ로 재구성 → backscatter(Θ≥120°)서 음수·진동·비수렴 감지. exit 0=PASS, 1=FAIL.
- moment-only면 **FAIL이 정상**(결함 노출). 예 M80C@412: L40 min=−0.375, L80 min=−0.0199(음수=Gibbs), GATE FAIL.
- ⚠ 한계: native μ 격자 sparse라 고차 βₗ 정량 부정확(정성 음수 판정은 견고). dense-θ PCHIP 내삽 후 βₗ 재계산이 TODO.

**층 B — 복사휘도 (풀 RT, aerosol-only)**
- `compare_ocrt_osoaa_aer.py <ocrt_value.csv> <ocrt_moment.csv> <osoaa_LUM_vsVZA.txt> [golden_d%]`: VZA별 REFL·moment/value·d%·MAPE·nadir gate·1:1 scatter(linear).
- **반드시 aerosol-only**(`--pressure 0` / `AP.Pressure 0`).
- ⚠ `VZA_COL=0, REFL_COL=3`은 가정 → 실제 `LUM_vsVZA.txt` 헤더로 확인.

**워크플로우**: ①신규 작업 전 층 A → ②수정 후 층 A2(effective phase, TODO) → ③층 B 전 VZA 수렴+scatter → ④gate 통과=체크포인트(golden 수치파일+스크립트가 유일한 비교, prose 수치는 baseline 아님).

---

## 6. 결과 테이블 (template + 현 known values)

> **단일 golden 결과테이블은 아직 없다.** harness가 run별로 생성한다. 아래는 현 시점 known reference + 이번 세션 측정. **OSOAA 참조값은 재확인 필요**(§7 참조 — config 정합성 미확정).

### 6.1 VZA sweep (Case A: M80C, 412nm, aerosol-only, sza30, raa90, aod0.2, n-layers400)
약어: VZA=view zenith angle(°), d%=100·(OCRT/OSOAA−1), moment=Legendre kernel, value=angle-space kernel(2026-06-27 수정).

| VZA | OSOAA REFL | OCRT moment | d% | OCRT value | d% |
|---:|---:|---:|---:|---:|---:|
| 0 | 2.5061e-2 | 2.0033e-2 | −20.1 | 2.0087e-2 | −19.8 |
| 10 | 2.2669e-2 | 1.9681e-2 | −13.2 | 1.9642e-2 | −13.4 |
| 20 | 1.8888e-2 | 1.8157e-2 | −3.9 | 1.8158e-2 | −3.9 |
| 40 | 1.3054e-2 | 1.2956e-2 | −0.8 | 1.2957e-2 | −0.7 |
| 60 | 1.8181e-2 | 1.8015e-2 | −0.9 | 1.8014e-2 | −0.9 |
| 80 | 5.1645e-2 | 5.0331e-2 | −2.5 | 5.0337e-2 | −2.5 |
| MAPE | | | 6.89% | | 6.88% |

### 6.2 nadir 분해 (Case A, AOD=0.2, 2026-06-27 세션) — OSOAA-독립 진단
약어: SS=single scattering(해석적, OCRT 일치 확정), atm_MS=atmosphere multiple scattering, surface=Cox-Munk Fresnel reflection + 결합 MS.

| 항 | rho_I | cox_full 대비 |
|---|---:|---:|
| SS (해석적, 정확) | 1.0789e-2 | 53.9% |
| atm_MS | 2.3908e-3 | 11.9% |
| surface | 6.8531e-3 | 34.2% |
| **OCRT cox_full** | **2.0033e-2** | 100% |
| OSOAA ref | 2.5061e-2 | — |
| gap | 5.03e-3 | OSOAA +25.1% |

### 6.3 결과테이블 권장 컬럼 (신규 비교 시)
`case | particle | wl_nm | sza | raa | wind | aod | vza | OSOAA_REFL | OCRT_moment | OCRT_value | d%_moment | d%_value` + 행별 d% + MAPE + nadir gate(golden 대비).

---

## 7. ⚠ 이번 세션(2026-06-27) 정정사항 — 반드시 반영

**handoff의 "−20% = Legendre Gibbs, value kernel이 처방" framing은 오진이다.** 증명:
1. value kernel 2버그 수정(descending theta, beta0 정규화) 후 mid/high VZA 폭발 해소, **value≈moment ≤0.3%** (§6.1).
2. nadir −20%는 n_mu(24→48)·m_max(16→32)·L(80→120) 전부 고정 = **kernel-독립·수렴**. φ-Fourier phase moment는 L80서 이미 수렴 → Legendre Gibbs가 −20% 원인 아님.
3. **OCRT single-scatter는 해석적 SS와 thin-AOD서 ≤0.25% 일치**(정확).
4. gap은 **surface/MS 영역**(surface항 6.85e-3 규모와 commensurate), SS·phase·해상도 아님.

→ raw-LUT/value-kernel 전환은 **Legendre 진동-robustness 개선**(g 높은 phase용)이지 −20% 처방이 아니다. 둘을 conflate 금지.

**−20% 확정 다음 실험(Jae 로컬 OSOAA)**: 흡수면(Cox-Munk 제거, 완전흡수)으로 같은 조건 OSOAA nadir 추출 → OCRT black_full(1.318e-2)과 대조. OSOAA black ≈ 1.318e-2면 차이는 순수 surface(Cox-Munk 구현으로 국소화), 높으면 atm_MS 또는 aerosol 자체(ExtData).

---

## 부록: 재발 이슈 빠른 체크 (비교 안 맞을 때)
1. VZA 면 틀렸나? (Case B인데 Snell 안 했나 / Case A인데 했나) → nadir 외 어긋남이 신호.
2. REFL vs I 혼동? (μ_sun=0.866 = 13.4% 오차로 나타남)
3. glint 오염? (wind>0인데 RAA=0 썼나)
4. coupled인데 `--trunc-aer-loglin` 빠졌나? (Tronca 불일치)
5. n-layers/IGmax under-converged? (sweep해서 flat=수렴 확인)
6. depth/chl이 frozen 값 벗어났나? (Case B)
