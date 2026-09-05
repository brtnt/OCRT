# OCRT ↔ OSOAA 비교 하네스 문서 (통합 · 권위본)
버전 v1.08 · 2026-06-27 · 위치 `OCRT_AC_LUT_bundle/harness/`
출처 통합: `CONVENTIONS_OCRT_vs_OSOAA.md`(frozen 규약) + harness `README.md`(툴) + `COMPARISON_GUIDELINE`(config·방법론·테이블) + 2026-06-27 세션 findings

---

## 이 문서의 역할 (재발 방지 — 반드시 읽기)

매 세션 시작마다 하네스를 복구하느라 한 세션을 통째로 쓰는 일을 막기 위한 **단일 권위 문서**다. **마이그레이션 시 4종을 동시에 준비한다:**
1. 변경된 **OCRT 코드** (tarball)
2. 변경된 **OSOAA 코드/patch** (tarball)
3. **이 하네스 문서** (마이그레이션 문서와 독립적으로 갱신)
4. **마이그레이션 문서** (체크포인트/복구 절차)

새 세션은 비교 전 이 문서를 읽고, **§1 규약**과 **§7 체크리스트**를 그대로 따른다. 같은 숫자가 같은 기하/물리량을 뜻한다고 **절대 가정 금지**.

---

## 1. 규약 — OCRT와 OSOAA는 정의가 다르다 (frozen, 재유도 금지)

### ★ 1.0 RAA(상대방위각) — 가장 잘 놓치는 함정
- **사실**: OCRT `--raa`와 OSOAA `-OSOAA.View.Phi`는 **정의가 다르다.** RAA=0/180에서 같은 숫자를 직접 비교하면 **틀린 기하를 비교**한다.
- **두 코드가 `90°`에서는 일치(convention-invariant)** — φ↔(180−φ) 관계의 고정점이 90°라서. 이게 **모든 비교를 RAA=90에서 하는** 이유: (a) glint-free, (b) 규약 invariant라 변환 불필요. memory "RAA=90 preferred"의 근거.
- **OCRT 내부(코드 검증)**: 전 경로에 **+π(180°) azimuth shift** — 물 `phi_arg=raa+π`(rt_water_rt.c:139), `phi_v=π−RAA`(:3001), 대기 TOA `V3 raa=(raa_AF+180) mod 360`(rt_solver.c:863). 6SV 검증된 **의도적 규약 — 코드 정상, 버그 아님**.
- **⚠ glint azimuth 정정 (2026-06-27 실측)**: 첨부 `CONVENTIONS` 구판은 "OCRT raa=0 = glint"라 했으나 **실측 결과 반대다.** OCRT TOA `sunglint_up_dir`은 **raa=180, vza=sza(specular)에서 피크** — 즉 **OCRT glint는 raa=180, raa=0은 glint-free.** 코드(rt_solver.c:859 "specular at raa=180")와 일치. 표준이 RAA=90이라 이 오류가 여태 안 물렸을 뿐. **→ 0/180 비교 시 raa=180이 glint측.**
- **non-90 매핑 (출력 레벨별로 다름, 경험적 확인 전제)**: in-water rrs(0−) air-side = OCRT --raa = OSOAA Phi(flip 없음) [v1.06]; above-water Rrs(0+) = 180−Phi [v1.06]; in-water Stokes Q = Phi−180 [v1.07-2]; 대기 coupled = v1.07-5 "identity" 보고(empirical View.Phi=45→Θ=115.86°). 표면상 충돌 = 레벨 의존 → **단일 규칙 단정 금지.**
- **0/180을 꼭 써야 하면 — glint/Θ 검증 먼저**: 라벨 암기 금지. ① OCRT를 wind>0, vza=sza로 raa=0·180 돌려 `sunglint_up_dir` 피크측 = glint측 확인(OCRT는 raa=180). ② OSOAA도 Phi=0·180 동일 시험. ③ 양 코드 산란각 Θ 일치 확인 후 매칭. 확정 전 0/180 raw 비교 금지.
- **역사적 교훈**: 이 규약 무시로 과거 OSOAA가 "OCRT outlier(버그)"로 오판(v1.07-2). 원인은 OCRT측 azimuth 라벨링 오류 — 코드·OSOAA 판정 모두 정상이었음.

### 1.1 VZA(뷰 천정각) — air-side vs water-side
- **OCRT `--vza` = air-side.** 항상.
- **OSOAA**: `View.Level 1`(TOA) = air-side. `View.Level 4`(0−, 수중) = water-side.
- **함정**: 0− 비교 시 OSOAA water-side VZA를 OCRT air-side VZA와 같은 숫자로 비교하면 틀림. nadir(0°)만 우연히 일치, off-nadir 어긋남 = VZA 해석오류 신호.
- **변환(0− 비교 시)**: OCRT air VZA = `arcsin(n_water · sin(OSOAA water VZA))`, **n_water=1.34** (Snell).
- **TOA 비교(aerosol, Level 1)**: 양쪽 air-side → **Snell 불필요**, 같은 VZA 비교 OK.
- **검증**: OSOAA `SCA_ANG`로 확인 — nadir SCA≈158°이면 water-side.

### 1.2 출력 복사량(reflectance) — 정의가 다르다
- **OCRT `rho_I`** = π·L/Ed = 표준 reflectance ρ.
- **OSOAA `REFL` 컬럼** = π·L/Ed = ρ. → **OCRT rho_I ↔ OSOAA REFL** (맞는 짝).
- **함정**: OSOAA `I` 컬럼 = π·L/Esun = ρ·μ_sun. `I`를 `rho_I`와 비교하면 **μ_sun=cos(sza) 배만큼 틀림**(sza30 → 0.866배 = −13.4%). 이 차이는 **전 RAA에 존재**(azimuth invariance와 무관).
- **규칙**: 항상 OSOAA **REFL** 컬럼. `compare_ocrt_osoaa_aer.py`의 `REFL_COL`을 출력 헤더로 확인.

### 1.3 풍속(wind)
- **OSOAA wind=0 미구현.** wind 3/5/10만(Cox-Munk 행렬 사전계산 DATABASE/SURF_MATR: wind{3,5,7,10}×SZA{30,60}). OCRT는 wind=0 가능하나 비교 시 wind≥3.

### 1.4 수심(depth)
- **OSOAA depth=200**(rrs 비교). **depth=5000 금지**(거친 layer 이산화로 rrs +0.8%). aerosol TOA(black water)는 depth=100 충분.

### 1.5 Chl / black water
- **OSOAA Chl=0.0**(순수). **0.00001 금지**(+0.13% artifact).
- **black water(aerosol-only TOA)**: `-YS.Abs440 1000` + `-PHYTO.Chl 0 -SED.Csed 0 -DET.Abs440 0`.
- OSOAA chl=1 0− with atmosphere는 신뢰 불가(1.82× hard bound).

### 1.6 다중산란 차수(IGmax)
- **`-SOS.IGmax 100`(또는 unlimited) = full MS.** `IGmax=1`은 single-scatter가 **아님**(−22% deficit) — single-scatter 목적이어도 쓰지 말 것.

### 1.7 OSOAA hydrosol truncation
- 비교 기준선 **`OSOAA_HYD_TRUNCATION=0`**(안전). 켜면 별도 정합.

### 1.8 Rayleigh — bw는 이미 일치
- OSOAA patch로 **Z09(Zhang 2009) bw 항상 사용** → OCRT와 Rayleigh 산란계수 동일. → 잔여 Rayleigh @412 mid-VZA +2% 차이는 **bw 아님, depolarization factor + Rayleigh phase matrix** 차이(어느 쪽 옳은지 미해결). Rayleigh 디버깅은 여기부터. aerosol-only(`--pressure 0`/`AP.Pressure 0`)엔 무관.

### 1.9 ExtData / 변환기
- OSOAA AER.ExtData(IMOD=4): 컬럼 `ANGLE F11 −F12/F11 F22/F11 F33/F11`, col3 = **−F12/F11**(부호반전), 각도 **0→180 오름차순**(.mie는 180→0 내림차순이라 역순). F22=F11=1(구형). 변환기 `harness/mie_to_osoaa_extdata.py`(검증됨)가 처리. 361각도 → `inc/OSOAA.h` CTE_MAXNB_ANG_EXT≥400.

---

## 2. 하네스 툴 3종 (`OCRT_AC_LUT_bundle/harness/`)

| 툴 | 역할 | 사용 |
|---|---|---|
| `mie_to_osoaa_extdata.py` | .mie → OSOAA AER.ExtData (양 코드 동일 phase) | `python … <file.mie> <wl_nm> <out.extdata>` |
| `mie_phase_consistency.py` | 층 A gate (Legendre 진동 감지, 풀 RT 불필요) | `python … <file.mie> <wl_nm> [L…]` (L에 production 80 포함) |
| `compare_ocrt_osoaa_aer.py` | 층 B gate (복사휘도 비교 + scatter) | `python … <ocrt_value.csv> <ocrt_moment.csv> <osoaa_LUM_vsVZA.txt> [golden_d%]` |

> **하네스 복구 우선순위**: 위 3개 스크립트가 없으면 비교 불가. 마이그레이션 tarball에 **반드시 포함**(과거 rev3가 "포함"이라 적고 누락 → 한 세션 낭비). 없으면 README §툴 스펙 + `mie_io.c` 파서 기준 재구성.

---

## 3. Matched configuration + 실행 명령 (frozen)

### 3.1 Case A: aerosol-only TOA (현 active)
**OCRT** (value/moment 각각):
```
./build/v2_solver_vk --surface coxmunk --wind-speed 3 --sza 30 --vza 0 --raa 90 \
  --wavelength 412 --pressure 0 --mie inputs/M80C.mie --aod 0.2 \
  --aer-l-max 80 --m-max 16 --aer-h-km 2.0 --n-layers 400 --trunc-aer-loglin \
  --aer-phase-kernel value --vector --lut --lut-vza-step 10 --lut-vza-max 80 \
  --lut-raa-step 30 --output-full-grid /tmp/vk_value.csv
```
**필수 3플래그**(coupled 정합): `--aer-h-km 2.0`(OCRT 기본 an23 V=23km ↔ OSOAA AP.HA 2.0), `--n-layers 400`(기본 40 부족), `--trunc-aer-loglin`(coupled은 자동 안 켜짐, main.c:1349는 `aerosol_on && !rayleigh_on`만 → 명시해야 OSOAA Tronca 1과 일치). **aerosol-only는 `--pressure 0`**(coupled@412는 Rayleigh가 nadir 지배 → aerosol 결함 가림).

**OSOAA**:
```
export OSOAA_ROOT=$(pwd)
python harness/mie_to_osoaa_extdata.py M80C.mie 412 /tmp/M80C_412.extdata
exe/OSOAA_MAIN.exe -OSOAA.Wa 0.412 -OSOAA.View.Level 1 -OSOAA.View.Phi 90 -SOS.IGmax 100 \
  -AER.Model 4 -AER.ExtData /tmp/M80C_412.extdata -AER.Tronca 1 -AER.Waref 0.412 -AER.AOTref 0.2 \
  -AP.HA 2.0 -AP.HR 8.0 -AP.Pressure 0.0 \
  -SEA.Dir $OSOAA_ROOT/DATABASE/SURF_MATR -SEA.Ind 1.34 -SEA.Wind 3 -SEA.Depth 100 \
  -YS.Abs440 1000 -PHYTO.Chl 0 -SED.Csed 0 -DET.Abs440 0 \
  -OSOAA.Log run.log -OSOAA.ResRoot /tmp/OSOAA_AERONLY_M80C_412
```
coupled은 양쪽 `Pressure 1013.25`. `-OSOAA.Log`는 relative만(absolute → ERROR_900).

### 3.2 Case B: pure-water 0− (frozen baseline / regression guard)
- 555nm, wind 3, SZA 30, RAA 90, **depth 200**, **chl 0.0**, VZA air-side Snell.
- **OSOAA rrs(0−) = Adv_UP level27 의 I / Flux level27 의 Ed** (출력 profile 0− 레벨 행; `View.Level 4` 인자와 별개).
- IGmax 무제한. **IOP parity 확인**: ω 0.02728, b_w 1.67e-3, a_w 5.96e-2.
- baseline 산출물(FROZEN 2026-06-25): `OCRT_purewater_atm_BASELINE_golden_v1.08_2026-06-25.csv`, `ocrt_purewater_atm_REGRESSION_GUARD.py`, `OCRT_BASELINE_RULES_purewater_atm_2026-06-25.md`.

---

## 4. 2-layer 검증 방법론

과거 "풀 RT 복사휘도만" 비교는 **coupled@412서 Rayleigh가 aerosol backscatter 결함을 가려** 통과했다(과거 "aerosol 0.43% 검증"이 결함 못 봄). 두 층으로 막는다.

**층 A — phase-consistency (풀 RT 불필요, 신규 aerosol 작업 전 필수 gate)**
raw P11을 OCRT 규약 Σ(2l+1)βₗPₗ로 재구성 → backscatter(Θ≥120°) 음수·진동·비수렴 감지. exit 0=PASS,1=FAIL. **moment-only면 FAIL이 정상**(결함 노출). 예 M80C@412: L40 min=−0.375, L80 min=−0.0199(음수=Gibbs), GATE FAIL. ⚠ native μ sparse라 고차 βₗ 정량 부정확(정성 음수 판정은 견고) — dense-θ PCHIP 내삽 βₗ 재계산이 TODO.

**층 B — 복사휘도 (풀 RT, 반드시 aerosol-only)**
VZA별 REFL·moment/value·d%·MAPE·nadir gate·1:1 scatter(linear). ⚠ `VZA_COL=0,REFL_COL=3`은 가정 → 실제 헤더 확인.

**워크플로우**: ①층 A → ②수정 후 층 A2(effective phase, TODO) → ③층 B 전 VZA 수렴+scatter → ④gate 통과=체크포인트(golden 수치파일+스크립트가 유일한 비교, prose 수치는 baseline 아님).

---

## 5. OSOAA 측 셋업 (patches)
pristine OSOAA V2.0 대비 수정, 재컴파일 필수:

| Patch | 파일 | 내용 | 거동 |
|---|---|---|---|
| P1 | inc/OSOAA.h | CTE_MAXNB_ANG_EXT 200→**400** | 배열 상한(361각도 수용). 물리 불변 |
| P2 | OSOAA_TRPHI.F | env `OSOAA_NO_DIRECT_GLINT` glint 제외 | **default-OFF**. 진단용. RAA=90 무관 |
| P3a | OSOAA_PROFILE.F | Z09 Rayleigh bw 항상-on | **물리 변경**(§1.8). aerosol-only 무관 |
| P3b | OSOAA_PROFILE.F | user-depth J-1 선형 보간 | depth 정합 |
| P4 | OSOAA_MAIN.F | seawater index 기본 1.34 | `-SEA.Ind` 명시 시 무관 |

⚠ sandbox 빌드 불가(현 패키지엔 수정 4파일+harness만, OSOAA_SOS.F/SURFACE.F/Makefile/DATABASE 본체 없음) → **OSOAA 실행은 Jae 로컬**. surface는 precomputed Cox-Munk BRDF(`DATABASE/SURF_MATR`).

---

## 6. 결과 테이블 (template + 현 reference)

> 단일 golden 결과테이블 없음. harness가 run별 생성. 아래는 현 reference + 세션 측정. **OSOAA 참조값 재확인 필요**(§8).

### 6.1 VZA sweep (Case A: M80C 412nm aerosol-only, sza30 raa90 aod0.2, n-layers400)
약어: d%=100·(OCRT/OSOAA−1).

| VZA | OSOAA REFL | OCRT moment | d% | OCRT value | d% |
|---:|---:|---:|---:|---:|---:|
| 0 | 2.5061e-2 | 2.0033e-2 | −20.1 | 2.0087e-2 | −19.8 |
| 10 | 2.2669e-2 | 1.9681e-2 | −13.2 | 1.9642e-2 | −13.4 |
| 20 | 1.8888e-2 | 1.8157e-2 | −3.9 | 1.8158e-2 | −3.9 |
| 40 | 1.3054e-2 | 1.2956e-2 | −0.8 | 1.2957e-2 | −0.7 |
| 60 | 1.8181e-2 | 1.8015e-2 | −0.9 | 1.8014e-2 | −0.9 |
| 80 | 5.1645e-2 | 5.0331e-2 | −2.5 | 5.0337e-2 | −2.5 |
| MAPE | | | 6.89% | | 6.88% |

### 6.2 nadir 분해 (Case A AOD=0.2) — OSOAA-독립
약어: SS=single scattering(해석적, OCRT 일치 확정), atm_MS=atmosphere MS, surface=Cox-Munk Fresnel+결합 MS.

| 항 | rho_I | cox_full 대비 |
|---|---:|---:|
| SS (해석적, 정확) | 1.0789e-2 | 53.9% |
| atm_MS | 2.3908e-3 | 11.9% |
| surface | 6.8531e-3 | 34.2% |
| **OCRT cox_full** | **2.0033e-2** | 100% |
| OSOAA ref | 2.5061e-2 | — |
| gap | 5.03e-3 | OSOAA +25.1% |

### 6.3 권장 컬럼 (신규 비교)
`case | particle | wl_nm | sza | raa | wind | aod | vza | OSOAA_REFL | OCRT_moment | OCRT_value | d%_moment | d%_value` + MAPE + nadir gate.

---

## 7. 비교 전 체크리스트 (매번)
1. **RAA=90** 쓰는가? (아니면 §1.0 glint/Θ 검증 — glint는 OCRT raa=180측)
2. **VZA**: TOA면 air-side 그대로; 0−면 Snell(n=1.34) 변환했는가?
3. OSOAA **REFL** 컬럼인가? (I 아님 — I면 μ_sun 배 틀림)
4. wind≥3, depth(0−는 200), Chl=0.0?
5. IGmax≥100, HYD_TRUNCATION=0?
6. aerosol coupled면 3플래그(`--aer-h-km 2.0 --n-layers 400 --trunc-aer-loglin`)? aerosol 검증이면 aerosol-only(`--pressure 0`)?
7. 층 A gate 통과(또는 moment-only FAIL이 예상대로)?

---

## 8. Known issues / 2026-06-27 세션 findings (반드시 반영)

1. **value kernel 2버그 수정 완료**(rt_aerosol_runtime.c): ROOT CAUSE 1 descending theta(오름차순 저장), ROOT CAUSE 2 beta0 정규화(P11/P12/P33을 ½∫P11 dμ로 나눔, 4096-grid PCHIP bit-consistent). → mid/high VZA 폭발 해소, value≈moment ≤0.3%.
2. **−20% nadir gap은 Legendre Gibbs 아님**(핸드오프 framing 오진). 증명: value≈moment, nadir −20%가 n_mu·m_max·L 전부 고정(수렴·kernel-독립), **OCRT single-scatter는 해석적 SS와 thin-AOD서 ≤0.25% 일치**. gap은 **surface/MS 영역**(surface항 6.85e-3 규모). raw-LUT/value-kernel은 진동-robustness 개선이지 −20% 처방 아님.
3. **−20% 다음 실험(Jae 로컬 OSOAA)**: 흡수면으로 OSOAA nadir 추출 → OCRT black_full(1.318e-2)과 대조. ≈1.318e-2면 surface 문제, 높으면 atm_MS/aerosol.
4. **glint azimuth 정정**: OCRT glint는 **raa=180**(실측), 구 CONVENTIONS 문서의 "raa=0"은 오류(§1.0).
5. **phase LUT interpolation**: 기록 표준 = P11 PCHIP log10(θ)–log10(P11)(offline .mie 구축), P12/P33 linear PCHIP. value kernel runtime은 현재 linear → moment와 일관성 위해 linear-θ PCHIP 정렬 권장(저위험, 2차).
6. **미해결**: Rayleigh @412 depolarization/phase matrix 차이(§1.8); LUT driver `ocrt_ac_lut.py` 3플래그 미반영; native simple-chl delta-M 절단 버그(412nm rrs +48%).
