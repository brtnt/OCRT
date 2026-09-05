# OCRT ↔ OSOAA in-water rrs(0−) 정합성 하네스 — 문서

**버전:** 2026-06-30 KST
**대상:** OCRT v1.08 (v2_solver), OSOAA V2.0 (HYD.Model 3)
**스크립트:** `ocrt_osoaa_consistency_harness.py`

---

## 0. 약어 (표·그림 독립 이해용)

- **OCRT**: KIOST 결합 대기-해양 vector RT 코드 (SOS법)
- **OSOAA**: CNES vector RT 참조코드 (Ordres Successifs Ocean Atmosphere)
- **rrs(0−)**: 해수면 직하 상향 복사반사율 Lu(0−)/Ed(0−) [sr⁻¹]
- **TSM**: 부유광물 (Total Suspended Mineral)
- **Csed**: 부유광물 질량농도 [g/m³]
- **IOP**: 고유광학특성 (a 흡수, b 산란, bb 후방산란) [m⁻¹]
- **a\*, b\***: 질량별 IOP [m²/g] (a=Csed·a\*, b=Csed·b\*)
- **P11**: 위상함수 (1,1) 성분 (스칼라 위상)
- **ExtData**: OSOAA 외부 위상 파일 (IMOD=4 포맷)
- **SOS**: Successive Orders of Scattering (산란차수 순차법)
- **MAPE**: 평균절대백분율오차
- **n_mu**: 물속 Gauss quadrature 노드 수
- **sza/vza/raa**: 태양천정각/관측천정각/상대방위각 [deg]

---

## 1. 목적·범위

세션 시작 시 OCRT 마이그레이션/빌드가 **frozen된 물속 RT 동작을 재현하는지** OSOAA 대비 검증하는 정합성 테스트다. mineral .mie 하나에서 전 입력을 재생성하므로 샌드박스 리셋과 무관하다. 단일 산란차수 검증경로(fixed-bulk)로 RT 엔진 정확성을 격리 확인한다 — production native(simple-chl) 경로 검증과는 별개다.

**Gate = 회귀가드**: golden MAPE 1.90% 대비 materially 나쁘면 FAIL(마이그레이션 실패 의심), 좋으면 IMPROVEMENT(golden 재freeze), 같으면 PASS.

---

## 2. 비교 워크플로

```
mineral.mie ──변환기──▶ OSOAA ExtData (벡터 위상) + bulk ext/sca
    │
    ├─▶ OCRT side : ExtData의 F11(mineral P11) + 내장 Rayleigh P11 블렌드
    │               → fixed-bulk total IOP + 블렌드 위상 LUT
    │               (블렌드 = OCRT↔OSOAA parity 핵심)
    │
    └─▶ OSOAA side: HYD.Model 3 = UserProfile(a,b mineral) + ExtData(위상)
                    + OSOAA 내부 순수해수 Z09
    ↓
  양쪽 실행 (1코어 순차, resumable) → 비교 → MAPE + scatter
```

**블렌드 위상 (양측 parity의 핵심):**
OSOAA는 hydrosol 위상과 순수해수 Rayleigh를 **내부에서** 블렌딩한다. 따라서 OCRT의 단일매질 fixed-bulk에는 같은 블렌드를 주입해야 한다:

```
P_blend(θ) = (b_min·P_min(θ) + b_w·P_Ray(θ)) / (b_min + b_w)
```

- `P_min` = ExtData F11 (mineral 위상)
- `P_Ray` = depolarized Rayleigh, δ=0.039: c₂=(1−δ)/(1+δ), P=(1+c₂cos²θ)/(1+c₂/3). P(0°)=1.4713
- `b_min` = Csed·b\*, `b_w` = 순수해수 산란 (Z09 6밴드 frozen)
- 재구성 검증: 기존 검증 LUT 대비 최대오차 0.0001%

---

## 3. OCRT side 레시피 (fixed-bulk)

**실행:**
```
OMP_NUM_THREADS=1 OCRT_ADVANCED=1 ./build/v2_solver_vk --surface ocean --wind-speed 3 \
  --sza 30 --vza 0 --raa 90 --wavelength <band> --pressure 0 --aod 0 \
  --fixed-bulk-iop <a> <b> <bb> --fixed-bulk-phase-lut <blend_P11.csv> \
  --n-mu-water 48
```
stdout의 `rrs0minus=` grep.

**IOP:** a = Csed·a\* + a_w, b = Csed·b\* + b_w, bb = Csed·(bb/b)·b\* + 0.5·b_w

**Trap (절대 단순화 금지):**
- **옵션 게이트 (2026-06-30)**: 정확도 quadrature(`--n-mu-water`, `--fixed-bulk-phase-nphi`)는 `OCRT_ADVANCED=1` 없이는 거부된다(기본값 48/720이 검증됨). 검증 도구 `--integration-method constant`는 `OCRT_DEBUG=1`이 필요하다. Rayleigh reference-model 및 compatibility-clamp 옵션은 production core에서 제거됐다. **하네스는 OCRT_ADVANCED=1을 설정**하므로 --n-mu-water 48 명시 가능. ocean 표면 shadowing은 강제 on(끄기 옵션 없음).
- **LUT 모드 n_mu 기본 96 (2026-06-30)**: `--output-full-grid`에서는 in-water n_mu_water 기본이 96으로 승격된다(대화형 기본 48은 최고-ω에서 +0.44% 미수렴, n_mu~80 수렴, 96=phyto 전방위상 margin). 명시적 `OCRT_ADVANCED=1 --n-mu-water`가 LUT 모드에서도 우선한다. **LUT 생성 권장 = 96.**
- **fixed-bulk = SINGLE MEDIUM**: rt_water_rt.c:1638이 순수해수를 ZERO화하고 주입 total IOP+위상만 씀. 그래서 블렌드 위상 주입이 필수(물 Rayleigh가 코드 안에 없음).
- **`--fixed-bulk-iop`의 bb는 INERT**: SOS는 위상에 박힌 bb/b를 씀. bb 값은 loader 제약(bb < 0.5·b)만 통과하면 됨.
- **n_mu ≥ 48 권장 (LUT 모드 포함)**: forward-peaked mineral 위상의 floor. **n_mu=32는 수렴값(n_mu=64) 대비 ~1.5% 미수렴**이다 (n_mu=24 ~3%, n_mu=48 ~0.2%; Brown_earth Csed5/555nm). floor는 single-scattering albedo 따라 상승하므로 고-ω 밴드는 >48이 필요할 수 있다 — production LUT는 case별 수렴 확인 필수. (기존 기록의 "n_mu=32 → 0.6%"는 phyto chl3(g=0.97)용이었고, Brown_earth는 다른 phase·ω라 ~1.5%. OCRT 코드 기본값도 24→48로 상향됨 2026-06-30 — default 24에 의존하던 frozen baseline은 재freeze 필요.)
- **`--pressure 0 --aod 0`**: 물속 격리(대기 OFF). 기본 run은 대기 ON임.
- **`--fixed-bulk-phase-lut`(model 3)** = OSOAA식 smooth forward cap. `--water-mie-phase`(벡터 moment, cap 없음)는 lmax200 timeout — 금지. scalar P11만 — nadir rrs엔 충분, 벡터 Q/U는 별도.
- **`OCRT_DEBUG=1` 금지**: verbose 로깅이 수십배 느림.

---

## 4. OSOAA side 레시피 (HYD.Model 3)

**HYD.Model 3** = 사용자 a,b 프로파일(UserProfile) + 외부 위상(ExtData) + OSOAA 내부 순수해수 Z09. Csed/Chl은 무시됨(경고는 정상).

**UserProfile 포맷:** 5 헤더줄 + `(depth_m, a_abs, b_sca)` 행 (mineral만), depth=0 시작·증가순. 균질해양 = 2행 (0,a_min,b_min)+(200,a_min,b_min).

**ExtData:** 변환기 `python3 mie_to_osoaa_extdata.py <mie> <band_nm> <out>` (임의 밴드 PCHIP 보간). 헤더 EXTINCTION_COEF/SCATTERING_COEF + `ANGLE F11 −F12/F11 F22/F11 F33/F11` 361각 0→180. (OSOAA.h CTE_MAXNB_ANG_EXT ≥ 400, 패치됨 700)

**실행 플래그 핵심:** `-HYD.Model 3 -HYD.ExtData <ext> -HYD.UserProfile <prof> -AP.MOT 0.01 -AER.AOTref 0.0 -SEA.Depth 200 -SEA.Ind 1.34 -SEA.Wind 3 -OSOAA.View.Phi 90 -OSOAA.View.Level 1` + env `{OSOAA_ROOT, HOME, OSOAA_NO_DIRECT_GLINT=1}`.

**Trap (절대 단순화 금지):**
- **OSOAA atm OFF 불가**: AP.MOT 1e-4 → **전 복사장 NaN**(TOA 포함). degenerate 대기가 in-water 경계조건을 오염시켜 물속까지 NaN 전파. **AP.MOT 0.01(thin)이 floor.** thin 대기 skylight 효과는 band-independent·무시가능(Csed5/555 anchor 0.06% 일치).
- **OSOAA_ROOT 환경변수 필수**: 누락 시 ERROR_4000.
- **0− 레벨 = Flux.txt의 numeric level 27** (Z=−0.0, Direct_Up=0). Adv_UP의 "0−"는 텍스트 줄("0- is level 27")일 뿐, 데이터 행은 numeric 27.
- **Adv_UP 컬럼**: `level Z vza sca_angle I Q U` → vza=s[2], I=s[4].
- **rrs(0−,nadir) = I(27, vza=0) / Total_Down(27)**. AP.MOT 0.01이면 Total_Down 유효.

---

## 5. Golden 결과 (Brown_earth)

매트릭스: Csed{0.5, 5, 50} × band{412,443,490,555,660,865}, nadir, sza30, raa90, wind3.

**MAPE = 1.90%, max = 5.08%(Csed=50/865nm), n=18 (재수립 2026-06-30, n_mu=48).** Csed별 MAPE: 0.5=2.39%, 5=1.60%, 50=1.71%. anchor Csed=5/555nm = −1.15%.

> **재수립 경위:** 직전 golden(1.60%, max 4.99%)은 n_mu=32(미수렴, 수렴값 대비 +1.5%) 기준이라, "Csed5/555 +0.06%" 앵커가 미수렴과 truncation 차의 우연한 상쇄였다. n_mu=48(수렴, +0.2%)에서 재측정하니 앵커가 −1.15%로 이동하였고, 이는 사전 예측(수렴 시 −1.3%대)과 부합한다. 즉 새 golden 1.90%는 OCRT 수렴오차를 제거한 **실제 OCRT–OSOAA 차이**에 해당한다.

**잔차 구조 (마이그레이션 실패 아님 — 알려진 baseline):**
- 18개 전 케이스에서 OCRT < OSOAA (systematic 음의 편차).
- Blue/green(412–490): −0.3~−2%.
- Red/NIR(660,865): −2~−5%, Csed 클수록 커져 Csed=50/865nm에서 최악(−5.08%).
- **유력 후보 원인(미검증)**: model3-cap(OCRT) vs OSOAA-truncation의 mineral 위상 truncation 불일치. red/NIR은 b_w가 극소(b_w(865)=2.76e-4)라 블렌드가 mineral 위상에 지배되어 truncation 불일치가 노출되고, rrs∝bb/a라 bb 차이가 직결된다고 해석된다. **단 이는 물리적으로 그럴듯한 가설일 뿐 아직 격리·검증되지 않았다. Rule 3에 따라 MAPE를 줄이기 위한 튜닝은 하지 않는다.**
- sub-1% 경로(제안, 미실행): 양측 위상 truncation 정합(외부 pre-truncate 또는 cap delta_f 정합) + 선택적 scalar-equivalent ExtData. 실행 전 red/NIR 편차가 truncation에 기인함을 먼저 격리·검증해야 한다.

golden 수치 파일: `golden_ocrt_rrs_nmu48_2026-06-30.csv` (case별 bit-level reference). 회귀가드는 MAPE 1.90% ±5% 밴드로 PASS/IMPROVEMENT/REGRESSION 판정.

---

## 6. 실행·의존성

**실행:** `python3 ocrt_osoaa_consistency_harness.py` (config 블록 상단 수정으로 mineral/bands/Csed/geometry 조정). resumable — 중단 후 재실행시 완료분 skip. 산출: `consistency_ocrt_vs_osoaa.png` + CSV.

**마이그레이션 후 재빌드 의존성:**
- OCRT: `gcc -std=c11 -O2 -fopenmp -Isrc $(find src -name '*.c') -o build/v2_solver_vk -lm`
- OSOAA: CNES/RadiativeTransferCode-OSOAA + 4 패치파일, `OSOAA.h` CTE_MAXNB_ANG_EXT≥400 후 recompile
- 변환기 `mie_to_osoaa_extdata.py`
- mineral .mie (Brown_earth / Yellow_clay / Calcareous_sand / Red_clay)
- 나머지(ExtData, 블렌드 LUT, Rayleigh, UserProfile)는 하네스가 재생성

**참고:** nproc=1 환경이면 순차 실행(~11분). 병렬화는 코어 추가시에만 의미.

---

## 7. 세션 시작 체크리스트

1. OCRT 빌드 → `build/v2_solver_vk` 존재 확인
2. OSOAA 빌드 → `exe/OSOAA_MAIN.exe` 존재 확인
3. config 블록의 mineral .mie 경로·Z09 a_w/b_w 확인
4. 하네스 실행 → MAPE 산출
5. 회귀가드 판정: golden 1.90% 대비 PASS/IMPROVEMENT/REGRESSION
6. scatter(linear+log+residual) 확인 — red/NIR 잔차 구조가 golden과 일치하는지

---

## 8. 비교 실수 방지 (재발 방지 — 매 세션 확인)

- OCRT --vza = air-side, OSOAA 0− VZA = water-side. nadir(0=0)만 직접비교 가능, off-nadir는 Snell(air θ→arcsin(sinθ/1.34)) 필요.
- OSOAA wind 3/5/10만 유효(wind=0 미구현). RAA=90/180만 glint-free.
- OSOAA depth=200 (5000은 +0.8% 층 artifact). IGmax=1은 −22%(단일산란 아님).
- a_w/b_w는 formula 금지(λ^−4.32는 11% 오차). frozen Z09 6밴드 사용, 다른 밴드는 inputs/water_iop/water_coef_z09_1nm.txt에서 읽기.
