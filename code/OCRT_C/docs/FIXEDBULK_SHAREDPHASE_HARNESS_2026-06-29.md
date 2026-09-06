# FIXED-BULK SHARED-PHASE 검증 하네스 (2026-06-29)
## HARNESS_OCRT_OSOAA_v1.08_2026-06-27.md 보강 — in-water constituent(Chl/TSM/aDOM) shared-phase 검증
## 매 세션 비교 전 §A–§I 필독. 모든 항목은 2026-06-29 세션에서 실측·코드확인됨.

---

### §A. a_w, b_w (순수해수 IOP) — 공식 사용 절대 금지
**약어**: a_w = 순수해수 흡수계수 [m⁻¹]; b_w = 순수해수 산란계수 [m⁻¹]; bb_w = 후방산란 = 0.5·b_w (Rayleigh 대칭).
- **출처는 파일 하나뿐**: `inputs/water_iop/water_coef_z09_1nm.txt` (OSOAA-exact-matched Z09; a_w=Pope&Fry1997, b_w=OSOAA molecular scattering n=4.29 Morel-like).
- OSOAA도 **동일 테이블** 사용: `OSOAA_PROFILE.F:1877` "OCRT parity default (2026-06-22): use the tabulated OCRT/Z09 table" 패치. → 양측 bit-동일 Z09.
- 6-band 값 (이 값만 쓸 것):

| λ(nm) | a_w (m⁻¹) | b_w (m⁻¹) |
|---|---|---|
| 412 | 4.5506e-3 | 6.650e-3 |
| 443 | 7.069e-3 | 4.872e-3 |
| 490 | 1.500e-2 | 3.164e-3 |
| 555 | 5.960e-2 | 1.859e-3 |
| 660 | 4.100e-1 | 8.875e-4 |
| 865 | 4.605e0 | 2.763e-4 |

- **오류 기록(2026-06-29)**: b_w를 λ⁻⁴·³² 공식(b_w(500)=0.00257)으로 계산하면 파일값 대비 **11–12% 낮음**(412: 0.00593 vs 0.00665; 555: 0.001637 vs 0.001859). a_w는 Pope&Fry라 우연히 일치하나 b_w는 틀림. **반드시 파일을 읽어 쓸 것.**

---

### §B. fixed-bulk = SINGLE MEDIUM (물 ZERO화) — 검증 전제
- `rt_water_rt.c:1638` `a_w = b_w = bb_w = 0.0`. fixed-bulk는 순수해수를 ZERO화하고 주입 total a,b,bb만 단일매질로 사용.
- native(`--simple-chl`)는 §A 파일을 직접 사용하지만 **fixed-bulk는 사용하지 않음** → shared-phase(fixed-bulk)는 total에 §A의 Z09를 직접 더해 주입해야 함.
- **frozen pure-water baseline은 native 경로**라 fixed-bulk의 물 주입을 검증한 적 없음.
- → **fixed-bulk pure-water 체크가 constituent 스윕 전 전제**: chl=0(=물만), §A의 Z09 주입 + Rayleigh phase로 fixed-bulk 실행 → frozen native baseline 재현 확인. (§A b_w 11% 오류처럼 주입 setup 오류가 실재하므로 이 체크가 필수.)

---

### §C. 블렌드 phase (물 Rayleigh + 입자)
**약어**: P_part = 입자 phase; P_Ray = 탈분극 Rayleigh phase; b_part/b_w = 입자/물 산란계수.
- fixed-bulk는 단일 phase라 total b 전체에 적용 → **블렌드 주입**: `(b_part·P_part + b_w·P_Ray)/(b_part+b_w)` (밴드별, chl/Csed별).
- P_Ray (δ_w=0.039, Δ=(1-δ)/(1+δ/2)=0.9426; .mie 정규화 ½∫P11 sinθdθ=1):
  - P11 = Δ·0.75·(1+cos²θ) + (1-Δ)
  - P12 = -Δ·0.75·sin²θ
  - P33 = Δ·1.5·cosθ
- **일관성 검증 필수**: 블렌드 P11 후방반구 적분 = 주입 bb_tot/b_tot. (불일치 시 phase·bb 모순.)
- 물 가중치(b_w/b_total) 무시 불가: chl=0.03에서 9.3%(412)→1.6%(865). 적색/NIR 물흡수(a_w(865)=4.6) 때문에 단일밴드 근사도 금지.

---

### §D. n_mu ≥ 32 (forward-peaked phase 필수) — 정확도 전제
**약어**: n_mu = in-water 이산종좌표(streams) 수; rrs0- = water-side remote-sensing reflectance [sr⁻¹].
- **n_mu=6 절대금지** — forward-peaked(g~0.94-0.97) under-resolve로 다중산란 과대누적 → **6× 과대**.
- chl3/555/sza0 nadir rrs0- 수렴열: n_mu{12,16,24,32,48} = {1.777, 1.670, 1.626, 1.618, 1.611}e-3. 수렴값 ≈1.608e-3.
- n_mu=32 → 0.6% 오차, n_mu=48 → 0.2%. **sub-1%엔 n_mu≥32.**
- **교차검증**: Gordon QSSA rrs0- ≈ G·bb/(a+bb), G=0.0949 → chl3/555 nadir 1.70e-3 (n_mu수렴값과 정합, n_mu=6값과 6× 불일치로 적발).

---

### §E. 실행 = 단일 solve (--output-full-grid 금지) — 속도 전제
- `--output-full-grid`는 출력점(**n_mu vza Gauss노드 × n_raa**)마다 `rt_solve_case_ocean` 재호출 → N배 느려 타임아웃(`main.c:454-467`). 출력 vza는 `--lut-vza-step`이 **아니라** n_mu Gauss노드.
- **단일 런**(`--vza V --raa R`, full-grid 없음)이 빠름 — stdout에 전부:
  - 첫 토큰 = TOA rho_I; 그리고 `rrs0minus=`(in-water), `Rrs0plus=`(above-water), `Lu0minus`, `Ed0minus`, `orders`, `conv`, `a_total/b_total/bb_total`.
- **표준 실행**:
```
OMP_NUM_THREADS=4 ./build/v2_solver_vk --surface ocean --wind-speed 3 \
  --sza S --vza V --raa 90 --wavelength W --pressure 0 --aod 0 \
  --fixed-bulk-iop A B BB --fixed-bulk-phase-lut P11.csv --n-mu-water 32
```
  → `grep -oE 'rrs0minus=[0-9.eE+-]+'`. (atm OFF = --pressure 0 --aod 0, skylight 버그 격리.)

---

### §F. 금지/주의 플래그 (전부 2026-06-29 실측 함정)
- **OCRT_DEBUG=1 금지** — verbose 로깅이 런을 수십배 느리게 해 타임아웃. debug-gated 플래그(`--debug-water-layer-dtau` 등)도 OCRT_DEBUG 필요 → 같이 느려짐. production AUTO(WATER_LAYER n_layers~369 dtau0.05, depth AUTO로 z_max~27.6m τ18.4 truncate, orders~50 conv)가 빠름.
- **--water-mie-phase 금지(forward-peaked)** — vector moment(betal/gammal/alphal/zetal), cap 없음 → lmax200 모먼트로 타임아웃 + Gibbs. 대신 **model 3 = --fixed-bulk-phase-lut**.
- **OMP_NUM_THREADS 미설정 금지** — 전코어 oversubscribe로 느림. =4 권장.
- **--fixed-bulk-phase-lut 포맷**: 헤더 `theta_deg,P11` 필수(csv_find_col로 컬럼명 탐색, log-theta·log-P11 보간). 헤더 없으면 rc=-4 로드실패.
- **--lut-raa-max 무효 플래그** — exit 2. (--lut-vza-step/--lut-vza-max/--lut-raa-step만 유효.)

---

### §G. model 3 cap vs OSOAA truncation (sub-1% 정합)
- model 3 = OSOAA식 smooth log-linear forward cap(T1=6/T2=3, CAP FIX 2026-06-17, delta_f~0.026 weak) + moment-GSF kernel, b*=b(1-A)로 bb보존, **radiance-preserving**.
- OSOAA HYD_TRUNCATION=0(full phase)와 **~few% 차이**(과거 FF bb030 3-way ±3.2% 검증).
- **sub-1% 엄밀비교엔 양측 동일 truncation 필요** — 외부 pre-truncation(delta-M) 후 양측 동일 phase 주입, or OSOAA cap 정합. 먼저 model3 vs OSOAA full로 일치도 보고 >1%면 정합.

---

### §H. mineral(TSM) IOP — Ahn H5.0 mass-specific (4종 모두 준비됨)
**약어**: Csed = 부유광물 농도 [g/m³]; a*/b*/bb* = 질량별(mass-specific) 흡수/산란/후방산란 [m²/g]; bb/b = 후방산란비(밴드무관 상수, Ahn Table 3-II).
- **벡터 phase .mie (프로젝트 4종)**: Brown_earth(g=0.938, 가장 forward-peaked, 모먼트 Lm=80 필요), Yellow_clay, Calcareous_sand, Red_clay.
- **bb/b (Ahn 1999/1990, HydroLight H5.0)**:

| mineral | bb/b |
|---|---|
| Yellow_clay | 0.0092 |
| Brown_earth | 0.0099 |
| Calcareous_sand | 0.0122 |
| Red_clay | 0.0067 |

- **농도 환산**: a_min(λ)=Csed·a*(λ), b_min(λ)=Csed·b*(λ), bb_min(λ)=Csed·(bb/b)·b*(λ).
- 측정 400-700nm; 300-400 & 700-1000nm spline; >1000nm 선형외삽-0충전(비물리, 기계적). 검증밴드 ≤748nm은 측정/spline 범위.
- **Bukata "total suspended mineral" 제외**: a*가 NIR로 증가(비물리), bb/b 없음.
- 4종 전부 milestone #0(v1.07) phase-wiring 검증됨(backscatter Q/I 1-4% 일치). **단 TSM constituent Rrs 검증(matrix #5)은 미완료** — 최근 다성분 작업은 CDOM부터 시작했음.

---

### §I. OSOAA HYD shared-phase (외부 phase 주입)
- `-HYD.Model 2 또는 3` + `-HYD.ExtData <phase파일>` + `-HYD.UserProfile <a,b 프로파일>` (OSOAA_MAIN.F:848).
- ExtData = IMOD=4 포맷 (변환기 `mie_to_osoaa_extdata.py`: col ANGLE/F11/−F12·F11/F22·F11/F33·F11, 0→180 오름차순). 361각이면 `OSOAA.h CTE_MAXNB_ANG_EXT≥400`(현재 700 패치됨).
- OSOAA 런: EXE=/home/user/osoaa_build/exe/OSOAA_MAIN.exe, env OSOAA_NO_DIRECT_GLINT=1. 출력 Standard/RESLUM_vsVZA.txt(REFL air-side), Advanced/RESLUM_Adv_UP.txt(level0=TOA, level27=0- water-side), Flux.txt(Ed).

---

### 비교 규약 (불변, 매번 적용)
- OCRT rho_I = OSOAA REFL(=π·L/Ed), **NOT** OSOAA I(=ρ·μ_sun).
- RAA=90 = glint-free + 규약invariant.
- VZA: OCRT --vza=air-side; 0- 비교는 water-side Snell(air θ → arcsin(sinθ/1.34), n=1.34).
- OCRT_U = −OSOAA_U; nadir Q = |Q| (meridian 특이점).
- 보간 linear/PCHIP만, nearest 금지.
- atm OFF(--pressure 0 / OSOAA AP.MOT 0)로 in-water RT 격리(skylight Ed/Lu 버그 회피).
