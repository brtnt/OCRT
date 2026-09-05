# 세션 시작 정합성 검증 가이드라인 (OCRT v1.08 + OSOAA V2.0)
버전: 2026-07-04 KST | 이 문서는 매 세션(특히 마이그레이션 직후) 가장 먼저 읽고 그대로 수행한다.

## 0. 읽기 순서
1. 본 문서 (검증 절차 + golden 수치)
2. `HARNESS_OCRT_OSOAA_v1.08_2026-06-27.md` §1 규약 · §7 체크리스트 (OSOAA 비교 전 필독)
3. `HARNESS_OCRT_OSOAA_CONSISTENCY_2026-06-30.md` (harness 상세 + trap)
4. 비교 작업 시 `COMPARISON_PITFALLS_AND_PROCEDURE_v1.08.md`

## 1. 재구축 (샌드박스 신규 시)
```sh
cd ocrt && gcc -std=c11 -O2 -fopenmp -Isrc $(find src -name '*.c') -o build/v2_solver_vk -lm
cd osoaa && export OSOAA_ROOT=$PWD && (cd gen && make)   # OSOAA 필요 시
```
env: OCRT는 `OMP_NUM_THREADS=1`(단일코어), 명시적 `--n-mu-water`는 `OCRT_ADVANCED=1` 필요.
OSOAA는 `OSOAA_ROOT` 필수 + `OSOAA_NO_DIRECT_GLINT=1`.

## 2. 정합성 검증 시퀀스 (Tier 순서대로, 각 단계 golden 일치 확인 후 다음 진행)

### Tier 0 — OCRT bit 확인 (~2분)
```sh
OMP_NUM_THREADS=1 ./build/v2_solver_vk --surface ocean --wind-speed 3 --sza 30 --vza 0 \
  --raa 90 --wavelength 555 --pressure 0 --aod 0 \
  --fixed-bulk-iop 0.176200 2.618309 0.038578 \
  --fixed-bulk-phase-lut <Csed5/555 blend P11 CSV>
```
- **golden: `rrs0minus=2.887685e-02` (bit 일치, n_mu=48 기본)**
- blend P11 CSV가 없으면 harness가 Brown_earth.mie로부터 재생성한다(§Tier 3).

### Tier 1 — 가스흡광 (~2분; inputs/xsec 존재 시)
```sh
(위 명령) --wavelength 555 --use-absorption
```
- **golden: `o3: column=9.321e18 mol/cm²(=346.9 DU), τ_abs_col=3.2292e-02`**
- 보조 golden: 602nm τ_O3=4.8015e-02, 440nm τ_NO2=1.7402e-03(0.106 DU)
- "WARNING all 6 xsec files missing" → inputs/xsec 미설치. o3/no2는 continuum 병합본
  (헤더 `CONTINUUM MERGED` 표식)이어야 한다 — LBL-only 판이면 가시광 O3/NO2가 0이다.

### Tier 2 — 옵션 게이트 + LUT 모드 (~1분)
- `--n-mu-water 48` (OCRT_ADVANCED 없이) → **거부되어야 정상**
- `--integration-method constant` (OCRT_DEBUG 없이) → **거부되어야 정상**
- ocean LUT(`--output-full-grid`) 실행 → stderr `n_mu_water default raised to 96` **1건**
- black 표면 LUT → 위 메시지 **0건** (ocean 한정)

### Tier 3 — OCRT↔OSOAA consistency harness (~15–20분, resumable)
```sh
python3 scripts/ocrt_osoaa_consistency_harness.py
```
- **golden: MAPE=1.90% (±5% 밴드 PASS), n=18, max=5.08% @Csed50/865nm, anchor Csed5/555 = −1.15%**
- case별 bit 기준: `golden/golden_ocrt_rrs_nmu48_2026-06-30.csv`
- MINERAL_MIE 기본 경로는 `/mnt/project/Brown_earth.mie` — 새 환경은 `aux/phase_mie/`로
  config 수정 또는 해당 위치 복사.
- 백그라운드 실행은 bash 호출 간 생존하지 않는다 — 포그라운드 + resumable 반복.

### Tier 4 (선택) — n_mu 수렴 spot check
- `scripts/nmu_converge.py 660` (Csed50 고정): n_mu=48 → +0.47%, 64 → +0.14%, 80 → +0.03% (vs 96)

## 3. 알려진 open 결함 — 회귀로 오인 금지
| 결함 | 크기 | 상태 |
|---|---|---|
| skylight Ed/Lu 불일치 (atm ON rrs) | −1.04% | 원인 확정, 수정 대기 |
| aerosol phase moment kernel (glory 누락) | −20% nadir | value kernel 미구현 — **production 에어로졸 LUT 실행 금지** |
| native(simple-chl) delta-M 절단 | +48% 고-ω 412 | TSM 검증 후 수정 예정 |
| red/NIR OCRT<OSOAA (truncation 가설) | −2~−5% | **미검증 가설** — 검증은 n_mu≥80에서 |
| Rayleigh mid-VZA @412 | +2% | 어느 쪽이 옳은지 미확정 |
| 고-ω에서 n_mu=48 미수렴 | +0.44% | 정상(설계) — LUT 모드는 96 기본 |

## 4. 함정 요약
- OSOAA 대기 OFF 불가: AP.MOT 0.01이 하한(1e-4는 전체 NaN)
- OSOAA 0− = Flux.txt level 27; OCRT `--vza`는 air-side(0− off-nadir는 Snell 변환 air=arcsin(sinθ/1.34))
- wind>0은 glint-free RAA(90/180)에서만 비교; OSOAA wind=0 미구현; depth=200 고정
- `OCRT_DEBUG=1`은 생산 런 금지(수십 배 감속). n_mu/n_phi 변경은 `OCRT_ADVANCED=1`
- 샌드박스 shell은 /bin/sh — brace expansion·bashism 금지, python 또는 임시파일 사용
- np.trapz → np.trapezoid (NumPy 2.x)
- 내삽·외삽에 nearest 금지(선형/PCHIP), 외삽은 원칙 금지(clamp + 경고)

## 5. 소실 자산 고지 (찾아 헤매지 말 것)
- 순수해수 native baseline 3종(2026-06-25 동결분) — 소실 확정, 재수립 필요 시 새 정의로
- June-13 last-good 체크포인트 — 부재(−1% ambiguity 판별 불가 지속)
- HydroLight gfortran 포트 — 이 패키지에 없음(로컬 확인)

## 6. 검증 통과 후
전 Tier PASS 시 이 환경은 이전 세션과 정합하다고 판정하고 본 작업을 시작한다.
어느 Tier든 실패 시: **STOP → last-good 대비 diff만 수행** — trap 재조사·재유도 금지,
임의 튜닝 금지(물리 메커니즘 기반 수정만).
