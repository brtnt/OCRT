# OCRT 상세 사용 문서 (v1.09, 2026-07-05)
> **[폐판 안내 2026-07-15]** 본 문서는 `OCRT_USAGE_v1.12.md` 로 대체되었다.
> M1~M4(벡터 기본화·실행 모드 단일화·게이트 확정·도움말 재작성) 이전 기준이므로
> 실행 모드·게이트·`--vector`·`--n-water` 기술이 현행과 다르다. 이력 참조 전용.

본 문서는 v1.09 소스(커밋 #1~#5)와 실검증 커맨드 기준으로 작성됐다.
`--help`(190행)가 1차 레퍼런스이나 aerosol/batch/LUT 옵션군이 누락돼 있어 본 문서가 보완한다.
기호·수식 상세는 소스 주석(문서화 수준, 전 수식 커버)이 원본이다.

## 1. 빌드와 실행 기본
```sh
gcc -std=c11 -O2 -fopenmp -Isrc $(find src -name '*.c') -o build/v2_solver_vk -lm
cd ocrt && OMP_NUM_THREADS=1 ./build/v2_solver_vk [옵션...]
```
- ★함정: solver 는 `inputs/`(water_iop, xsec, afgl_atm)를 **cwd 상대경로**로 연다.
  반드시 `ocrt/`에서 실행. 다른 cwd → `cannot open inputs/...` 전면 실패.
- OMP_NUM_THREADS=1 은 재현성 기준 조건(golden 은 단일스레드 실측).

## 2. 출력 규약
단일 실행 stdout 1행: `rho_I rho_Q rho_U ... 필드=값 ...`
- rho = 반사도 정규화 I/μ_s (F_sun=π 규약, rt_solver.h 참조)
- ocean 케이스는 rrs0minus/Rrs0plus(sr⁻¹), Ed/Lu, Kd/Ku(진단용) 등 부가 필드 출력
- Tier 0 golden: 아래 §7 커맨드에서 `rrs0minus=2.887685e-02` bit 일치가 정상 상태

## 3. 시나리오별 핵심 옵션
### 3.1 대기 (Rayleigh/기체흡수)
- Rayleigh 는 `--pressure`로 제어(0 = off). 모델: `--rayleigh-model bodhaine`(기본,
  v1.09 first-principles) | `6sv`(검증용, DEBUG) | `hansen-travis`(pre-v1.09 재현, DEBUG)
- 기체흡수: `--use-absorption` + `--gas-column-o3 DU` 등. xsec 는 continuum-merged
  판이어야 함(가이드라인 Tier 1 참조)

### 3.2 에어로졸 (--help 누락분, 소스 확정 의미)
- `--mie F` : 6SV형 다중밴드 .mie 로드 (M50C 는 게이트됨 — §5)
- `--aod X` : 실행 파장에서의 AOD
- `--aer-h-km H` : 수직분포 지수형 scale height [km] (0=6SV an23 표, OSOAA 대조 시 2.0)
- `--aer-l-max N` / `--m-max N` : Legendre/Fourier 차수 (OSOAA 대조 프로토콜: 80/16)
- `--trunc-aer-loglin` : **production 절단** (OSOAA 식 log-linear cap, 자체 정합 보상)
- `--aer-phase-kernel value` : 각도공간 value kernel (moment 절단 ringing 회피, production)
- 검증된 OSOAA 대조 프로토콜(MAPE 0.48% 실측): coxmunk wind3 + loglin + value +
  l-max 80 + m-max 16 + n-layers 400 + --vector
- **옵션 등급 (v1.09 커밋 #13, 2026-07-07: 검증 구성 기본값 승격)**
  에어로졸 실행(--aod>0)은 loglin+value+L80+m16+400층이 **자동 기본** —
  최소 커맨드(--mie X --aod Y [--aer-h-km H] --vector)가 곧 검증 프로토콜이다
  (bit 동일 검증 완료). 아래 표는 그 위의 세부 등급:
  | 등급 | 옵션 | 취급 |
  |---|---|---|
  | 프로토콜 상용구 | --trunc-aer-loglin, --aer-phase-kernel value, --vector | 검증 구성 선택 플래그. 기본값 아님 — 에어로졸 검증 실행 시 세트로 복사 |
  | 수렴 노브 | --aer-l-max, --n-layers, --m-max | 자동 기본(80/400/16). 기본값과 다른 값은 OCRT_ADVANCED=1 필요; 같은 값 재명시는 무게이트 |
  | Rayleigh-only | (m=2 등 종전 기본) | 무영향 — 자동승격은 에어로졸 실행 한정 |
  | 전문가 게이트 | --n-mu-water 등 | OCRT_ADVANCED=1 필요 (게이트 일람 §5) |
  | 입력 지정 | --mie, --water-mie-phase, --fixed-bulk-iop | 일반 접근 유지 (게이트 대상 아님) |
  일반 사용: §7 예제를 통째로 복사하고 개별 플래그를 바꾸지 말 것.
- 연구 전용(DEBUG 게이트): `--trunc-aer-fit CUT` `--trunc-aer-m N` `--nt-tau`
  (Wiscombe 표준식으로 v1.09 수정됨; 그래도 절단형상 잔차 ±10% — S-001 기록 참조)

### 3.3 해수 (--surface ocean)
- native 구성성분: `--simple-mode --simple-chl X --simple-min X --simple-adom440 X`
- 검증/디버그 주입: `--fixed-bulk-iop A B BB --fixed-bulk-phase-lut F`
  (LUT CSV: `theta_deg,P11`, **[0,180] 완전 커버 필수** — §5 가드)
- `--wind-speed U` [m/s] : Cox-Munk 거칠기 (0=flat Fresnel)
- `--n-water X` : 해수 굴절률 (기본 1.34; harness/OSOAA 정합값)
- LUT 생산: `--output-full-grid F` (n_mu_water 자동 96 승격, ocean 한정)

### 3.4 batch/LUT
- `--batch IN.csv --output OUT.csv`, `--lut --lut-vza-step S --lut-raa-step S
  --lut-vza-max M [--lut-output DIR]` (`--lut` 는 `--vector` 필수)
- `--q-convention N` / `--sigma-type N` : Stokes Q 부호·Cox-Munk 분산 모델 정수 선택자
  — 값 의미는 소스(rt_types.h) 참조, OSOAA 대조 시 harness 기본값 유지 권장

## 4. 환경변수
| 변수 | 효과 |
|---|---|
| OCRT_DEBUG=1 | 게이트 옵션 해제(§5) + 진단 출력 |
| OCRT_ADVANCED=1 | 수렴 관련 오버라이드 허용(--n-mu-water, --fixed-bulk-phase-nphi) |
| OMP_NUM_THREADS | 재현성 기준 1 |
| OCRT_DUMP_SKY=1 | skylight Ed/Lu 진단 덤프 (결함 조사용, default off) |

## 5. 게이트 일람 (기본 거부 → 사유 출력)
| 대상 | 게이트 | 사유 |
|---|---|---|
| --rayleigh-model 6sv / hansen-travis | DEBUG | 검증/회귀 재현 전용 |
| --integration-method constant | DEBUG | v1 grid-divergence 병리 재현 전용 |
| --n-mu-water, --fixed-bulk-phase-nphi | ADVANCED | 수렴 검증된 기본값 보호 |
| --trunc-aer-fit / -m / --nt-tau | DEBUG | 비검증 연구 절단 (S-001) |
| --sos-acceleration geometric | 항상 거부 | 미구현 dead option (S-003) |
| --mie 경로에 M50C | DEBUG(+경고) | lineage 오염, diagnostic-only (정책 문서) |
| 위상표 [0,180] 미커버 | 항상 거부 | moment 적분기 불일치 (S-002) |
| --debug-water-* 계열 | DEBUG | production 자동선택 보호 |

## 6. 입력 파일 형식
- .mie (6SV형): 헤더(파장별 Ext/Sca/ssa/asym 표) + P11/P12/P33 각도블록(0~180° 필수).
  생성: mie_generator/vrt_solver `--mie-gen X.inp out.mie --dtheta 0.5`
  (★.inp 원본 세트는 현 패키지에 부재 — MANIFEST 참조)
- 위상 LUT CSV: `theta_deg,P11` (+선택 case_id,wavelength_nm)
- AFGL: inputs/afgl_atm/ (userdef 는 afgl_userdef.dat 직접 작성)

## 7. 검증된 예제 커맨드 (전부 본 세션 실측)
```sh
# Tier 0 golden (bit: rrs0minus=2.887685e-02)
OMP_NUM_THREADS=1 ./build/v2_solver_vk --surface ocean --wind-speed 3 --sza 30 --vza 0 \
  --raa 90 --wavelength 555 --pressure 0 --aod 0 \
  --fixed-bulk-iop 0.176200 2.618309 0.038578 --fixed-bulk-phase-lut <blendP11.csv>

# 에어로졸 OSOAA 대조 프로토콜 (실측 −0.3~−1.3% @555/aot0.3/sza40)
OMP_NUM_THREADS=1 ./build/v2_solver_vk --surface coxmunk --wind-speed 3 --sza 40 --vza 30 \
  --raa 90 --wavelength 555 --pressure 0 --mie inputs/M80C.mie --aod 0.3 \
  --aer-l-max 80 --m-max 16 --aer-h-km 2.0 --n-layers 400 \
  --trunc-aer-loglin --aer-phase-kernel value --vector

# 기체흡수 (Tier 1 golden: O3 346.9 DU, τ=3.2292e-02 @555)
... --use-absorption   # 555nm 기본 프로파일

# 대기 LUT (black, --vector 필수)
... --surface black --pressure 1013.25 --vector --lut --lut-vza-step 10 \
    --lut-raa-step 30 --lut-vza-max 60
```

## 8. 주의사항 요약
cwd=ocrt 필수(§1) / 위상표 [0,180] / harness 는 `cd ocrt && python3 scripts/...` /
LUT ocean 은 n_mu_water=96 자동(무겁다: 단일코어 수 분) / 결함·이력은
BUG_SUSPECTS.md 와 가이드라인이 원본 — docs/ 의 v1.08 결함 표는 감사로 일부 철회됨

## --batch-full-grid (v1.09 #19)
행당 1개의 ocean full-grid 를 OpenMP 로 병렬 실행한다(OMP_NUM_THREADS).
CSV 헤더 필수, 빈 칸은 CLI base 상속:
`out`(필수), `wavelength`, `sza`, `wind`, `aod`, `aer_h_km`,
`iop_a`,`iop_b`,`iop_bb`(3종 동반), `water_mie_phase`, `mie`,
`n_mu_water`(OCRT_ADVANCED=1 필요), `raa_step`.
행별 출력 파일이 분리되어 경합이 없다. base CLI 에는 단일 실행과 동일하게
--surface ocean --sza --vza --raa --wavelength 등 필수 인자가 있어야 한다
(base solve 1회 선행 후 폐기). wavelength 를 행에서 바꾸면 n_water 는
자동 재파생된다(--n-water 명시 시 고정).
예:
```
out,wavelength,iop_a,iop_b,iop_bb,water_mie_phase,n_mu_water
out/c03_443.csv,443,0.29,2.37,0.105,phase/case03.mie,24
out/c14_443.csv,443,0.18,4.10,0.093,phase/case14.mie,48
```

## full-grid 직접 추출 (v1.09 #20, 기본 ON)
full-grid 는 모든 vza 를 zero-weight 수중 노드로 추가해 각 (vza,raa) 값을
풀이에서 직접 읽는다(보간 없음, OSOAA UserAngFile 동형). 종전 선형보간으로
되돌리려면(회귀 비교용) `OCRT_GRID_VIEW_INTERP=1`.

## 병렬 실행 성능 노트 (v1.09 #21)
--batch-full-grid 의 solve 경로에서 getenv/할당 경합을 제거했다(env 1회
스냅샷 + 단일 아레나). Windows/MinGW 에서 스레드 스케일링이 회복되는지는
로컬 실측으로 확인할 것. 스냅샷은 프로세스 시작 후 최초 solve 전에 1회
고정되므로, **실행 도중 환경변수 변경은 반영되지 않는다**(원래도 대부분
static 캐시라 동일했음).

## --water-shared-grid (v1.09 #22, OSOAA-parity)
수중 quadrature 가 대기 노드 세트를 공유한다(단일 mu 격자 — OSOAA 방식).
n_mu_water 는 파라미터로서 소멸하며 --n-mu-water / batch n_mu_water 컬럼과
상호 배타다. 각도수는 --n-mu 로 통제(OSOAA parity = 48). 강흡수(NIR) 밴드는
24 로 +0.5% 잔차 실측 — 48 권장. 수중 SOS 는 극소 노드(n_mu<8급) 미지원.
