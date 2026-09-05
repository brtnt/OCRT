# OCRT 상세 사용 문서 (v1.12, 2026-07-15)
본 문서는 v1.12 소스(M1~M4 + W-PCHIP + S7 계열 수정 반영) 기준으로 작성됐다.
v1.09판(2026-07-05)을 대체한다. `--help`(일반 옵션)와 `--examples`(검증 예시 6종)가
1차 레퍼런스이고, 본 문서는 게이트 옵션·배치·입력 형식 등 심화 내용을 다룬다.
기호·수식 상세는 소스 주석(문서화 수준, 전 수식 커버)이 원본이다.

유지보수 규칙: 옵션/게이트/골든 값을 바꾸면 main.c 게이트 일람 주석, `--help`,
`--examples`, 본 문서를 **같은 커밋에서** 갱신한다. `--examples`의 예시를 바꾸면
반드시 실측 재검증한다.

v1.09 → v1.12 주요 변경 요약(상세는 CHECKPOINT_2026-07-14.md):
- M1: 스칼라 솔버 삭제. **벡터 I/Q/U 가 항상 기본**.
- M2: 실행 모드 3종으로 단일화(§2). 충돌 플래그 조합은 오류 종료.
- M3: 수치 노브 ADVANCED 게이트 + 진단 옵션 DEBUG 게이트 확정(§6).
- M4/2026-07-15: `--help` 는 일반 옵션 전용, 게이트 옵션은 미수록(본 문서 §6 이 목록).
  `--examples` 신설(검증된 실행 예시 6종 출력).
- 2026-07-15 옵션 정리 1차: `--lut-output` 삭제(`--output-full-grid` 와 완전 중복),
  `--output-mode {simple|debug}` 삭제 → `--output-advanced` 플래그로 대체
  (기본=간단 열: I/Q/U 반사도만; 플래그 지정 시 `--batch` CSV 에 기준 대비 오차·
  수렴 차수·계산 시간 열 추가), `--view-as-node` 삭제(단일 기하 모드가 자동
  활성하므로 순수 선언용이었음).  삭제된 옵션 지정 시 미지 인자 오류로 종료.
- 2026-07-15 옵션 정리 6차(G5 1차 물성·CDOM 확정): `--f-sun` **삭제**(TOA
  조도는 내부 π 규약 -- 결합 모델에서 대기가 항상 선행 정의: pressure 기본
  1013.25 + 기체흡수 기본 활성), `--water-aw-lut`/`--water-psi-T-lut` **삭제**
  (inputs/water_iop/ 고정), **CDOM 3종 세트 규칙** 신설. adom 대조 스크립트에
  ref-lambda 440 보충(물리 불변).
- 2026-07-15 옵션 정리 5차(G4 해면·계면 확정): `--surface` **필수화**,
  `--wind-speed` ocean/black_fresnel_ocean **필수화**(0 명시 허용), `--n-water` 기본
  **1.34 고정**(Quan-Fry 자동 폐지), `--sigma-type`/`--q-convention` →
  **ADVANCED**, 선글린트 기본 **포함**으로 반전(`--decouple-sunglint` 신설,
  구 옵션 삭제). 골든 러너·Tier-0 프로토콜은 동결 물리 명시 핀으로 bit 보존
  (골든 2·3·5·6 각 5/5 재통과 실증).
- 2026-07-15 옵션 정리 4차(G3 에어로졸 확정): `--no-aerosol`(무동작)
  `--aer-phase-kernel` `--trunc-aer-loglin` **삭제**(자동 dispatch 일원화 +
  결합 케이스 확장). `--aer-h-km` 기본 0→**2.0**. 가드 신설(aod>0 은 mie 필수,
  aod=0 에서 에어로졸 옵션 오류). 러너·하니스의 loglin 토큰 제거 후 골든
  2·3·5·6 비트 재통과.
- 2026-07-15 옵션 정리 3차(G2 대기 확정): `--use-absorption`(기본 활성 전환)
  `--tau-r` `--tau-r-from-input` `--dump-atm` `--rayleigh-model` `--xsec-dir`
  `--afgl-dir` **삭제**. 러너 5종·OSOAA 대조 5종·하니스 2종·LUT 생산 1종에 6기체 0 주입
  (골든 러너 2·3·5·6 비트 재통과로 물리 보존 실증).
- 2026-07-15 옵션 정리 2차(1안 확정): `--vector` `--lut` **완전 삭제**. 패키지
  내부 스크립트(골든 러너 5, OSOAA 대조 5, 하니스 2, LUT 생산 1, ps1 도구 2)의
  해당 토큰을 같은 커밋에서 제거했고, 골든 러너 비트 재통과로 결과 불변을
  실증했다. ★패키지 밖 개인 스크립트에 두 토큰이 남아 있으면 미지 인자
  오류로 멈추므로 직접 제거해야 한다. 아카이브 사본(`*_ORIGINAL_*`)은 이력
  보존을 위해 수정하지 않았다.
- W-PCHIP/S7: 격자 산출과 단일 기하가 rrs(0-)·Rrs(0+) 7유효자리로 일치.

## 1. 빌드와 실행 기본
```sh
cd ocrt && gcc -std=c11 -O3 -march=native -ffp-contract=fast -fassociative-math \
  -fno-signed-zeros -fno-trapping-math -DOCRT_FAST_KERNELS -fopenmp -Isrc \
  $(find src -name '*.c') -o build/v2_solver_vk -lm
```
- ★함정: solver 는 `inputs/`(water_iop, xsec, afgl_atm)를 **cwd 상대경로**로 연다.
  반드시 `ocrt/`에서 실행. 다른 cwd → `cannot open inputs/...` 전면 실패.
- OMP_NUM_THREADS=1 은 재현성·속도 측정 기준 조건(golden 은 단일스레드 실측).

## 2. 실행 모드 (M2, 2026-07-14 단일화)
모드는 세 가지뿐이며, 모드 정책이 충돌하는 플래그 조합은 오류로 즉시 종료한다.

| 모드 | 선택 조건 | 동작 |
|---|---|---|
| 단일 기하 | `--vza`+`--raa` **둘 다** 지정 | 그 한 방향만 계산. 관측 방향을 0-가중 절점으로 자동 추가(view-as-node)해 보간 없이 추출 |
| 격자(기본) | `--vza`/`--raa` 미지정 | 균일 vza/raa 격자(기본 2.5°, vza 0..85°)를 단일 solve 로 산출. 간격은 `--lut-vza-step`/`--lut-raa-step`(해양-대기 공통) |
| 배치 | `--batch` 또는 `--batch-full-grid` | CSV 행별 연속/병렬 실행(§10). `--batch` 결과 상세 열은 `--output-advanced` |

규칙(전부 실측 확인, 2026-07-15):
- `--vza` 하나만 지정 = 오류. `--raa` 도 동일.
- 단일 기하 플래그 + 격자/배치 플래그 혼용 = 오류
  (`conflicting mode flags: ... Pick one.`).
- `--sza` 와 `--wavelength` 는 **모든 모드에서 필수**(배치에서도 행 미지정 시의
  기본값 역할로 커맨드라인에 있어야 함).
- 격자 출력 파일 미지정 시 실행 디렉터리에 `result<N>.csv` 자동 명명(N 최대+1).
  확장자 없으면 `.csv` 자동 부여, `.csv` 외 확장자는 오류.
- `--lut` 와 `--vector` 는 2026-07-15 삭제되었다(문서 서두 변경 요약 참조).
  격자 모드는 기본이므로 별도 선언이 필요 없다.

## 3. 출력 규약
단일 기하 stdout 1행: `rho_I rho_Q rho_U ... 필드=값 ...`
- rho = 반사도 정규화 I/μ_s (F_sun=π 규약, rt_solver.h 참조). 항상 벡터(I/Q/U).
- ocean 케이스는 rrs0minus/Rrs0plus(sr⁻¹), Ed/Lu, Kd/Ku(진단용) 등 부가 필드 출력.
- 격자 모드 CSV: 헤더 1행 + (vza,raa) 행. ocean 은 Rrs(0+) I/Q/U 와 rrs(0-) 수록.
- Tier 0 golden: §8 커맨드에서 `rrs0minus=2.887685e-02` bit 일치가 정상 상태.

## 4. 시나리오별 핵심 옵션
### 4.1 대기 (Rayleigh/기체흡수)
- Rayleigh 는 `--pressure`로 제어(0 = off, 기본 1013.25). 모델은 Bodhaine 1999
  first-principles **고정**이다(2026-07-15: `--rayleigh-model` 삭제, hansen-travis
  구현 제거 -- 이로써 τ_R 기준 미결 사안은 Bodhaine 으로 종결. 레거시 H-T 델타
  수치는 v1.09 changelog 에 보존). `--tau-r`/`--tau-r-from-input`/`--dump-atm` 도
  삭제 -- τ_R 은 사용자가 신경 쓰지 않는다.
- 기체흡수는 **기본 활성**이다(2026-07-15 확정; 구 `--use-absorption` 삭제).
  HITRAN 단면적(continuum-merged, `inputs/xsec` 고정 -- `--xsec-dir` 삭제) +
  AFGL 프로파일 층별 적분. 총량 지정(`--gas-column-o3 DU` 등)은 형상 유지·총량
  스케일. **완전 off = 6기체 총량 전부 0**(레일리 전용 결과와 bit 동일, 실측).
  압력 0 이어도 기체흡수는 독립 적용되며, 비율량(rrs/Rrs/Kd)은 불변이라
  Tier-0 앵커는 유지된다(실측). OSOAA 는 기체흡수가 없으므로 **모든 OCRT↔OSOAA
  대조는 6기체 0 필수** -- 대조 스크립트·골든 러너 전 호출에 주입 완료.

### 4.2 에어로졸 (2026-07-15 G3 확정)
- **on/off 는 `--aod` 하나다**: 기본 0 = 에어로졸 없음이며, 이때 에어로졸
  옵션을 하나라도 지정하면 오류다. `--aod > 0` 에는 `--mie` 가 필수다(가드
  실재; 배치는 행 단위로 검사·SKIP). 구 `--no-aerosol` 은 무동작이라 삭제.
- `--mie F` : 6SV형 다중밴드 .mie 로드 (M50C 는 DEBUG 게이트 — §6).
  ★.inp 직접 입력 인터페이스는 없다 — mie_generator 2단계
  (`--mie-gen X.inp out.mie`)로 .mie 를 먼저 만든다(§7).
- `--aer-h-km H` : 수직분포 scale height [km]. **기본 2.0(지수형)**;
  0 = 6SV an23 표 프로파일 선택지.
- **최소 커맨드가 곧 검증 프로토콜이다**: `--aod>0` 이면 loglin 절단 + value
  커널 + aer-l-max 80 + m-max 16 + n-layers 400 + scale height 2.0 이 자동
  기본이다. 2026-07-15 부로 **결합(레일리+에어로졸) 케이스까지 자동**이다
  (종전에는 에어로졸 단독만 자동이라 결합 최소 명령이 loglin 없이 돌았다 —
  명시 프로토콜 대비 I 0.51% @555 실측). 신최소 명령이 구 명시 프로토콜과
  bit 동일함을 실측 검증했다. `--aer-phase-kernel`/`--trunc-aer-loglin` 은
  자동화에 따라 삭제(moment 커널은 CLI 도달 불가로 전환).
- 수렴 노브(--aer-l-max/--n-layers/--m-max)는 ADVANCED 게이트(기본값 재명시는
  무게이트).
- 연구 전용(DEBUG 게이트): `--trunc-aer-fit CUT` `--trunc-aer-m N` `--nt-tau`
  출력 U 부호는 OCRT meridian-basis convention으로 고정된다.

### 4.3 해수 (--surface ocean)
- CDOM 은 **3종 세트**다(2026-07-15): `--cdom-a440 --cdom-slope
  --cdom-ref-lambda` 를 셋 다 지정(활성)하거나 전부 미지정(비활성)한다.
  부분 지정은 오류(숨은 기본값 혼입 방지).
- native 구성성분(SIMPLE 모드): `--simple-mode --simple-chl X --simple-min X
  --simple-adom440 X`. `--simple-*` 하위 옵션만으로도 모드는 활성화되지만,
  **`--surface ocean`은 반드시 명시**한다. SIMPLE은 수중 구성성분→IOP/위상 변환만
  선택하며 대기 aerosol, wind, 선글린트 정책을 바꾸지 않는다. native `--cdom-*`,
  `--fixed-bulk-*`, `--water-mie-*`와 혼용하면 오류다. `--ccrr-*` 별칭은 경고 후
  호환 수용한다. ★제약: `--simple-chl > 0` 은 `aph_bricaud_1998.txt` 부재로 현
  패키지에서 **실행 실패**한다(이월 항목 — 2026-07-15 실측 재확인). 엽록소
  케이스는 고정 IOP 주입 경로(아래)를 쓴다(골든 항목3 방식).
- 검증/디버그 주입: `--fixed-bulk-iop A B BB --fixed-bulk-phase-lut F`
  (LUT CSV: `theta_deg,P11`, **[0,180] 완전 커버 필수** — §6 가드)
- water-Mie 위상 경로: `--water-mie-truncation`은 S-009 정책 확정 전 opt-in,
  `--water-mie-ss-mode 0|1|2`는 각각 보정 없음/IMS/NT-TMS다. 구현된 public CLI이며
  `--water-mie-moment-nmu` 기본은 400이다.
- `--wind-speed U` [m/s] : 거친 Fresnel 해면의 경사분산 (0=flat Fresnel).
  **ocean/black_fresnel_ocean 에서 필수다(2026-07-15; 0 도 명시, 배치는 행 단위 검사).**
- `--n-water X` : 해수 굴절률. **기본 1.34 고정**(2026-07-15; Quan-Fry(λ) 자동
  기본 폐지). 종전 자동 기본은 OSOAA(1.34 고정) 대조에서 굴절률 불일치를
  만들고 있었다 -- 새 기본이 대조 정합이다. 골든 러너는 동결 물리를
  밴드별 Quan-Fry 값으로 명시 핀했다.
- 해양 격자: `--output-full-grid F` (n_mu_water 자동 96 승격, 단일코어에서 무겁다:
  15°/30° 격자 72행 실측 약 2분, 기본 2.5° 는 훨씬 무겁다)

### 4.4 계면 규약 (2026-07-15 G4 확정)
- `--surface` 는 **필수**다(black|flat|black_fresnel_ocean|ocean; 조용한 black 기본 폐지,
  ccrr-compare 모드만 예외). lambert 는 미구현 명시 차단.
- `--q-convention`(0=legacy (Rs−Rp)/2, 1=Hansen/Mishchenko=기본)과
  `--sigma-model`(nakajima-tanaka 또는 ocrt-floor=기본; 구 `--sigma-type` 숫자 별칭 유지)은 **ADVANCED
  게이트**다(기본 1 재명시는 통과).
- 선글린트: **기본은 직달 단일반사 항 포함**이다(2026-07-15 반전).
  `--decouple-sunglint` 지정 시 표적 출력에서 사후 분리·제외한다(AF1982
  기준이 분리형이므로 대조 시 지정). 구 `--no-decouple-sunglint` 삭제.

## 5. 환경변수
| 변수 | 효과 |
|---|---|
| OCRT_DEBUG=1 | 게이트 옵션 해제(§6) + 진단 출력 |
| OCRT_ADVANCED=1 | 수렴 파라미터 오버라이드 허용(§6 ADVANCED 행) |
| OMP_NUM_THREADS | 재현성 기준 1. --batch-full-grid 병렬 폭 |
| OCRT_DUMP_SKY=1 | skylight Ed/Lu 진단 덤프 (결함 조사용, default off) |
| OCRT_GRID_VIEW_INTERP=1 | full-grid 를 구식 선형보간 추출로 회귀(회귀 비교 전용, §11) |

## 6. 게이트 일람 (M3 확정판; 기본 거부 → 사유 출력)
`--help` 에는 아래 옵션이 실리지 않는다(2026-07-15 방침). 본 표가 목록이다.

| 대상 | 게이트 | 사유 |
|---|---|---|
| --n-mu --n-layers --max-orders --sos-max-orders --sos-tolerance --conv-tol --l-max --m-max --water-m-max --n-mu-water --fixed-bulk-phase-nphi + 에어로졸 노브(--aer-l-max 등) | ADVANCED (기본값 재명시는 통과) | 수렴 검증된 기본값 보호 |
| --integration-method constant | DEBUG | v1 grid-divergence 병리 재현 전용 |
| --trunc-aer-fit / -m / --nt-tau | DEBUG | 비검증 연구 절단 (S-001) |
| --scalar-ff-* / --ccrr-particle-phase-* / --simple-compare | DEBUG | 진단 경로 |
| --sos-save-orders | DEBUG | 산란차수 진단 전용 |
| --debug-water-* 계열 / --debug-fixed-bulk-* | DEBUG | production 자동선택 보호 |
| --mie 경로에 M50C | DEBUG(+경고) | lineage 오염, diagnostic-only (정책 문서) |
| --sos-acceleration geometric | 항상 거부 | 미구현 dead option (S-003) |
| 위상표 [0,180] 미커버 | 항상 거부 | moment 적분기 불일치 (S-002) |

## 7. 입력 파일 형식
- .mie (6SV형): 헤더(파장별 Ext/Sca/ssa/asym 표) + P11/P12/P33 각도블록(0~180° 필수).
  생성: mie_generator/vrt_solver `--mie-gen X.inp out.mie --dtheta 0.5`
  (★.inp 원본 세트는 현 패키지에 부재 — MANIFEST 참조)
- 위상 LUT CSV: `theta_deg,P11` (+선택 case_id,wavelength_nm)
- AFGL: inputs/afgl_atm/ **경로 고정**(`--afgl-dir`/`--xsec-dir` 삭제, 2026-07-15).
  userdef 는 inputs/afgl_atm/afgl_userdef.dat 에 직접 작성

## 8. 검증된 예제 커맨드 (전부 2026-07-15 v1.12 빌드 실측)
일상 예시 6종(해양 단일/격자, 대기, 에어로졸, 기체흡수, 배치)은
`./build/v2_solver_vk --examples` 가 소요 시간과 함께 출력한다. 아래는 그 외
검증·프로토콜용이다.
```sh
# Tier 0 golden (bit: rrs0minus=2.887685e-02; --n-mu-water 명시 없음, 기본 48)
# 2026-07-15 G4 핀: 동결(07-12) 당시 물리 명시 -- Quan-Fry n(555)=1.341224933502,
# 선글린트 분리.  핀 없이 신규 기본(1.34, 포함)으로 돌리면 앵커가
# 2.887871e-02 로 이동한다(+0.0064%; 차기 재동결 시 채택 여부 결정 항목).
OCRT_ADVANCED=1 OMP_NUM_THREADS=1 ./build/v2_solver_vk --surface ocean --wind-speed 3 \
  --sza 30 --vza 0 --raa 90 --wavelength 555 --pressure 0 --aod 0 \
  --n-water 1.341224933502 --decouple-sunglint \
  --fixed-bulk-iop 0.176200 2.618309 0.038578 \
  --fixed-bulk-phase-lut ../validation/assets/phase_lut/tier0_cs5.0_555.csv

# 에어로졸 OSOAA 대조 프로토콜 — 최소 커맨드가 곧 프로토콜(자동 기본, §4.2)
OMP_NUM_THREADS=1 ./build/v2_solver_vk --surface black_fresnel_ocean --wind-speed 3 --sza 40 \
  --vza 30 --raa 90 --wavelength 555 --pressure 1013.25 \
  --mie inputs/M80C.mie --aod 0.3 --aer-h-km 2.0

# 기체흡수 제어 (기본 활성; 총량 오버라이드 예 / off 는 6기체 0 -- --examples [5])
OMP_NUM_THREADS=1 ./build/v2_solver_vk --surface black_fresnel_ocean --wind-speed 3 --sza 40 \
  --vza 30 --raa 90 --wavelength 555 --pressure 1013.25 --gas-column-o3 200

# 대기 격자 (black 표면)
OMP_NUM_THREADS=1 ./build/v2_solver_vk --surface black --pressure 1013.25 --sza 40 \
  --wavelength 555 --lut-vza-step 30 --lut-raa-step 90 --lut-vza-max 60 \
  --output-full-grid atm_grid.csv
```

## 9. 주의사항 요약
cwd=ocrt 필수(§1) / 위상표 [0,180] / `--sza`+`--wavelength` 전 모드 필수 /
모드 플래그 혼용 금지(§2) / LUT ocean 은 n_mu_water=96 자동(무겁다: 단일코어
수 분) / native chl 경로는 aph 파일 부재로 불가(§4.3) / harness 는
`cd ocrt && python3 scripts/...` / 결함·이력은 BUG_SUSPECTS.md 와 가이드라인이
원본 — docs/ 의 v1.08 결함 표는 감사로 일부 철회됨

## 10. --batch-full-grid (v1.09 #19; M2 갱신)
행당 1개의 ocean full-grid 를 OpenMP 로 병렬 실행한다(OMP_NUM_THREADS).
CSV 헤더 필수, 빈 칸은 CLI base 상속:
`out`(필수), `wavelength`, `sza`, `wind`, `aod`, `aer_h_km`,
`iop_a`,`iop_b`,`iop_bb`(3종 동반), `water_mie_phase`, `mie`,
`n_mu_water`(OCRT_ADVANCED=1 필요), `raa_step`.
행별 출력 파일이 분리되어 경합이 없다.

★M2 변경(2026-07-15 실측): base CLI 에 `--vza`/`--raa` 를 넣으면 모드 충돌
**오류**다(v1.09 문서의 "필수" 안내는 폐기). 필수는 `--surface ocean` 과
`--sza`/`--wavelength`(행 미지정 시 기본값 역할) 뿐이다.
wavelength 를 행에서 바꾸면 n_water 는 자동 재파생된다(--n-water 명시 시 고정).
예:
```
out,wavelength,iop_a,iop_b,iop_bb,water_mie_phase,n_mu_water
out/c03_443.csv,443,0.29,2.37,0.105,phase/case03.mie,24
out/c14_443.csv,443,0.18,4.10,0.093,phase/case14.mie,48
```

## 11. full-grid 뷰 추출 (v1.09 #20 + v1.12 갱신)
full-grid 는 모든 vza 를 zero-weight 노드로 추가해 각 (vza,raa) 값을 풀이에서
직접 읽는다(보간 없음, OSOAA UserAngFile 동형). v1.12 에서 S7-relax/S7-boa
수정으로 격자 산출과 단일 기하가 rrs(0-)·Rrs(0+) 7유효자리로 일치한다
(2026-07-14 게이트). 종전 선형보간으로 되돌리려면(회귀 비교용)
`OCRT_GRID_VIEW_INTERP=1`.

## 12. 병렬 실행 성능 노트 (v1.09 #21 존치)
--batch-full-grid 의 solve 경로에서 getenv/할당 경합을 제거했다(env 1회
스냅샷 + 단일 아레나). 스냅샷은 프로세스 시작 후 최초 solve 전에 1회
고정되므로, **실행 도중 환경변수 변경은 반영되지 않는다**.

## 13. --water-shared-grid (v1.09 #22 존치, OSOAA-parity)
수중 quadrature 가 대기 노드 세트를 공유한다(단일 mu 격자 — OSOAA 방식).
n_mu_water 는 파라미터로서 소멸하며 --n-mu-water / batch n_mu_water 컬럼과
상호 배타다. 각도수는 --n-mu 로 통제(OSOAA parity = 48). 강흡수(NIR) 밴드는
24 로 +0.5% 잔차 실측 — 48 권장. 수중 SOS 는 극소 노드(n_mu<8급) 미지원.
