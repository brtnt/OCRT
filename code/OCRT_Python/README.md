# OCRT 파이썬 배치 자료생산 패키지 (GOCI-III 편광 논문)

C `--water-model ocrt` 경로를 파이썬으로 재현하고, GPU(CuPy)/CPU(numpy)
배치로 grid 전체 자료를 생산한다. 대기 SOS·수중 SOS·공기-물 결합·구성모델·
가스흡수·에어로졸 value 커널·Rrs까지 배치화되어 있으며, 단일 케이스 대비
비트 일치(또는 7유효자리) 수준으로 검증되었다.

## 구성

- `ocrt_py/`         : 복사전달 파이썬 모듈. 배치 엔진은 `atmos_batch.py`,
                       오케스트레이터는 `batch_driver.py`.
- `data/`            : 입력 자료.
    - 에어로졸 `.mie` 96종 = OPAC 16종(C50/T50/M80C/M95C/M98C/O99 등)
      + AccuRT 유래 Ahmad 80종(`r<RH>f<FINE>v01`).
    - 순수수/식물플랑크톤/무기물 IOP, AFGL 대기, 기체 흡수 단면적.
- `produce_grid.py`  : **배치 자료생산 진입점.** grid CSV → 12밴드 × R1/R2/R3.
- `gpu_gate.py`      : GPU 무결성 게이트(numpy vs CuPy 7유효자리 대조).
- `ocrt_solve.py`    : 단일 케이스/grid CLI(검증·소규모용, 기존).
- `verify.py`        : C 대비 검증 재확인.
- 실행 스크립트      : `run_produce.bat/.sh`(CPU), `run_produce_gpu.bat/.sh`(GPU),
                       `run_gpu_gate.bat`(게이트).
- `example_grid.csv` : grid CSV 형식 예시.

## 요구 환경

- CPU 실행: Python 3.8+ 와 NumPy. `pip install numpy`
- GPU 실행: 위에 더해 CuPy. CUDA 버전에 맞게 설치.
  예: `pip install cupy-cuda12x` (CUDA 12.x) 또는 `cupy-cuda11x`.

## grid CSV 형식

열: `case_id, sza, vza, raa, wind, aerosol, aod865, chl, tsm, acdom440`
(각도는 도, aod865는 865 nm 기준 AOD, chl은 mg/m³, tsm은 g/m³, acdom440은 /m).
`example_grid.csv` 참조. 밴드별 AOD는 `.mie` 스펙트럴 소광비로 자동 스케일한다.

## 포함된 grid 및 재생성

- `full_grid_design_v1.csv` : 1000케이스 설계 grid(자료생산 입력).
- `iop_grid_design_v1.csv`  : 해수 IOP만 담은 부속(검증·기록용).
- `make_grid.py` / `make_grid.bat` / `make_grid.sh` : grid 재생성 스크립트.

설계 사양: 해수 IOP(TSM 0.1-20, Chl 0.1-5, aCDOM440 0.01-0.1)는 로그공간 상관
다변량 정규 + 범위 truncation으로 대양→연안 covary(목표 상관 TSM-Chl 0.85,
Chl-aCDOM 0.75, TSM-aCDOM 0.77). 에어로졸 6종 균등, AOD865 0.05-0.3 균등,
풍속 1-10 균등, 기하는 sunglint-free(Cox-Munk 직접 glint rho_g < 5e-4).

```bash
# grid 재생성 (시드 고정 = 재현 가능)
python make_grid.py --n 1000 --seed 20260718
# 또는 make_grid.bat / ./make_grid.sh
```

**주의**: 포함된 `full_grid_design_v1.csv`는 원본 설계 스크립트/시드가 소실되어
사양(범위·상관·sunglint 기준)으로 재생성한 것이다. truncation 때문에 실현 상관은
목표보다 약간 낮다(0.85→0.81, 0.75→0.70, 0.77→0.72). 원본과 개별 케이스는
다르나 설계 사양은 동일하다.


## 산출 반사도 (행/케이스마다 세 실행, 차감 분리)

- R1 `rho_TOA`        : ocean + 레일리 대기(기압 1013.25) + 기체흡수
                        + 에어로졸(파장별 스케일) + 구성모델 water-leaving.
- R2 `rho_R`          : Cox-Munk 레일리 대기(에어로졸 없음), 선글린트 decouple.
- R3 `rho_(R+A,black)`: Cox-Munk 레일리+에어로졸(black ocean), 선글린트 decouple.
- 도출: `rho_RC = R1-R2`, `rho_A+rho_RA = R3-R2`, `t*rho_w = R1-R3`,
        `Rrs = R1의 Rrs(0+)`.  항등 `rho_RC = (rho_A+rho_RA) + t*rho_w` 성립.

출력 CSV 열: `case_id, band_nm, rho_TOA, rho_R, rho_RpA, rho_RC, rho_A_RA,
t_rho_w, Rrs`.

## 사용법

### CPU (numpy)

```bash
# 리눅스/맥
./run_produce.sh full_grid_design_v1.csv result.csv

# Windows
run_produce.bat full_grid_design_v1.csv result.csv

# 직접 호출(옵션 조정)
python produce_grid.py --grid full_grid_design_v1.csv --out result.csv \
    --data data --n-mu-water 24 --nt-atm 400 --m-max 16 --max-it-water 500
```

### GPU (CuPy)

```bash
# 리눅스/맥
./run_produce_gpu.sh full_grid_design_v1.csv result.csv

# Windows
run_produce_gpu.bat full_grid_design_v1.csv result.csv
```

내부적으로 `OCRT_PY_GPU=1`을 설정한다. 각 행이 독립이므로 배치 축(B)으로
스택해 GPU에서 병렬 처리한다.

### 재개·밴드 선택

- 재개: 출력 CSV에 이미 있는 `(case_id, band)`는 자동 건너뛴다. 같은 `--out`을
  다시 지정하면 중단 지점부터 이어서 실행한다.
- 밴드 선택: `--bands 443,555,865` (기본 12밴드 전체).
- 가스흡수 끄기: `--no-gas` (기본 켜짐).

## GPU 무결성 게이트 (로컬 필수)

배치의 부동소수 연산 순서가 CPU와 GPU에서 다를 수 있어, 결과가 **비트 일치가
아니라 7유효자리(상대오차 1e-7)** 로 일치하는지 게이트로 확인한다. 이 게이트는
CuPy가 설치된 로컬 CUDA에서만 실행 가능하다.

```bash
# 소규모 grid(4~16행)로 CPU와 GPU 결과를 대조
python gpu_gate.py --grid small_grid.csv --data data --bands 443,555,865
# 또는 Windows: run_gpu_gate.bat small_grid.csv
```

`rho_TOA/rho_R/rho_RpA/Rrs` 각각의 최대 상대오차와 PASS/FAIL, CPU/GPU 속도를
출력한다. 자료생산 전에 이 게이트를 통과시키는 것을 권장한다.

## 옵션 기본값 (자료생산 조건)

- `--n-mu-water 24`   : 수중 GL 절점.
- `--m-max 16`        : 대기 SOS 푸리에 모드 수.
- `--nt-atm 400`      : 대기 SOS 층 수.
- `--max-it-water 500`: **수중 산란차수 상한. 혼탁수(고산란) 밴드는 320차수
  이상 필요하므로 기본 160으로는 수렴 실패한다. 500 권장(기본값).**
- `--chunk 64`        : 한 배치의 행 수. GPU 메모리에 맞춰 조정.
- `--phyto`를 생략하면 Stage-2 검증 기본 Chl 흡수 경로를 사용한다.
  흡수는 `data/water_iop/phyto_absorption_default.csv`, 식물플랑크톤
  `b=bb=0`, 입자 위상은 detritus/mineral만 사용한다.
- Chl>0에서 `--phyto`를 명시하면 종 이름과 무관하게 fail-loud 한다.
  17종 EAP catalog는 생성기/향후 검증용이며 constituent 생산 RT에 자동 투입되지 않는다.

## 검증 요약 (배치 vs 단일)

- R2(레일리+black Fresnel ocean surface): rel 0 (완전 비트 일치).
- R3(레일리+에어로졸, value 커널): rel 0~2e-16 (C50/M95C/Ahmad).
- R1(3패스 결합 ρ_TOA): rho_I rel ~5e-16 (수렴 케이스), 7유효자리 이내.
- Rrs(0+): 배치 vs 단일 rel ~2e-16 (Lu_0plus·Ed_0plus 모두 일치).
- 가스흡수 ON: R3 rel 5e-16, R1 rel 0.
- 항등층 패딩(층수 다른 케이스 한 배치): rel ~2e-16.

## 성능 주의

파이썬 CPU 결합 1행은 자료생산 조건에서 케이스당 수십 초~수 분이다. 1000행 ×
12밴드 = 12000 케이스는 CPU로는 오래 걸리므로 **GPU 실행을 권장**한다. GPU
게이트를 먼저 통과시켜 무결성을 확인한 뒤 전체 자료를 생산한다.

## 에어로졸 모델 주의

`data/`의 Ahmad 80종(`r<RH>f<FINE>v01`)은 **AccuRT 유래 Ahmad 계열**이다
(AccuRT 건조 굴절률·미세물리 + Ahmad 2010 RH 상태·미세모드 분율·성장인자).
Ahmad et al. (2010) 논문 원본 Table-4 입자 반지름의 문자 그대로의 재현이
아니다. 논문 원본 버전을 별도로 쓰려면 그에 해당하는 `.mie` 파일을 추가해야
한다(현재 패키지에는 AccuRT 유래 버전만 포함).

## 2026-07-21 full-grid aerosol object fix

Version: `pyOCRT-v1.2-2026-08-09-chl-zero-extension-legacy-mie-cleanup`

The atmosphere-ocean angular full-grid path is available through:

```bash
python ocrt_fullgrid.py --out grid.csv --wl 490 --sza 30 --wind 3 \
  --chl 1.2 --tsm 3.5 --acdom440 0.04 --aod865 0.2 --aer C50 \
  --vza-values 0,30 --raa-values 0,90,180,270 --no-gas
```

The `.mie` file is read once and a frozen `AerosolRuntime` is prepared once before the geometry loop. The same object is explicitly passed to every cell and both atmospheric passes. Supplying `AOD > 0` without aerosol data now raises an error instead of silently using a Rayleigh-only/default state.

Detailed results are in `validation/python_fullgrid_aerosol_object_fix/VALIDATION_REPORT.md`.


## 2026-07-21 black Fresnel ocean surface naming

The canonical high-level Python surface mode is `black_fresnel_ocean`. It means
a rough Fresnel interface over a black ocean, with no water-leaving radiance.
The legacy `coxmunk` mode and function names remain deprecated aliases. The
default slope variance is documented separately as the OCRT floor law,
`sigma^2 = 0.003 + 0.00512 * max(0.01, W)`, not as an exact Cox-Munk model.


## EAP phase data update (2026-07-26)

The obsolete `pico_Synechococcus_EAP.mie`, `nano_Haptophytes_EAP.mie`, and `Diatoms_centric_EAP.mie` files were removed. Seventeen regenerated representative EAP files are under `data/water_iop/eap/`. The legacy `pico`, `nano`, and `micro` names remain compatibility aliases to the new Synechococcus, Prymnesiaceae and centric-diatom files. The original source data and reproducible C generator are included under `data/eap_source/` and `tools/eap_phase_generator/`.


## EAP 식물플랑크톤 Stage-2 계약

Chl>0에서 종 옵션을 생략하면 검증된 기본 흡수 스펙트럼을 사용하고
식물플랑크톤 산란은 0으로 둔다. `--phyto`/`--phyto-group`을 명시하면
`OCRT_ADVANCED=1` 여부와 무관하게 fail-loud 한다. 17종 EAP phase 파일과
생성기는 패키지에 유지되지만 현재 constituent 생산 RT에는 연결하지 않는다.


## Optional direct water value kernels (2026-07-27)

These paths are opt-in and do not change the default coefficient kernel.

- `OCRT_VALUE_PHASE_SPLINE=1`: scalar direct value kernel using cubic-Hermite-style spline interpolation in `mu=cos(theta)` and Gauss normalization.
- `OCRT_WATER_VALUE_KERNEL_POL=1`: direct polarized value kernel (`pfm, gr, gt, arr, art, att`) when raw `value_phase` data are supplied by the caller.
- `OCRT_WATER_MIE_MOMENT_N_MU=400`: normalization quadrature count.
- `OCRT_WATER_VALUE_NPHI=N`: explicit azimuth integration count.

Exact VZA=0 is a known polarized value-kernel output-extraction singularity. Use VZA=0.001 deg (validated safe range 0.001--0.01 deg); do not use angles below about 1e-6 deg as a nadir workaround. No automatic geometry substitution is performed.

The OCRT constituent Chl gate is extended through 1100 nm. The phytoplankton absorption table explicitly covers 350–1100 nm: source-backed values end at 850 nm and the approved 850–1100 nm extension is exactly zero. Pure-water a_w and total b_w are now tabulated from 200 through 2449 nm. The separate temperature-correction LUT remains limited to 400–900 nm, so non-20 °C calculations outside that interval require a future psi_T extension.
