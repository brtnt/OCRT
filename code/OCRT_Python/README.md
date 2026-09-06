# OCRT Python batch data-production package (GOCI-III polarization study)

This package reproduces the C `--water-model ocrt` path in Python and produces full-grid data
in GPU (CuPy) or CPU (NumPy) batches. The atmospheric SOS, the in-water SOS, the air–water
coupling, the constituent model, gas absorption, the aerosol value kernel and Rrs are all
batched. Batch results agree with single-case results to the bit (or to 7 significant digits).

## Contents

- `ocrt_py/`: radiative-transfer Python modules. The batch engine is `atmos_batch.py`; the
  orchestrator is `batch_driver.py`.
- `data/`: input data.
    - aerosol `.mie` tables (96 models used by this package = 16 OPAC models such as
      C50/T50/M80C/M95C/M98C/O99 + 80 AccuRT-derived Ahmad models `r<RH>f<FINE>v01`);
      the tables are installed from the GitHub Release `data-v1` (see the root `README.md`).
    - pure-water / phytoplankton / mineral IOPs, AFGL atmospheres, gas cross sections.
- `produce_grid.py`: **entry point of batch data production.** grid CSV → 12 bands × R1/R2/R3.
- `gpu_gate.py`: GPU integrity gate (NumPy vs CuPy, 7 significant digits).
- `ocrt_solve.py`: single-case / grid CLI (validation and small runs).
- `verify.py`: re-check against the C code.
- run scripts: `run_produce.bat/.sh` (CPU), `run_produce_gpu.bat/.sh` (GPU),
  `run_gpu_gate.bat` (gate).
- `example_grid.csv`: example of the grid CSV format.

## Requirements

- CPU runs: Python 3.8+ and NumPy (`pip install numpy`).
- GPU runs: in addition CuPy, matching the CUDA version, e.g. `pip install cupy-cuda12x`
  (CUDA 12.x) or `cupy-cuda11x`.

## Grid CSV format

Columns: `case_id, sza, vza, raa, wind, aerosol, aod865, chl, tsm, acdom440`
(angles in degrees; aod865 is the AOD at 865 nm; chl in mg/m³; tsm in g/m³; acdom440 in 1/m).
See `example_grid.csv`. The AOD of each band is scaled automatically with the spectral extinction
ratio of the `.mie` table.

## Included grids and regeneration

- `full_grid_design_v1.csv`: 1000-case design grid (input of data production).
- `iop_grid_design_v1.csv`: companion file with the water IOPs only (validation and record).
- `make_grid.py` / `make_grid.bat` / `make_grid.sh`: grid regeneration scripts.

Design: the water IOPs (TSM 0.1–20, Chl 0.1–5, aCDOM440 0.01–0.1) follow a correlated
multivariate normal distribution in log space with range truncation, so that open-ocean and
coastal cases co-vary (target correlations TSM–Chl 0.85, Chl–aCDOM 0.75, TSM–aCDOM 0.77). Six
aerosol models uniform, AOD865 0.05–0.3 uniform, wind 1–10 uniform, and sun-glint-free geometry
(Cox–Munk direct glint rho_g < 5e-4).

```bash
# regenerate the grid (fixed seed = reproducible)
python make_grid.py --n 1000 --seed 20260718
# or make_grid.bat / ./make_grid.sh
```

**Note.** The included `full_grid_design_v1.csv` was regenerated from the design specification
(ranges, correlations, glint criterion) because the original script and seed were lost. Because
of the truncation the realized correlations are a little below the targets (0.85→0.81,
0.75→0.70, 0.77→0.72). Individual cases differ from the original, but the design is the same.

## Output reflectances (three runs per row/case, separated by subtraction)

- R1 `rho_TOA`: ocean + Rayleigh atmosphere (1013.25 hPa) + gas absorption + aerosol (scaled per
  wavelength) + constituent-model water-leaving radiance.
- R2 `rho_R`: Cox–Munk Rayleigh atmosphere (no aerosol), sun glint decoupled.
- R3 `rho_(R+A,black)`: Cox–Munk Rayleigh + aerosol over a black ocean, sun glint decoupled.
- Derived: `rho_RC = R1 − R2`, `rho_A + rho_RA = R3 − R2`, `t·rho_w = R1 − R3`,
  `Rrs = Rrs(0+) of R1`. The identity `rho_RC = (rho_A + rho_RA) + t·rho_w` holds.

Output CSV columns: `case_id, band_nm, rho_TOA, rho_R, rho_RpA, rho_RC, rho_A_RA, t_rho_w, Rrs`.

## Usage

### CPU (NumPy)

```bash
# Linux / macOS
./run_produce.sh full_grid_design_v1.csv result.csv

# Windows
run_produce.bat full_grid_design_v1.csv result.csv

# direct call (adjust options)
python produce_grid.py --grid full_grid_design_v1.csv --out result.csv \
    --data data --n-mu-water 24 --nt-atm 400 --m-max 16 --max-it-water 500
```

### GPU (CuPy)

```bash
# Linux / macOS
./run_produce_gpu.sh full_grid_design_v1.csv result.csv

# Windows
run_produce_gpu.bat full_grid_design_v1.csv result.csv
```

Internally `OCRT_PY_GPU=1` is set. Rows are independent, so they are stacked along a batch axis
(B) and processed in parallel on the GPU.

### Resume and band selection

- Resume: `(case_id, band)` pairs already present in the output CSV are skipped. Running again
  with the same `--out` continues from the point of interruption.
- Band selection: `--bands 443,555,865` (default: all 12 bands).
- Gas absorption off: `--no-gas` (default on).

## GPU integrity gate (required locally)

Floating-point operation order can differ between CPU and GPU. The gate therefore checks that
the results agree **to 7 significant digits (relative error 1e-7)**, not to the bit. It runs
only on a local CUDA machine with CuPy.

```bash
# compare CPU and GPU on a small grid (4–16 rows)
python gpu_gate.py --grid small_grid.csv --data data --bands 443,555,865
# or on Windows: run_gpu_gate.bat small_grid.csv
```

It prints the maximum relative error of `rho_TOA/rho_R/rho_RpA/Rrs`, PASS/FAIL and the CPU/GPU
speed. Passing this gate before data production is recommended.

## Default options (data-production conditions)

- `--n-mu-water 24`: in-water Gauss–Legendre nodes.
- `--m-max 16`: number of Fourier modes of the atmospheric SOS.
- `--nt-atm 400`: number of atmospheric SOS layers.
- `--max-it-water 500`: **upper limit of in-water scattering orders. Turbid (strongly scattering)
  bands need more than 320 orders, so the old default 160 fails to converge. 500 is recommended
  (default).**
- `--chunk 64`: rows per batch. Adjust to the GPU memory.
- Without `--phyto`, the Stage-2 validated default Chl absorption path is used: absorption from
  `data/water_iop/phyto_absorption_default.csv`, phytoplankton `b = bb = 0`, particle phase from
  detritus / minerals only.
- With Chl > 0 and an explicit `--phyto`, the run fails loudly regardless of the species name.
  The 17-species EAP catalog is for the generator and for future validation; it is not fed into
  the constituent production RT automatically.

## Validation summary (batch vs single case)

- R2 (Rayleigh + black Fresnel ocean surface): relative difference 0 (bit-identical).
- R3 (Rayleigh + aerosol, value kernel): 0 to 2e-16 (C50 / M95C / Ahmad).
- R1 (three-pass coupled ρ_TOA): rho_I ~5e-16 (converged cases), within 7 significant digits.
- Rrs(0+): ~2e-16 (both Lu_0plus and Ed_0plus agree).
- Gas absorption ON: R3 5e-16, R1 0.
- Identity-layer padding (cases with different layer counts in one batch): ~2e-16.

## Performance note

One coupled row on the Python CPU path takes tens of seconds to minutes per case under
data-production conditions. 1000 rows × 12 bands = 12 000 cases take long on the CPU, so
**GPU execution is recommended**. Pass the GPU gate first, then produce the full data set.

## Aerosol model note

The 80 Ahmad models in `data/` (`r<RH>f<FINE>v01`) are the **AccuRT-derived Ahmad family**
(AccuRT dry refractive index and microphysics + Ahmad 2010 RH states, fine-mode fractions and
growth factors). They are not a literal reproduction of the Table-4 particle radii of Ahmad et
al. (2010). To use the paper-based version, add the corresponding `.mie` files (this package
contains the AccuRT-derived version only).

## 2026-07-21 full-grid aerosol object fix

Version: `pyOCRT-v1.2-2026-08-09-chl-zero-extension-legacy-mie-cleanup`

The atmosphere–ocean angular full-grid path is available through:

```bash
python ocrt_fullgrid.py --out grid.csv --wl 490 --sza 30 --wind 3 \
  --chl 1.2 --tsm 3.5 --acdom440 0.04 --aod865 0.2 --aer C50 \
  --vza-values 0,30 --raa-values 0,90,180,270 --no-gas
```

The `.mie` file is read once and a frozen `AerosolRuntime` is prepared once before the geometry
loop. The same object is passed explicitly to every cell and to both atmospheric passes. Giving
`AOD > 0` without aerosol data now raises an error instead of silently using a Rayleigh-only or
default state.

Details: `validation/python_fullgrid_aerosol_object_fix/VALIDATION_REPORT.md`.

## 2026-07-21 black Fresnel ocean surface naming

The canonical high-level Python surface mode is `black_fresnel_ocean`: a rough Fresnel interface
over a black ocean, with no water-leaving radiance. The legacy `coxmunk` mode and function names
remain as deprecated aliases. The default slope variance is documented separately as the OCRT
floor law, `sigma² = 0.003 + 0.00512 · max(0.01, W)`, not as an exact Cox–Munk model.

## EAP phase data update (2026-07-26)

The obsolete `pico_Synechococcus_EAP.mie`, `nano_Haptophytes_EAP.mie` and
`Diatoms_centric_EAP.mie` files were removed. Seventeen regenerated representative EAP files are
under `data/water_iop/eap/`. The legacy names `pico`, `nano` and `micro` remain as compatibility
aliases of the new Synechococcus, Prymnesiaceae and centric-diatom files. The original source data
and the reproducible C generator are under `data/eap_source/` and `tools/eap_phase_generator/`.

## EAP phytoplankton Stage-2 contract

With Chl > 0 and no species option, the validated default absorption spectrum is used and
phytoplankton scattering is set to zero. An explicit `--phyto` / `--phyto-group` fails loudly
regardless of `OCRT_ADVANCED=1`. The 17-species EAP phase files and the generator stay in the
package but are not connected to the constituent production RT.

## Optional direct water value kernels (2026-07-27)

These paths are opt-in and do not change the default coefficient kernel.

- `OCRT_VALUE_PHASE_SPLINE=1`: scalar direct value kernel with cubic-Hermite-style spline
  interpolation in `mu = cos(theta)` and Gauss normalization.
- `OCRT_WATER_VALUE_KERNEL_POL=1`: direct polarized value kernel (`pfm, gr, gt, arr, art, att`)
  when raw `value_phase` data are supplied by the caller.
- `OCRT_WATER_MIE_MOMENT_N_MU=400`: normalization quadrature count.
- `OCRT_WATER_VALUE_NPHI=N`: explicit azimuth integration count.

Exactly VZA = 0 is a known output-extraction singularity of the polarized value kernel. Use
VZA = 0.001° (validated safe range 0.001–0.01°); do not use angles below about 1e-6° as a nadir
workaround. No automatic geometry substitution is performed.

The OCRT constituent Chl gate is extended to 1100 nm. The phytoplankton absorption table covers
350–1100 nm explicitly: source-backed values end at 850 nm and the approved 850–1100 nm
extension is exactly zero. Pure-water a_w and total b_w are tabulated from 200 to 2449 nm. The
separate temperature-correction table is still limited to 400–900 nm, so calculations at
temperatures other than 20 °C outside that interval need a future psi_T extension.

---

# OCRT Python 배치 자료생산 패키지 (GOCI-III 편광 연구)

이 패키지는 C 의 `--water-model ocrt` 경로를 Python 으로 재현하고, GPU(CuPy) 또는 CPU(NumPy) 배치로
전체 격자 자료를 생산한다. 대기 SOS, 수중 SOS, 공기–물 결합, 성분 모델, 기체 흡수, 에어로졸 value 커널,
Rrs 가 모두 배치화되어 있다. 배치 결과는 단일 케이스 결과와 비트 수준(또는 유효숫자 7자리)으로 일치한다.

## 구성

- `ocrt_py/`: 복사전달 Python 모듈. 배치 엔진은 `atmos_batch.py`, 오케스트레이터는 `batch_driver.py`.
- `data/`: 입력 자료.
    - 에어로졸 `.mie` 표(이 패키지가 쓰는 96 모델 = OPAC 16종(C50/T50/M80C/M95C/M98C/O99 등) +
      AccuRT 유래 Ahmad 80종 `r<RH>f<FINE>v01`). 표는 GitHub Release `data-v1` 에서 설치한다
      (루트 `README.md` 참조).
    - 순수해수/식물플랑크톤/광물 IOP, AFGL 대기, 기체 흡수 단면적.
- `produce_grid.py`: **배치 자료생산 진입점.** 격자 CSV → 12 밴드 × R1/R2/R3.
- `gpu_gate.py`: GPU 무결성 게이트(NumPy 와 CuPy 를 유효숫자 7자리로 대조).
- `ocrt_solve.py`: 단일 케이스/격자 CLI(검증·소규모용).
- `verify.py`: C 코드 대비 검증 재확인.
- 실행 스크립트: `run_produce.bat/.sh`(CPU), `run_produce_gpu.bat/.sh`(GPU), `run_gpu_gate.bat`(게이트).
- `example_grid.csv`: 격자 CSV 형식 예시.

## 요구 환경

- CPU 실행: Python 3.8 이상과 NumPy (`pip install numpy`).
- GPU 실행: 위에 더해 CuPy. CUDA 버전에 맞게 설치한다. 예: `pip install cupy-cuda12x`(CUDA 12.x)
  또는 `cupy-cuda11x`.

## 격자 CSV 형식

열: `case_id, sza, vza, raa, wind, aerosol, aod865, chl, tsm, acdom440`
(각도는 도, aod865 는 865 nm 기준 AOD, chl 은 mg/m³, tsm 은 g/m³, acdom440 은 1/m).
`example_grid.csv` 참조. 밴드별 AOD 는 `.mie` 표의 분광 소광비로 자동 조정한다.

## 포함된 격자와 재생성

- `full_grid_design_v1.csv`: 1000 케이스 설계 격자(자료생산 입력).
- `iop_grid_design_v1.csv`: 해수 IOP 만 담은 부속 파일(검증·기록용).
- `make_grid.py` / `make_grid.bat` / `make_grid.sh`: 격자 재생성 스크립트.

설계: 해수 IOP(TSM 0.1–20, Chl 0.1–5, aCDOM440 0.01–0.1)는 로그 공간의 상관 다변량 정규분포에 범위
절단을 적용해 대양에서 연안까지 함께 변하도록 했다(목표 상관 TSM–Chl 0.85, Chl–aCDOM 0.75,
TSM–aCDOM 0.77). 에어로졸 6종 균등, AOD865 0.05–0.3 균등, 풍속 1–10 균등, 기하는 선글린트 없는
조건(Cox–Munk 직접 글린트 rho_g < 5e-4)이다.

```bash
# 격자 재생성 (시드 고정 = 재현 가능)
python make_grid.py --n 1000 --seed 20260718
# 또는 make_grid.bat / ./make_grid.sh
```

**주의.** 포함된 `full_grid_design_v1.csv` 는 원본 설계 스크립트와 시드가 소실되어 설계 사양(범위·상관·
글린트 기준)으로 다시 만든 것이다. 절단 때문에 실현 상관은 목표보다 조금 낮다(0.85→0.81, 0.75→0.70,
0.77→0.72). 개별 케이스는 원본과 다르지만 설계 사양은 같다.

## 산출 반사도 (행/케이스마다 세 번 실행, 차감으로 분리)

- R1 `rho_TOA`: 해양 + Rayleigh 대기(1013.25 hPa) + 기체 흡수 + 에어로졸(파장별 조정) + 성분 모델
  수출광.
- R2 `rho_R`: Cox–Munk Rayleigh 대기(에어로졸 없음), 선글린트 분리.
- R3 `rho_(R+A,black)`: 흑색 해양 위 Cox–Munk Rayleigh + 에어로졸, 선글린트 분리.
- 유도: `rho_RC = R1 − R2`, `rho_A + rho_RA = R3 − R2`, `t·rho_w = R1 − R3`, `Rrs = R1 의 Rrs(0+)`.
  항등식 `rho_RC = (rho_A + rho_RA) + t·rho_w` 가 성립한다.

출력 CSV 열: `case_id, band_nm, rho_TOA, rho_R, rho_RpA, rho_RC, rho_A_RA, t_rho_w, Rrs`.

## 사용법

### CPU (NumPy)

```bash
# Linux / macOS
./run_produce.sh full_grid_design_v1.csv result.csv

# Windows
run_produce.bat full_grid_design_v1.csv result.csv

# 직접 호출 (옵션 조정)
python produce_grid.py --grid full_grid_design_v1.csv --out result.csv \
    --data data --n-mu-water 24 --nt-atm 400 --m-max 16 --max-it-water 500
```

### GPU (CuPy)

```bash
# Linux / macOS
./run_produce_gpu.sh full_grid_design_v1.csv result.csv

# Windows
run_produce_gpu.bat full_grid_design_v1.csv result.csv
```

내부적으로 `OCRT_PY_GPU=1` 을 설정한다. 각 행이 독립이므로 배치 축(B)으로 쌓아 GPU 에서 병렬 처리한다.

### 재개와 밴드 선택

- 재개: 출력 CSV 에 이미 있는 `(case_id, band)` 는 자동으로 건너뛴다. 같은 `--out` 으로 다시 실행하면
  중단 지점부터 이어서 실행한다.
- 밴드 선택: `--bands 443,555,865` (기본값은 12 밴드 전체).
- 기체 흡수 끄기: `--no-gas` (기본값 켜짐).

## GPU 무결성 게이트 (로컬 필수)

CPU 와 GPU 는 부동소수 연산 순서가 다를 수 있다. 그래서 게이트는 결과가 비트 동일이 아니라 **유효숫자
7자리(상대오차 1e-7)** 로 일치하는지 확인한다. 이 게이트는 CuPy 가 설치된 로컬 CUDA 환경에서만 실행할 수
있다.

```bash
# 소규모 격자(4–16행)로 CPU 와 GPU 결과 대조
python gpu_gate.py --grid small_grid.csv --data data --bands 443,555,865
# 또는 Windows: run_gpu_gate.bat small_grid.csv
```

`rho_TOA/rho_R/rho_RpA/Rrs` 각각의 최대 상대오차와 PASS/FAIL, CPU/GPU 속도를 출력한다. 자료생산 전에 이
게이트를 통과시키기를 권한다.

## 옵션 기본값 (자료생산 조건)

- `--n-mu-water 24`: 수중 Gauss–Legendre 절점 수.
- `--m-max 16`: 대기 SOS 푸리에 모드 수.
- `--nt-atm 400`: 대기 SOS 층 수.
- `--max-it-water 500`: **수중 산란차수 상한. 혼탁(고산란) 밴드는 320차 이상이 필요하므로 옛 기본값 160
  으로는 수렴하지 않는다. 500 을 권한다(기본값).**
- `--chunk 64`: 한 배치의 행 수. GPU 메모리에 맞춰 조정한다.
- `--phyto` 를 생략하면 Stage-2 검증 기본 Chl 흡수 경로를 쓴다: 흡수는
  `data/water_iop/phyto_absorption_default.csv`, 식물플랑크톤 `b = bb = 0`, 입자 위상은 쇄설물/광물만.
- Chl > 0 에서 `--phyto` 를 명시하면 종 이름과 무관하게 fail-loud 한다. 17종 EAP 목록은 생성기와 향후
  검증용이며 성분 생산 RT 에 자동으로 들어가지 않는다.

## 검증 요약 (배치 vs 단일)

- R2(Rayleigh + 흑색 Fresnel 해면): 상대차 0 (비트 동일).
- R3(Rayleigh + 에어로졸, value 커널): 0~2e-16 (C50/M95C/Ahmad).
- R1(3패스 결합 ρ_TOA): rho_I 약 5e-16 (수렴 케이스), 유효숫자 7자리 이내.
- Rrs(0+): 약 2e-16 (Lu_0plus·Ed_0plus 모두 일치).
- 기체 흡수 ON: R3 5e-16, R1 0.
- 항등 층 패딩(층수가 다른 케이스를 한 배치에): 약 2e-16.

## 성능 주의

Python CPU 결합 계산은 자료생산 조건에서 케이스당 수십 초에서 수 분이 걸린다. 1000행 × 12밴드 =
12,000 케이스는 CPU 로는 오래 걸리므로 **GPU 실행을 권한다**. GPU 게이트를 먼저 통과시킨 뒤 전체 자료를
생산한다.

## 에어로졸 모델 주의

`data/` 의 Ahmad 80종(`r<RH>f<FINE>v01`)은 **AccuRT 유래 Ahmad 계열**이다(AccuRT 건조 굴절률·미세물리 +
Ahmad 2010 의 상대습도 상태·미세모드 분율·성장인자). Ahmad et al. (2010) 표 4 의 입자 반지름을 그대로
재현한 것이 아니다. 논문 원본 버전을 쓰려면 해당 `.mie` 파일을 추가해야 한다(이 패키지에는 AccuRT 유래
버전만 있다).

## 2026-07-21 전체격자 에어로졸 객체 수정

버전: `pyOCRT-v1.2-2026-08-09-chl-zero-extension-legacy-mie-cleanup`

대기–해양 각도 전체격자 경로는 다음으로 쓴다.

```bash
python ocrt_fullgrid.py --out grid.csv --wl 490 --sza 30 --wind 3 \
  --chl 1.2 --tsm 3.5 --acdom440 0.04 --aod865 0.2 --aer C50 \
  --vza-values 0,30 --raa-values 0,90,180,270 --no-gas
```

`.mie` 파일은 한 번 읽고, 기하 루프 전에 고정된 `AerosolRuntime` 을 한 번 준비한다. 같은 객체를 모든 셀과
두 대기 패스에 명시적으로 전달한다. 에어로졸 자료 없이 `AOD > 0` 을 주면 Rayleigh 단독/기본 상태를 조용히
쓰지 않고 오류를 낸다.

상세: `validation/python_fullgrid_aerosol_object_fix/VALIDATION_REPORT.md`.

## 2026-07-21 흑색 Fresnel 해면 이름 정리

Python 의 표준 상위 표면 모드는 `black_fresnel_ocean` 이다. 흑색 해양 위의 거친 Fresnel 계면이며 수출광은
없다. 옛 `coxmunk` 모드와 함수 이름은 폐기 예정 별칭으로 남긴다. 기본 경사분산은 정확한 Cox–Munk 모델이
아니라 OCRT 하한식 `sigma² = 0.003 + 0.00512 · max(0.01, W)` 로 별도 문서화한다.

## EAP 위상 자료 갱신 (2026-07-26)

옛 `pico_Synechococcus_EAP.mie`, `nano_Haptophytes_EAP.mie`, `Diatoms_centric_EAP.mie` 를 제거했다.
재생성한 대표 EAP 파일 17개는 `data/water_iop/eap/` 에 있다. 옛 이름 `pico`, `nano`, `micro` 는 새
Synechococcus, Prymnesiaceae, 중심 규조 파일의 호환 별칭으로 남긴다. 원자료와 재현 가능한 C 생성기는
`data/eap_source/` 와 `tools/eap_phase_generator/` 에 있다.

## EAP 식물플랑크톤 Stage-2 계약

Chl > 0 에서 종 옵션을 생략하면 검증된 기본 흡수 스펙트럼을 쓰고 식물플랑크톤 산란은 0 으로 둔다.
`--phyto`/`--phyto-group` 을 명시하면 `OCRT_ADVANCED=1` 여부와 무관하게 fail-loud 한다. 17종 EAP 위상
파일과 생성기는 패키지에 남기지만 성분 생산 RT 에는 연결하지 않는다.

## 선택적 직접 수중 value 커널 (2026-07-27)

이 경로들은 선택 사항이며 기본 계수 커널을 바꾸지 않는다.

- `OCRT_VALUE_PHASE_SPLINE=1`: `mu = cos(theta)` 에 대한 3차 Hermite 방식 스플라인 보간과 Gauss 정규화를
  쓰는 스칼라 직접 value 커널.
- `OCRT_WATER_VALUE_KERNEL_POL=1`: 호출자가 원시 `value_phase` 자료를 주면 쓰는 직접 편광 value 커널
  (`pfm, gr, gt, arr, art, att`).
- `OCRT_WATER_MIE_MOMENT_N_MU=400`: 정규화 구적 수.
- `OCRT_WATER_VALUE_NPHI=N`: 명시적 방위각 적분 수.

정확히 VZA = 0 은 편광 value 커널 출력 추출의 알려진 특이점이다. VZA = 0.001° 를 쓴다(검증된 안전 범위
0.001–0.01°). 천저 우회로 약 1e-6° 미만의 각도를 쓰지 않는다. 기하를 자동으로 바꾸지 않는다.

OCRT 성분 Chl 게이트는 1100 nm 까지 확장했다. 식물플랑크톤 흡수표는 350–1100 nm 를 명시적으로 덮는다.
원자료 기반 값은 850 nm 에서 끝나고, 승인된 850–1100 nm 확장은 정확히 0 이다. 순수해수 a_w 와 총 b_w 는
200–2449 nm 로 표화되어 있다. 별도의 온도 보정표는 아직 400–900 nm 로 제한되므로, 그 구간 밖에서 20 °C 가
아닌 계산은 향후 psi_T 확장이 필요하다.
