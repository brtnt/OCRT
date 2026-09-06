# OCRT polarization sensitivity run guide — three-part split

- Date: 2026-08-02
- Version: `OCRT-v1.2-2026-08-02-KST-dtpsign-phase-diagnostic`
- Target: Windows, 24 cores
- Files

| file | purpose |
|---|---|
| `ocrt_v1.2.exe` | the executable; place it in `<package root>\build\` |
| `run_star.py` | the runner |
| `run_part1.bat` `run_part2.bat` `run_part3.bat` | one batch file per part |
| `merge_parts.py` `merge_parts.bat` | merge the three parts |
| `ocrt_win_compat.c` `build_win.bat` | needed only for a rebuild |
| this document | |

---

## 1. What is computed

A star-shaped design of **615** coupled-grid runs. One run gives 18 viewing zenith angles × 72
azimuths, i.e. **1,296 directions**, at once.

| branch | scanned | fixed | runs |
|---|---|---|---|
| N | 27 water states (Chl 3 × TSM 3 × CDOM 3) | central atmosphere (r50f05v01, AOD 0.10) | 405 |
| S | 15 atmospheres (3 models × 5 AOD levels) | central water (Chl 1.0, TSM 1.5, CDOM 0.035) | 225 |
| shared | central atmosphere × central water | — | −15 |
| **total** | | | **615** |

Common conditions: SZA 25·50·75°, bands 412·555·667·748·865 nm, wind 5 m/s, sun glint decoupled,
pseudo-spherical approximation on.

Levels of the three water variables (full cross, so all corners that break the TSM–Chl
correlation are included):

- Chl: 0.15, 1.0, 4.0 mg/m³
- TSM: 0.15, 1.5, 15 g/m³
- CDOM (440 nm): 0.012, 0.035, 0.09 1/m

## 2. Preparation

**The input data are already available.** The uploaded `OCRT_v1_2_DTPSIGN_PHASE_DIAGNOSTIC_2026-08-02.zip`
contains the source and the `inputs\` data. Its unpacked folder is the package root below (the
folder that holds both `inputs` and `src`).

1. Python 3.8 or later with a working `python` command. No extra library is needed.
2. Copy `ocrt_v1.2.exe` into `<package root>\build\`.
3. Put `run_star.py` and `run_part*.bat` in the same `build\` folder. The script walks upward to
   find the package root, so no path needs to be written.
4. The default output folder is `<package root>\star_run`. Use `--out` to change it.

**Delete the results of the previous version.** The 2026-08-02 version fixed the sign of the
odd azimuthal modes, which changes the U component of the coupled runs. Old and new results must
not be mixed.

## 3. Run order

### 3.0 Check the executable (once)

```
.\run_part1.bat --verify
```

Two reference cases are run and compared with stored values. Both must report "pass".

### 3.1 Check the plan

```
.\run_part1.bat --dry-run
```

Shows the number of runs assigned and the output folder.

### 3.2 Measure the real duration (recommended)

```
.\run_part1.bat --limit 24
```

24 runs are one batch on 24 cores. The "average per run" printed at the end is the real value on
that machine. These 24 runs are skipped automatically in the main run, so nothing is wasted.

### 3.3 Main run

Run the three parts on three machines, or one after another on one machine.

```
.\run_part1.bat
.\run_part2.bat
.\run_part3.bat
```

The parts are exactly 205 runs each and do not overlap. Bands and SZA are mixed evenly across the
parts, so one finished part already shows the overall structure.

### 3.4 Merge

```
.\merge_parts.bat --copy
```

Edit the three folder paths at the top of the batch file to the folders where each part ran, then
execute. If all parts ran on one machine, the three paths may be identical. The merge

- joins the per-part records into one `index.csv`;
- finds and reports missing runs, failures, non-converged runs and duplicates;
- warns if the executable or the angle grid differs between parts;
- with `--copy`, also collects the result CSV files into one folder.

## 4. Reading the screen

One line is refreshed after every run.

```
[####......] 82/205  40.0% | left 00:21:33 | done 14:52 | avg 1.5 min | parallel 24
```

The remaining time is computed from the measured average per run and the number of parallel
runs. It starts from an assumed 9 minutes and updates as measurements accumulate. Failures or
non-converged runs are appended only when they occur.

## 5. Stop and resume

- `Ctrl+C` stops new runs, waits for the running ones to finish, cleans up and exits.
- Running the same batch file again skips finished runs and runs only the rest.
- Only completed files are moved to the final location. Files being computed stay in
  `out\ocn\_tmp\` and are moved only on success, so no half-written file remains as a result.

## 6. Output layout

```
<output folder>\
  index_part1of3.csv            per-part record; one row = one run
  runlist_*.csv                 planned run lists
  run_*.log                     detailed log of each run
  bin\                          copy of the executable and its sha256
  out\ocn\<name>.csv             results (one file = all 1,296 directions of one condition)
  out\ocn\_tmp\                 temporary files during computation (empty when done)
```

The file name carries every axis. For example

```
O_b555_s50_chl1_tsm1p5_cdm0p035_mr50f05v01_a0p1.csv
```

is 555 nm, SZA 50°, Chl 1.0, TSM 1.5, CDOM 0.035, model r50f05v01, AOD 0.10. A decimal point is
written as `p`. Because the AOD is in the name, adding other values later does not overwrite
anything.

Record columns:

```
run_id, branch, band_nm, sza_deg, chl, tsm, cdom, aer_model, ssa443, aod865,
wind_ms, pssa, glint, vza_step_deg, raa_step_deg, out_path, status, wall_s,
n_rows, n_unconverged, max_water_orders, code_version, binary_sha256,
run_tag, finished_at
```

`branch` is N, S, or NS for the central run shared by both.

## 7. Failures and non-convergence

- **Failure** (`status=failed`): no result file remains. Check the reason in `run_*.log` and run
  the same batch file again; only the failed runs are repeated.
- **Non-convergence** (`status=unconverged`): a file remains, but some rows hit the in-water
  scattering-order limit. Delete that result file and rerun with `--max-orders 800`.

The most turbid combination (Chl 4, TSM 15, CDOM 0.09) at 412 nm and SZA 75° converged in all
1,296 rows with a maximum order of 87, so the default should be enough.

## 8. Expected duration

On the measuring machine (1 core) one run takes about 88 s. A denser angle grid adds almost no
time, so all 1,296 directions are produced.

On 24 cores the 615 runs should take about 40 minutes. This is based on the measuring machine, so
measure the real value first with step 3.2. One part takes one third.

Output size is about 1.3 MB per run, about 0.8 GB for 615 runs.

## 9. Changing the design

Only the section "1. experiment design constants" at the top of `run_star.py` needs editing.
Bands, SZA, the levels of the three water variables, aerosol models and AOD, wind and the angle
grid spacing are all there.

The angle grid can also be changed with `--vza-step` and `--raa-step`; the values used are
recorded. To change the number of parts use e.g. `--part 1/5`. To add AOD values later use
`--n-aod 0.05 --n-aod 0.3`; existing files are untouched and only new files are added.

## 10. Note on the Windows build

This version does not build on Windows as is: `setenv` and `unsetenv` are missing at link time
(see the separate bug report). The supplied `ocrt_v1.2.exe` was built with `ocrt_win_compat.c`,
which provides those two functions.

To rebuild, put `ocrt_win_compat.c` and `build_win.bat` in the package root and run
`build_win.bat`. MSYS2 MinGW64 gcc is required. The architecture is set to `x86-64-v3` (AVX2). The
Core Ultra 9 275HX has no AVX-512, so do not build with `cascadelake`.

---

# OCRT 편광 민감도 실행 안내 — 3조각 분할판

- 작성일: 2026-08-02
- 판본: `OCRT-v1.2-2026-08-02-KST-dtpsign-phase-diagnostic`
- 대상: Windows, 24코어
- 파일 목록

| 파일 | 쓰임 |
|---|---|
| `ocrt_v1.2.exe` | 계산 실행 파일. `<패키지루트>\build\` 에 둔다 |
| `run_star.py` | 실행기 본체 |
| `run_part1.bat` `run_part2.bat` `run_part3.bat` | 조각별 실행 |
| `merge_parts.py` `merge_parts.bat` | 세 조각 합치기 |
| `ocrt_win_compat.c` `build_win.bat` | 직접 다시 빌드할 때만 필요 |
| 이 문서 | |

---

## 1. 무엇을 돌리는가

별 모양 설계로 결합 격자 계산 **615건**을 돌린다. 한 건이 관측천정각 18개 × 방위각 72개, 곧 **1,296방향**을
한 번에 낸다.

| 가지 | 무엇을 훑는가 | 무엇을 고정하는가 | 건수 |
|---|---|---|---|
| N | 해수 27종 (엽록소 3 × 총부유물 3 × CDOM 3) | 중심 대기 (r50f05v01, 광학두께 0.10) | 405 |
| S | 대기 15종 (3모델 × 광학두께 5단) | 중심 해수 (엽록소 1.0, 총부유물 1.5, CDOM 0.035) | 225 |
| 공유 | 중심 대기 × 중심 해수 | — | −15 |
| **합계** | | | **615** |

공통 조건은 태양천정각 25·50·75도, 밴드 412·555·667·748·865 nm, 풍속 5 m/s, 선글린트 분리, 의사구면근사
켜짐이다.

해수 세 값의 수준은 다음과 같다. 완전 교차이므로 총부유물·엽록소의 상관을 깨는 모서리가 모두 들어간다.

- 엽록소: 0.15, 1.0, 4.0 mg/m³
- 총부유물: 0.15, 1.5, 15 g/m³
- CDOM(440 nm): 0.012, 0.035, 0.09 1/m

## 2. 준비

**계산 자료는 이미 가지고 있다.** 업로드한 `OCRT_v1_2_DTPSIGN_PHASE_DIAGNOSTIC_2026-08-02.zip` 안에 소스와
`inputs\` 자료가 모두 들어 있다. 그 압축을 푼 폴더가 아래에서 말하는 패키지 루트다(`inputs` 와 `src` 가
함께 있는 곳).

1. Python 3.8 이상이 있고 `python` 명령이 동작해야 한다. 별도 라이브러리는 필요 없다.
2. `ocrt_v1.2.exe` 를 `<패키지루트>\build\` 에 복사한다.
3. `run_star.py` 와 `run_part*.bat` 을 같은 `build\` 폴더에 둔다. 스크립트가 위로 올라가며 패키지 루트를
   스스로 찾으므로 경로를 적을 필요가 없다.
4. 산출물 폴더는 기본값이 `<패키지루트>\star_run` 이다. 바꾸려면 `--out` 으로 지정한다.

**이전 판본 결과는 반드시 지운다.** 2026-08-02 판본에서 홀수 방위각 모드의 부호가 고쳐져 결합 계산의 U
성분이 달라졌다. 옛 결과와 섞으면 안 된다.

## 3. 실행 순서

### 3.0 실행 파일 확인 (한 번만)

```
.\run_part1.bat --verify
```

기준 사례 두 건을 돌려 미리 재둔 값과 대조한다. 두 건 모두 "통과" 가 나와야 한다.

### 3.1 계획 확인

```
.\run_part1.bat --dry-run
```

맡은 건수와 산출물 폴더 위치가 나온다.

### 3.2 실제 소요 재기 (권장)

```
.\run_part1.bat --limit 24
```

24건이면 24코어에서 한 묶음이 돈다. 끝날 때 나오는 "1건 평균" 이 그 기계의 실제 값이다. 여기서 돌린 24건은
본 실행에서 자동으로 건너뛰므로 낭비가 없다.

### 3.3 본 실행

세 조각을 세 기계에 나눠 돌리거나, 한 기계에서 차례로 돌린다.

```
.\run_part1.bat
.\run_part2.bat
.\run_part3.bat
```

조각은 205건씩 정확히 삼등분되고 서로 겹치지 않는다. 밴드와 태양천정각이 조각마다 고르게 섞이므로, 한
조각만 먼저 끝나도 그것만으로 전체 구조를 미리 살펴볼 수 있다.

### 3.4 합치기

```
.\merge_parts.bat --copy
```

배치 파일 위쪽의 세 폴더 경로를 각 조각이 돌아간 곳으로 고친 뒤 실행한다. 한 기계에서 차례로 돌렸다면 세
경로가 모두 같아도 된다. 합치기가 하는 일은 다음과 같다.

- 조각별 기록을 하나의 `index.csv` 로 합친다.
- 빠진 실행, 실패, 미수렴, 중복을 찾아 알린다.
- 조각마다 실행 파일이나 각도 격자가 다르면 경고한다.
- `--copy` 를 주면 결과 CSV 파일도 한 폴더로 모은다.

## 4. 화면 보는 법

한 건이 끝날 때마다 아래 한 줄이 갱신된다.

```
[####......] 82/205  40.0% | 남은 00:21:33 | 완료 14:52 | 평균 1.5분 | 동시 24
```

남은 시간은 지금까지 실측된 1건 평균과 동시 실행 수로 계산한다. 처음에는 9분으로 가정했다가 실측이 쌓이면
저절로 갱신된다. 실패나 미수렴이 생길 때만 뒤에 항목이 붙는다.

## 5. 중단과 재개

- `Ctrl+C` 를 누르면 새 실행을 멈추고, 돌고 있는 건이 끝날 때까지 기다린 뒤 정리하고 종료한다.
- 다시 같은 배치 파일을 실행하면 끝난 건은 건너뛰고 남은 건만 돌린다.
- 완성된 파일만 최종 위치로 옮긴다. 계산 중인 파일은 `out\ocn\_tmp\` 에 있다가 성공했을 때만 옮겨지므로,
  중간에 끊긴 반쪽 파일이 결과로 남지 않는다.

## 6. 산출물 구조

```
<산출물폴더>\
  index_part1of3.csv            조각별 기록. 한 행 = 한 실행
  runlist_*.csv                 계획된 실행 목록
  run_*.log                     한 건씩의 상세 기록
  bin\                          실행 파일 사본과 sha256
  out\ocn\<이름>.csv             결과 (한 파일 = 한 조건의 전 각도 1,296방향)
  out\ocn\_tmp\                 계산 중 임시 파일 (끝나면 비어 있다)
```

결과 파일 이름에 모든 축이 들어간다. 예를 들어

```
O_b555_s50_chl1_tsm1p5_cdm0p035_mr50f05v01_a0p1.csv
```

는 555 nm, 태양천정각 50도, 엽록소 1.0, 총부유물 1.5, CDOM 0.035, 모델 r50f05v01, 광학두께 0.10 이다.
소수점은 `p` 로 적는다. 광학두께가 이름에 들어 있으므로 나중에 다른 값을 더해도 덮어쓰지 않는다.

기록의 열은 다음과 같다.

```
run_id, branch, band_nm, sza_deg, chl, tsm, cdom, aer_model, ssa443, aod865,
wind_ms, pssa, glint, vza_step_deg, raa_step_deg, out_path, status, wall_s,
n_rows, n_unconverged, max_water_orders, code_version, binary_sha256,
run_tag, finished_at
```

`branch` 는 N, S, 또는 두 가지가 공유하는 중심 실행이면 NS 다.

## 7. 실패와 미수렴 처리

- **실패**(`status=failed`): 결과 파일이 남지 않는다. `run_*.log` 에서 사유를 확인한 뒤 같은 배치 파일을 다시
  돌리면 실패한 건만 재실행된다.
- **미수렴**(`status=unconverged`): 파일은 남지만 수중 계산이 산란차수 상한에 걸린 행이 있다는 뜻이다. 그
  결과 파일을 지우고 `--max-orders 800` 을 붙여 다시 돌린다.

가장 혼탁한 조합(엽록소 4, 총부유물 15, CDOM 0.09)을 412 nm, 태양천정각 75도에서 확인했을 때 1,296행 모두
수렴했고 최대 산란차수가 87이었다. 기본값으로 충분할 것으로 본다.

## 8. 예상 소요

측정 기계(1코어)에서 1건이 약 88초다. 각도 격자를 촘촘히 해도 시간이 거의 늘지 않는 구조이므로 1,296방향을
그대로 낸다.

24코어에서 615건이면 대략 40분 안팎으로 본다. 다만 이는 측정 기계 기준이므로 3.2단계로 실제 값을 먼저 재기
바란다. 조각 하나는 그 3분의 1이다.

산출물 용량은 1건에 약 1.3 MB, 615건이면 약 0.8 GB 다.

## 9. 설계를 바꾸려면

`run_star.py` 위쪽의 "1. 실험 설계 상수" 구역만 고치면 된다. 밴드, 태양천정각, 해수 세 값의 수준, 에어로졸
모델과 광학두께, 풍속, 각도 격자 간격이 모두 그 안에 모여 있다.

각도 격자는 `--vza-step` 과 `--raa-step` 으로도 바꿀 수 있고, 쓴 값이 기록에 남는다. 조각 수를 바꾸려면
`--part 1/5` 처럼 주면 된다. 나중에 광학두께를 더하려면 `--n-aod 0.05 --n-aod 0.3` 을 붙인다. 이때 기존
파일은 건드리지 않고 새 파일만 늘어난다.

## 10. Windows 빌드에 대한 단서

이 판본도 그대로는 Windows 에서 빌드되지 않는다. `setenv` 와 `unsetenv` 가 링크 단계에서 빠진다(별도 버그
보고서 참조). 함께 보낸 `ocrt_v1.2.exe` 는 `ocrt_win_compat.c` 로 그 두 함수를 채워 만든 것이다.

직접 다시 빌드하려면 `ocrt_win_compat.c` 와 `build_win.bat` 을 패키지 루트에 두고 `build_win.bat` 을
실행한다. MSYS2 의 MinGW64 gcc 가 필요하다. 아키텍처는 `x86-64-v3`(AVX2) 로 잡혀 있다. Core Ultra 9 275HX
는 AVX-512 를 지원하지 않으므로 `cascadelake` 로 빌드하면 안 된다.
