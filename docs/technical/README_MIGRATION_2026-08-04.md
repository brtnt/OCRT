# OCRT paper-work migration package

- Created: 2026-08-04
- Target version: `OCRT-v1.2-2026-08-02-KST-dtpsign-phase-diagnostic`
- Windows executable sha256: `6008d0a3e5b1afcc371ecfeea4f2ac11351d07e1eb49e0ec3fe6cff1d65e264e`

This one package lets the work on two papers continue: the **OCRT model description paper** and
the **GOCI-III polarization atmospheric-correction paper**.

---

## 0. What must be present as well

**The OCRT C source and input data are not in this package** (137 MB, already available). The
unpacked folder of the following archive is needed:

```
OCRT_v1_2_DTPSIGN_PHASE_DIAGNOSTIC_2026-08-02.zip
```

Its unpacked folder is called the "package root" below (the folder that holds both `inputs` and
`src`).

**The 8,045 raw runs of the polarization campaign (about 10 GB) are not included either.** They
are in the local campaign folder (`star_run`); this package contains a 55 MB reduced version,
which is enough to reproduce every analysis done so far.

---

## 1. Folder layout

```
01_paper/              OCRT model description paper
02_paper_figures/      its figures and source data
03_polarization/       data, scripts and figures of the polarization paper
04_run_tools/          simulation run tools
05_build/              Windows executable and build helpers
06_bugreports/         bug reports for the OCRT development session
07_docs/               earlier migration documents and raw-data bundles
```

---

## 2. OCRT model description paper (01_paper, 02_paper_figures)

### Status

- Manuscript: `01_paper/OCRT_reference_paper_draft.docx` (23 pages, 5 figures)
- Build script: `01_paper/build_ocrt_paper.js`
- Figure sources: `01_paper/paper_figs/` (5 files; the build reads them from here)

### Rebuild

```
NODE_PATH=<npm global path> node build_ocrt_paper.js
```

docx-js v9.6.1 is required. The output is `OCRT_reference_paper_draft.docx`. A changed figure must
be copied into `paper_figs/` to appear in the build.

### Done in this session

The pure-water validation figure (old Figure 4) was removed and the later figures were
renumbered. Two reasons: first, its data were produced before the 2026-07-23 correction of the
underwater output azimuth convention and no longer match the current code; second, the CDOM
series of the constituent comparison already contains the pure-water limit (Chl and TSM both 0).
One sentence in the text states this.

### Remaining

Three `[TO COMPLETE]` marks remain in the text:

1. the GOCI-III application example figure in Section 7 — **the results of 03_polarization go
   here**
2. the archive DOI in the data-availability section
3. the author list and acknowledgements

Three places where the text and the code disagree were not yet touched:

- the aerosol library is described as two sets (16 and 80 models), but the package holds a third
  set (`aerosol_ahmad2010_paper_mie`, A2010ver, 80 models)
- the surface option is written as `--surface coxmunk`, but the canonical name is
  `black_fresnel_ocean` and `coxmunk` is a deprecated alias
- the reproducible-build section does not yet say that `-march=native` was fixed to
  `-march=cascadelake`

### Open background items

- restoring phytoplankton scattering (Chl currently contributes absorption only)
- forward-peak truncation in the constituent path is not implemented (4–6 % bias in large grid
  production)

---

## 3. Polarization paper (03_polarization)

### What was found

The single-scattering albedo (absorption) and the particle size both move the spectral slope of
the intensity reflectance, so intensity alone cannot separate them. This is the degeneracy. The
study asks whether polarization resolves it.

**Main results**

1. The Ahmad 80-model library contains this degeneracy. With the fine-mode fraction fixed and only
   the relative humidity changed, the Ångström exponent stays almost the same while the albedo
   changes a lot. The 50 % fine-mode row is the purest degeneracy: Ångström 1.81–1.91 with albedo
   0.947–0.984.

2. In atmospheric correction of low-turbidity water with the two near-infrared bands (748, 865)
   there are two observations and three unknowns (size, absorption, amount), so the problem is
   under-determined. The required accuracy was derived backwards without assuming a noise level.

   | SZA | required in intensity reflectance | required in DoLP | relaxation factor |
   |---|---|---|---|
   | 25° | 1.5e-05 | 5.6e-03 | 381 |
   | 50° | 1.8e-05 | 4.9e-03 | 273 |
   | 75° | 8.8e-05 | 5.8e-03 | 66 |

   A realistic intensity accuracy is 2 % of the correction, about 5e-4, so intensity alone never
   reaches the requirement at any azimuth. Polarization needs 5e-3, which current technology can
   reach.

3. The azimuth dependence is clear. **The anti-solar side (135–180°) is the most forgiving**; at
   SZA 50° up to 1.1e-02 is allowed. At SZA 25° the requirement stays at 4–6e-03 for any azimuth,
   because at small SZA the azimuth loses its leverage in the scattering-angle formula.

4. The black-surface idealization is justified for low-turbidity water: adding clear water as an
   unknown changes the retrieval uncertainty by a factor of only 1.00.

### The method was replaced twice; the record is kept

**First method (withdrawn)** — the observations were approximated as linear in the unknowns and
the error was propagated by optimal estimation. Two problems:

- the AOD range spans a factor of 30, so the residual of the linear fit was 7.5 times the
  measurement noise; the values could not be trusted;
- as a result, the azimuth maps showed high-frequency speckle that is not physics. The roughness
  (second difference ÷ first difference) was 0.19 in the raw data but 0.56 in the metric.

**Second method (current)** — no approximation. Each of two models chooses its AOD freely so that
the 748 nm intensities agree; the remaining 865 nm intensity difference and the DoLP difference are
measured at that point. The required accuracy is that value divided by 3. No noise is assumed.

The diagnosis of the azimuth speckle is kept (`figures/diag_방위각_매끄러움.png`): the raw data are
smooth (modes of order 9 and higher: 1.6 %), the difference between the two models has real
structure with periods of 45–120°, and the speckle arose in the least-squares fitting step.

### Data (03_polarization/data)

| file | contents |
|---|---|
| `atm_reference.csv.gz` | all-angle reflectance of the atmosphere without ocean (black Fresnel): 24 aerosol models × 7 AOD levels × 5 bands × 3 SZA × 1,296 directions |
| `cases_thin.csv.gz` | raw values of the coupled-ocean runs on a coarse angle grid (viewing zenith 15°, azimuth 45°, 48 directions) |
| `metrics.csv.gz` | derived metrics (N, S, usefulness). **Made with the first method; do not use the quantitative values** |
| `runs_summary.csv` | status, duration, convergence, version and hash of the 8,045 runs |
| `smooth.csv` | final product of the current method: 3,888 directions × (required intensity accuracy, required polarization accuracy) |
| `aerlib.csv` | Ångström exponent, albedo and asymmetry factor of the 80 aerosol models |
| `README_reduced.md` | column guide of the reduced results |

### Pure-Rayleigh reference (03_polarization/rayleigh)

Reflectance of an atmosphere without aerosol: 748·865 nm × SZA 25/50/75° = 6 files, each with
1,296 directions. **They are required to build the Rayleigh-corrected reflectance and do not exist
locally.** They were made in this session; if lost, they must be recomputed.

The command used:

```
ocrt_v1.2 --sza <25|50|75> --wavelength <748|865> \
  --surface black_fresnel_ocean --wind-speed 5 --decouple-sunglint --pssa \
  --lut-vza-step 5 --lut-raa-step 5 --output-full-grid R_<band>_<sza>.csv
```

No aerosol option may be given at all; `--aod-865 0` causes an error.

### Scripts (03_polarization/scripts)

| file | purpose |
|---|---|
| `smooth_signal.py` | computes the current metric and writes `smooth.csv` |
| `make_final_figs.py` | figures F1–F3 |
| `make_map.py` | direction maps F4·F5 (filled format) |

The figures of the two scripts do not overlap; only `make_map.py` makes the direction maps.

They run in the folder that holds the reduced data and use the paths `/tmp/ray/`,
`/tmp/aerlib.csv` and `/tmp/smooth.csv`, so **the paths must be adapted to the local environment.**

### Figures (03_polarization/figures)

| file | contents |
|---|---|
| `F1_축퇴격자와_대표짝.png` | degeneracy structure and polarization of representative pairs |
| `F2_필요정확도.png` | required accuracy in intensity and polarization (key figure) |
| `F3_방위각의존.png` | required accuracy by azimuth |
| `F4_방향지도.png` | viewing-direction map (filled format) |
| `F5_방향지도_두기준.png` | half-of-pairs criterion and hardest-pair criterion |
| `diag_방위각_매끄러움.png` | basis of the method choice; not for the paper |

---

## 4. Run tools (04_run_tools)

Placed in the `build` folder of the package root, they find their location by themselves.

| file | purpose |
|---|---|
| `run_star.py` | polarization campaign runner; stages (base, p4a–p4c, p5a–p5, p6, p6a) are selectable |
| `run_part1~3.bat` | the first run in three parts |
| `run_part4.bat`, `run_part4-2.bat` | stage-4 extension |
| `run_part6.bat` | degeneracy grid and reference atmosphere |
| `reduce_grid.py` | reduces the 10 GB raw data to 55 MB |
| `merge_parts.py/.bat` | merges part results and finds gaps |
| `make_ac_set.py` | **atmospheric-correction validation data generator (next task)** |

The stage-5 batch files (wind and SZA extension) were withdrawn because the design changed when
the degeneracy grid was given priority. They can be revived with `run_star.py --stage p5a` etc.

---

## 5. Next steps

### (a) Build the atmospheric-correction validation data set — ready

Run with `04_run_tools/make_ac_set.py`.

```
python make_ac_set.py --dry-run     check the plan
python make_ac_set.py --n 20        small trial (about 2 min)
python make_ac_set.py               main run (24 cores, about 7 h)
```

1,000 cases × 12 GOCI-II bands (380–865 nm), five runs per combination. The result table holds the
geometry angles, Rho_TOA, Rho_R, Rho_C, Rho_Am, Td_R_s, Td_R_v, Td_Am_s, Td_Am_v and nadir Rrs.

Two automatic checks are attached: `Td_check` within 3e-7 of 1, and `Rho_C = Rho_Am + t·ρ_w` with a
relative difference of 1.6e-15.

**Two lines for the data description**

1. At 865 nm the chlorophyll absorption is treated as 0, because the table ends at 850 nm; pure
   water absorbs 4.6 /m there and dominates, so there is no effect. Detritus scattering is fully
   active (b equals the other wavelengths).
2. The detritus backscatter ratio keeps its 850 nm value; its spectral change over that interval is
   0.07 %, so the error is below 0.01 %.

### (b) Finish the polarization paper

Place F1–F5 in the text and describe the results above. Two remaining analysis candidates:

- band-combination optimization: which of the five bands separate most efficiently gives the basis
  for polarization band selection;
- wind and SZA extension, revivable with `run_star.py --stage p5a`.

### (c) Finish the model description paper

Resolve the three `[TO COMPLETE]` marks and the three inconsistencies of Section 2.

### (d) Bug-report follow-up

Two reports are in `06_bugreports/`.

1. **Chlorophyll wavelength limit** — only the limit constant was raised to 1100 nm; the table did
   not follow, so chlorophyll absorption silently becomes 0 above 850 nm. Judged harmless for the
   current use, but recorded.
2. **Windows build** — `setenv` / `unsetenv` are missing at link time. Worked around with
   `05_build/ocrt_win_compat.c`; not yet fixed.

---

## 6. Things to know

### Environment

- The CPU is a Core Ultra 9 275HX **without AVX-512.** A `-march=cascadelake` build does not run
  at all. `05_build/build_win.bat` uses `x86-64-v3` (AVX2) by default.
- Batch files must be pure ASCII with CRLF line endings. Korean characters make cmd miscount byte
  positions and execute the wrong line.
- In PowerShell a batch file needs the `.\` prefix.
- Python is used from the Anaconda prompt.

### Mismatched data fail silently

Replacing only the executable while keeping the old `inputs` has actually happened. `run_star.py`
therefore compares the hashes of the main data files; `--verify` and the start of the main run
catch it.

### Scattering-angle and sun-glint conventions

- Scattering angle: cos Θ = −cos(SZA)·cos(VZA) − sin(SZA)·sin(VZA)·cos(RAA). Confirmed by the
  Rayleigh limit, where the DoLP maximum falls exactly at 90°.
- Sun-glint tilt angle: cos β = (cos SZA + cos VZA) / √(2 + 2·s·v). Zero means the direction of
  direct glint on a flat surface. At SZA 40° it is zero at VZA 40°, RAA 180°, which coincides with
  the intensity maximum of a run with glint enabled.

### Upward transmittance

It is defined only when water-leaving radiance exists, so a black-Fresnel run gives `nan`. The
Rayleigh-only upward transmittance needs a separate **coupled run without aerosol**; this is the
fourth run of `make_ac_set.py`.

---

## 7. Migration check

Unpacked in an isolated folder and reproduced:

- paper rebuild: the body (`word/document.xml`) of the docx made by `build_ocrt_paper.js` has the
  same hash as the enclosed manuscript;
- polarization figures: F1–F5 and the diagnostic figure, six files, byte-identical.

The polarization scripts use absolute paths (`/tmp/ray/` etc.); after moving them, change the
paths to relative ones and run from the `data` folder:

```
sed -i "s|/tmp/aerlib.csv|../data/aerlib.csv|g; s|/tmp/ray/|../rayleigh/|g; \
        s|/tmp/smooth.csv|../data/smooth.csv|g; \
        s|/mnt/user-data/outputs/|../figures/|g" \
    ../scripts/*.py
```

The paper build script also needs `FIG_DIR` and the output path adapted.

---

## 8. Reproduction check

After the migration, the following shows whether the environment is right.

```
cd <package root>\build
python run_star.py --verify
```

The data-file check and two reference cases are run. Reference values:

| case | I | Q | U |
|---|---|---|---|
| atmosphere only (black Fresnel, 555 nm, SZA 25°, VZA 30°, RAA 90°, AOD 0.1, r50f05v01) | 4.4991861596e-02 | −9.0983939253e-04 | −7.5031808484e-03 |
| the same with coupled ocean (Chl 1.0, TSM 1.5, CDOM 0.035) | 1.2305024573e-01 | −6.3070476021e-04 | −1.1478355407e-02 |

Both must pass with a relative difference of 0.

---

# OCRT 논문 작업 이관 패키지

- 만든 날: 2026-08-04
- 대상 판본: `OCRT-v1.2-2026-08-02-KST-dtpsign-phase-diagnostic`
- Windows 실행 파일 sha256: `6008d0a3e5b1afcc371ecfeea4f2ac11351d07e1eb49e0ec3fe6cff1d65e264e`

이 패키지 하나로 두 논문의 작업을 이어받을 수 있다. 하나는 **OCRT 모델 기술 논문**이고, 다른 하나는
**GOCI-III 편광 대기보정 논문**이다.

---

## 0. 함께 있어야 하는 것

이 패키지에 **OCRT C 소스와 입력 자료는 들어 있지 않다.** 137 MB 로 크고 이미 가지고 있기 때문이다. 아래
압축 파일을 푼 폴더가 필요하다.

```
OCRT_v1_2_DTPSIGN_PHASE_DIAGNOSTIC_2026-08-02.zip
```

푼 폴더를 아래에서 "패키지 루트"라 부른다. `inputs` 와 `src` 가 함께 있는 곳이다.

편광 캠페인의 **원자료 8,045건(약 10 GB)도 들어 있지 않다.** 로컬 캠페인 폴더(`star_run`)에 있으며, 이
패키지에는 그것을 55 MB 로 줄인 축약본이 들어 있다. 축약본만으로 지금까지의 분석을 모두 재현할 수 있다.

---

## 1. 폴더 구성

```
01_paper/              OCRT 모델 기술 논문
02_paper_figures/      그 논문의 그림과 원자료
03_polarization/       편광 논문의 자료·스크립트·그림
04_run_tools/          시뮬레이션 실행 도구
05_build/              Windows 실행 파일과 빌드 보조
06_bugreports/         OCRT 개발 세션에 전달할 버그 보고서
07_docs/               이전 이관 문서와 원자료 묶음
```

---

## 2. OCRT 모델 기술 논문 (01_paper, 02_paper_figures)

### 상태

- 원고: `01_paper/OCRT_reference_paper_draft.docx` (23쪽, 그림 5개)
- 빌드 스크립트: `01_paper/build_ocrt_paper.js`
- 그림 원본: `01_paper/paper_figs/` (5개, 빌드가 여기서 읽는다)

### 다시 빌드하는 법

```
NODE_PATH=<npm 전역 경로> node build_ocrt_paper.js
```

docx-js v9.6.1 이 필요하다. 산출은 `OCRT_reference_paper_draft.docx` 다. 그림을 바꾸면 `paper_figs/` 에
복사해야 빌드에 반영된다.

### 이번 세션에서 한 일

순수해수 검증 그림(옛 Figure 4)을 삭제하고 이후 그림 번호를 다시 매겼다. 삭제한 이유는 두 가지다. 첫째, 그
그림의 자료가 2026-07-23 수중 출력 방위각 규약 보정 이전에 만들어져 지금 코드의 규약과 어긋난다. 둘째,
구성모델 비교의 CDOM 계열이 엽록소와 총부유물이 모두 0 이라 순수해수 극한을 이미 포함한다. 본문에 그 취지를
한 문장 넣었다.

### 남은 일

본문에 `[TO COMPLETE]` 가 세 곳 있다.

1. 7절 GOCI-III 응용 예시 그림 — **03_polarization 의 결과가 바로 이 자리에 들어간다**
2. 자료 공개 절의 아카이브 DOI
3. 저자 목록과 사사

그 밖에 본문과 코드가 어긋나는 지점이 세 군데 있다. 아직 손대지 않았다.

- 에어로졸 라이브러리를 두 벌(16종, 80종)로 기술했으나 패키지에는 세 번째 벌(`aerosol_ahmad2010_paper_mie`,
  A2010ver 80종)이 있다
- 표면 옵션을 `--surface coxmunk` 로 적었으나 정식 이름은 `black_fresnel_ocean` 이고 `coxmunk` 는 폐기
  예정 별칭이다
- 재현 가능 빌드 절에 `-march=native` 를 `-march=cascadelake` 로 고정한 사실이 반영되어 있지 않다

### 배경으로만 남은 미결

- 식물플랑크톤 산란 복원(현재 엽록소는 흡수만 기여)
- 구성모델 경로의 전방 피크 절단 미구현(대규모 격자 생산에 4~6 % 편향)

---

## 3. 편광 논문 (03_polarization)

### 무엇을 밝혔는가

단일산란알베도(흡수)와 입자 크기는 둘 다 세기 반사도의 분광 기울기를 움직인다. 그래서 세기만으로는 둘을
갈라낼 수 없다. 이것이 축퇴다. 편광이 이를 푸는지 보았다.

**주요 결과**

1. Ahmad 80종 라이브러리에 축퇴가 그대로 들어 있다. 미세입자 비율을 고정하고 상대습도만 바꾸면 옹스트롬
   지수는 거의 그대로인데 알베도만 크게 변한다. 미세입자 50 % 행이 가장 순수한 축퇴로, 옹스트롬 1.81~1.91 에
   알베도 0.947~0.984 다.

2. 탁도가 낮은 해수에서 근적외 두 밴드(748, 865)로 대기보정을 할 때, 관측이 둘인데 미지수가 셋(크기·흡수·
   농도)이라 풀리지 않는다. 잡음을 가정하지 않고 필요한 정확도를 거꾸로 구했다.

   | 태양천정각 | 세기 반사도에 필요 | 선형편광도에 필요 | 완화 배수 |
   |---|---|---|---|
   | 25도 | 1.5e-05 | 5.6e-03 | 381배 |
   | 50도 | 1.8e-05 | 4.9e-03 | 273배 |
   | 75도 | 8.8e-05 | 5.8e-03 | 66배 |

   세기 반사도의 현실적 정확도는 보정의 2 % 곧 약 5e-4 이므로, 세기만으로는 어느 방위각에서도 닿지 못한다.
   편광은 5e-3 이면 되고 이는 현재 기술로 닿는 값이다.

3. 방위각 의존이 뚜렷하다. **태양 반대쪽(135~180도)이 가장 관대하고**, 태양천정각 50도에서 1.1e-02 까지
   허용된다. 태양천정각 25도에서는 방위각을 어떻게 골라도 4~6e-03 에 머문다. 태양천정각이 작으면 산란각
   식에서 방위각의 지렛대 자체가 사라지기 때문이다.

4. 흑색 해면 이상화는 탁도가 낮은 해수에서 정당하다. 맑은 해수를 미지수로 넣어도 검색 불확실이 1.00배로
   거의 변하지 않는다.

### 방법을 두 번 갈아엎었다. 그 기록을 남긴다

**첫 번째 방법(폐기)** — 관측을 미지수의 일차식으로 근사하고 최적 추정으로 오차를 전파했다. 문제가 둘이었다.

- 광학두께 범위가 30배라 일차식의 잔차가 측정 잡음의 7.5배였다. 값을 신뢰할 수 없다.
- 그 결과 방위각 지도에 실제 물리가 아닌 고주파 얼룩이 생겼다. 굴곡도(2계 차분 ÷ 1계 차분)가 원자료 0.19 인데
  지표는 0.56 이었다.

**두 번째 방법(현행)** — 근사를 쓰지 않는다. 두 모델이 각자 광학두께를 자유롭게 골라 748 nm 세기가 같아지도록
맞춘 뒤, 그 지점에서 남는 865 nm 세기 차이와 선형편광도 차이를 잰다. 필요한 정확도는 그 값을 3 으로 나누면
된다. 잡음을 가정하지 않는다.

방위각 얼룩을 진단한 결과도 남겨 둔다(`figures/diag_방위각_매끄러움.png`). 원자료는 매끄럽고(9차 이상 모드
1.6 %), 두 모델의 차이에는 주기 45~120도의 실제 구조가 있으며, 얼룩은 최소제곱 맞춤 단계에서 생겼다.

### 자료 (03_polarization/data)

| 파일 | 내용 |
|---|---|
| `atm_reference.csv.gz` | 해수 없는 대기(흑색 Fresnel)의 전 각도 반사도. 에어로졸 24종 × 광학두께 7단 × 5밴드 × 태양천정각 3단 × 1,296방향 |
| `cases_thin.csv.gz` | 해수 결합 실행의 원값. 성긴 각도(관측천정각 15도, 방위각 45도 간격 48방향) |
| `metrics.csv.gz` | 유도 지표(N, S, 유효도). **첫 번째 방법으로 만든 것이라 정량값은 쓰지 말 것** |
| `runs_summary.csv` | 실행 8,045건의 상태·소요·수렴·판번호·해시 |
| `smooth.csv` | 현행 방법의 최종 산출. 3,888방향 × (필요 세기 정확도, 필요 편광 정확도) |
| `aerlib.csv` | 에어로졸 80종의 옹스트롬 지수·알베도·비대칭인자 |
| `README_reduced.md` | 축약 결과의 열 설명 |

### 순수 Rayleigh 기준 (03_polarization/rayleigh)

에어로졸이 없는 대기의 반사도다. 748·865 nm × 태양천정각 25/50/75도 = 6개 파일이며, 각 1,296방향이다.
**Rayleigh 보정 반사도를 만들려면 반드시 필요하고, 로컬에는 없다.** 이 세션에서 만든 것이므로 잃어버리면
다시 계산해야 한다.

만든 명령은 이렇다.

```
ocrt_v1.2 --sza <25|50|75> --wavelength <748|865> \
  --surface black_fresnel_ocean --wind-speed 5 --decouple-sunglint --pssa \
  --lut-vza-step 5 --lut-raa-step 5 --output-full-grid R_<band>_<sza>.csv
```

에어로졸 옵션을 아예 주지 않아야 한다. `--aod-865 0` 을 주면 오류가 난다.

### 스크립트 (03_polarization/scripts)

| 파일 | 하는 일 |
|---|---|
| `smooth_signal.py` | 현행 지표를 계산해 `smooth.csv` 를 만든다 |
| `make_final_figs.py` | F1~F3 그림 |
| `make_map.py` | F4·F5 방향 지도(면 형식) |

두 스크립트가 만드는 그림이 겹치지 않는다. 방향 지도는 `make_map.py` 만 만든다.

돌리는 곳은 축약 자료가 있는 폴더이고, `/tmp/ray/`, `/tmp/aerlib.csv`, `/tmp/smooth.csv` 경로를 쓰고
있으므로 **경로를 자기 환경에 맞게 고쳐야 한다.**

### 그림 (03_polarization/figures)

| 파일 | 내용 |
|---|---|
| `F1_축퇴격자와_대표짝.png` | 축퇴 구조와 대표 짝의 편광 |
| `F2_필요정확도.png` | 세기와 편광에 필요한 정확도(핵심 그림) |
| `F3_방위각의존.png` | 방위각별 요구 정확도 |
| `F4_방향지도.png` | 관측 방향 지도(면 형식) |
| `F5_방향지도_두기준.png` | 짝 절반 기준과 어려운 짝 기준 |
| `diag_방위각_매끄러움.png` | 방법 선택의 근거. 논문에는 안 들어간다 |

---

## 4. 실행 도구 (04_run_tools)

패키지 루트의 `build` 폴더에 두면 위치를 스스로 찾는다.

| 파일 | 하는 일 |
|---|---|
| `run_star.py` | 편광 캠페인 실행기. 단계(base, p4a~p4c, p5a~p5, p6, p6a)를 골라 돌린다 |
| `run_part1~3.bat` | 1차 실행을 세 조각으로 |
| `run_part4.bat`, `run_part4-2.bat` | 4단계 확장 |
| `run_part6.bat` | 축퇴 격자와 기준 대기 |
| `reduce_grid.py` | 원자료 10 GB 를 55 MB 로 줄인다 |
| `merge_parts.py/.bat` | 조각 결과를 합치고 빠짐을 찾는다 |
| `make_ac_set.py` | **대기보정 검증 자료 생성기(다음 작업)** |

5단계(풍속·태양천정각 확장) 배치 파일은 폐기했다. 축퇴 격자를 우선하기로 하면서 설계가 바뀌었기 때문이다.
필요하면 `run_star.py` 의 `--stage p5a` 등으로 되살릴 수 있다.

---

## 5. 다음에 할 일

### (가) 대기보정 검증 자료 만들기 — 준비 완료

`04_run_tools/make_ac_set.py` 로 돌린다.

```
python make_ac_set.py --dry-run     계획 확인
python make_ac_set.py --n 20        적은 수로 시험(약 2분)
python make_ac_set.py               본 실행(24코어 약 7시간)
```

사례 1,000개 × GOCI-II 12밴드(380~865 nm)이고, 조합마다 실행 다섯 번이 필요하다. 결과표에 기하각, Rho_TOA,
Rho_R, Rho_C, Rho_Am, Td_R_s, Td_R_v, Td_Am_s, Td_Am_v, 천저 방향 Rrs 가 들어간다.

검산 두 가지가 자동으로 붙는다. `Td_check` 는 1 에서 3e-7 이내, `Rho_C = Rho_Am + t·ρ_w` 는 상대차 1.6e-15
였다.

**자료 설명에 적어 둘 두 줄**

1. 865 nm 에서 엽록소 흡수는 0 으로 처리된다. 자료표가 850 nm 에서 끝나기 때문이며, 그 파장에서 순수해수
   흡수가 4.6 /m 로 지배적이라 영향이 없다. 쇄설물 산란은 온전히 살아 있다(b 가 다른 파장과 동일).
2. 쇄설물 후방산란 비율은 850 nm 값을 유지한다. 그 구간의 파장 변화가 0.07 % 라 오차가 0.01 % 아래다.

### (나) 편광 논문 마무리

F1~F5 를 본문에 배치하고, 위 결과를 서술한다. 남은 분석 후보는 둘이다.

- 밴드 조합 최적화. 다섯 밴드 중 어느 조합이 가장 효율적인지 보면 편광 밴드 선정 근거가 된다.
- 풍속과 태양천정각 확장. `run_star.py --stage p5a` 로 되살릴 수 있다.

### (다) 모델 기술 논문 마무리

2절의 `[TO COMPLETE]` 세 곳과 어긋난 지점 세 군데를 정리한다.

### (라) 버그 보고서 후속

`06_bugreports/` 에 둘이 있다.

1. **엽록소 파장 상한** — 상한 상수만 1100 nm 로 바뀌고 자료표가 따라오지 않았다. 850 nm 위에서 엽록소 흡수가
   조용히 0 이 된다. 지금 용도에는 무해하다고 판단했으나 기록은 남긴다.
2. **Windows 빌드** — `setenv`/`unsetenv` 가 링크에서 빠진다. `05_build/ocrt_win_compat.c` 로 우회하고 있다.
   아직 고쳐지지 않았다.

---

## 6. 알아 두어야 할 것

### 실행 환경

- CPU 는 Core Ultra 9 275HX 로 **AVX-512 가 없다.** `-march=cascadelake` 로 빌드하면 실행 자체가 실패한다.
  `05_build/build_win.bat` 은 `x86-64-v3`(AVX2) 를 기본으로 쓴다.
- 배치 파일은 반드시 순수 ASCII 에 CRLF 줄바꿈이어야 한다. 한글을 넣으면 cmd 가 바이트 위치를 잘못 계산해
  엉뚱한 줄을 명령으로 해석한다.
- PowerShell 에서는 배치 파일 앞에 `.\` 가 필요하다.
- Python 은 Anaconda 프롬프트에서 쓴다.

### 자료가 어긋나면 조용히 틀린다

실행 파일만 갈아 끼우고 `inputs` 를 그대로 두는 사고가 실제로 있었다. `run_star.py` 에 주요 자료 파일의 해시
대조를 넣어 두었으므로 `--verify` 와 본 실행 시작 때 걸러진다.

### 산란각과 선글린트 규약

- 산란각: cos Θ = −cos(태양천정각)·cos(관측천정각) − sin(태양천정각)·sin(관측천정각)·cos(상대방위각).
  Rayleigh 한계에서 선형편광도 최댓값이 정확히 90도에 오는 것으로 확인했다.
- 선글린트 기울기각: cos β = (cos 태양천정각 + cos 관측천정각) / √(2 + 2·s·v). 0 이면 잔잔한 수면에서 곧바로
  선글린트가 생기는 방향이다. 태양천정각 40도에서 관측천정각 40도·방위각 180도에 0 이 오고, 선글린트를 켠
  실행의 세기 최대점과 일치함을 확인했다.

### 상향 투과율

수출광이 있어야 정의되므로 흑색 Fresnel 실행에서는 `nan` 으로 나온다. Rayleigh 만의 상향 투과율을 얻으려면
**에어로졸 없는 결합 실행**이 따로 필요하다. `make_ac_set.py` 의 네 번째 실행이 이것이다.

---

## 7. 이관 검증 결과

격리 폴더에 풀어 실제로 재현되는지 확인했다.

- 논문 재빌드: `build_ocrt_paper.js` 로 만든 docx 의 본문(`word/document.xml`)이 동봉한 원고와 해시까지 같다.
- 편광 그림: F1~F5 와 진단 그림까지 여섯 장이 모두 바이트 단위로 같다.

편광 스크립트는 절대 경로(`/tmp/ray/` 등)를 쓰고 있으므로, 옮긴 뒤 아래처럼 상대 경로로 고치고 `data`
폴더에서 돌리면 된다.

```
sed -i "s|/tmp/aerlib.csv|../data/aerlib.csv|g; s|/tmp/ray/|../rayleigh/|g; \
        s|/tmp/smooth.csv|../data/smooth.csv|g; \
        s|/mnt/user-data/outputs/|../figures/|g" \
    ../scripts/*.py
```

논문 빌드 스크립트도 `FIG_DIR` 과 출력 경로를 자기 환경에 맞게 고쳐야 한다.

---

## 8. 재현 점검

이관받은 뒤 아래를 확인하면 환경이 맞는지 알 수 있다.

```
cd <패키지루트>\build
python run_star.py --verify
```

자료 파일 확인과 기준 사례 두 건이 차례로 나온다. 기준값은 이렇다.

| 사례 | I | Q | U |
|---|---|---|---|
| 대기만(흑색 Fresnel, 555 nm, 태양천정각 25도, 관측천정각 30도, 방위각 90도, 광학두께 0.1, r50f05v01) | 4.4991861596e-02 | −9.0983939253e-04 | −7.5031808484e-03 |
| 위에 해수 결합(엽록소 1.0, 총부유물 1.5, CDOM 0.035) | 1.2305024573e-01 | −6.3070476021e-04 | −1.1478355407e-02 |

둘 다 상대차 0 으로 통과해야 한다.
