# OCRT polarization sensitivity — stage 6 run guide (aerosol degeneracy grid)

- Date: 2026-08-03
- Version: `OCRT-v1.2-2026-08-02-KST-dtpsign-phase-diagnostic`
- Run after stages 1–4 (3,675 runs) are complete. The stage-5 batch file is withdrawn.

---

## 1. Why the design changes

The aim of this polarization study is to **separate aerosol models**. Both the single-scattering
albedo (absorption) and the particle size move the spectral slope of the intensity reflectance,
so intensity alone cannot separate absorption from size. This is the degeneracy.

The Ahmad 80-model library contains this degeneracy as it is. With the fine-mode fraction fixed
and only the relative humidity changed, the Ångström exponent stays almost the same while the
single-scattering albedo changes a lot.

| fine-mode fraction | Ångström range | SSA(865) range |
|---|---|---|
| 20 % | 0.37 | 0.0152 |
| 30 % | 0.27 | 0.0229 |
| **50 %** | **0.10** | **0.0364** |
| 80 % | 0.33 | 0.0524 |

The 50 % row is the purest degeneracy: the Ångström exponent is 1.81–1.91, practically the same,
while the albedo spans 0.947–0.984.

**The problem is that the 8 models used so far do not form this degeneracy.** They were chosen to
cover the albedo evenly, so they lie on a diagonal in the (size, absorption) plane; the two
effects move together and cannot be separated.

## 2. New grid

4 fine-mode fractions × 4 relative humidities = 16 models. **Rows are the size axis and columns are
the absorption axis**, so the two effects can be looked at separately.

| | RH 30 | RH 50 | RH 80 | RH 95 |
|---|---|---|---|---|
| **fine 20 %** | r30f20v01 | r50f20v01 | r80f20v01 | r95f20v01 |
| **fine 30 %** | r30f30v01 | r50f30v01 | r80f30v01 | r95f30v01 |
| **fine 50 %** | r30f50v01 | r50f50v01 | r80f50v01 | r95f50v01 |
| **fine 80 %** | r30f80v01 | r50f80v01 | r80f80v01 | r95f80v01 |

SSA(865) covers 0.9288–0.9890 and the Ångström exponent 0.87–2.39.

**The existing 8 models are kept.** Together they make a 24-model design, so all runs already
finished remain valid. The existing 8 models cover the albedo evenly and continue to serve the
general separability measure.

## 3. Two sub-stages

| sub-stage | what it does | new runs | 24 cores |
|---|---|---|---|
| `p6` | run the new 16 models with the coupled ocean | 1,680 | about 4.3 h |
| `p6a` | run all 24 models with the atmosphere only (no ocean) | 2,520 | about 0.4 h |

`p6a` is the **reference atmosphere without ocean** (black Fresnel surface). Subtracting it from
the coupled run under the same conditions gives exactly the ocean contribution in both intensity
and polarization. Without it one cannot judge whether the ocean polarization is negligible. An
atmosphere-only run takes 11 s, so it is very cheap.

Conditions: SZA 25·50·75°, 5 bands, wind 5 m/s, 7 AOD levels (0.02–0.60), sun glint decoupled,
pseudo-spherical approximation on.

## 4. Run order

Overwrite `run_star.py` with the new version; it contains the degeneracy grid and the surface
branch.

```
.\run_part6.bat --dry-run     should report "3675 already done"
.\run_part6.bat
```

Sub-stages can be run separately.

```
.\run_part6.bat p6a           reference atmosphere first (24 min)
.\run_part6.bat p6
```

The output folder must be the same as before. The default is `<package root>\star_run`.

## 5. Outputs

Reference-atmosphere results go to files starting with `A_`; coupled-ocean results start with `O_`.

```
A_b555_s50_mr30f50v01_a0p15.csv     atmosphere without ocean
O_b555_s50_chl1_tsm1p5_cdm0p035_mr30f50v01_a0p15.csv    coupled ocean
```

The record files are `index_p6.csv` and `index_p6a.csv`. A new `surface` column in the record
distinguishes the two kinds.

## 6. Reduction

After the runs, execute the reduction script; it now also processes the reference atmosphere.

```
python reduce_grid.py --run-dir D:\temp\star_run --out D:\temp\reduced --workers 12
```

A new file `atm_reference.csv.gz` appears: the all-angle reflectance of the atmosphere without
ocean, with columns `atm_I`, `atm_Q`, `atm_U`. Subtracting it from the coupled run gives the ocean
contribution in intensity and polarization.

Uploading only the reduced folder is enough to continue the analysis; it is expected to be about
20 MB.

## 7. Questions this data set answers

1. Does polarization separate aerosols that share the same intensity spectral slope? In the
   concept test (atmosphere without ocean, SZA 50°) the polarized reflectance of the 50 % row
   spread 9 times more than the intensity slope. The coupled runs will confirm this.
2. How much does the ocean blur that separation? Measured exactly by subtracting the reference
   atmosphere.
3. Which band and which viewing direction separate best?
4. Is the spectral slope of polarization (e.g. the 443 nm / 865 nm polarization ratio) better than
   single-band polarization? All five bands are polarization bands, so this is a usable observable.

## 8. Deferred

The extension to 3 wind speeds and 6 SZA levels is on hold. If the degeneracy break turns out to
depend strongly on the viewing direction, the SZA extension should be prioritized again. That
decision follows the stage-6 results.

---

# OCRT 편광 민감도 — 6단계 실행 안내 (에어로졸 축퇴 격자)

- 작성일: 2026-08-03
- 판본: `OCRT-v1.2-2026-08-02-KST-dtpsign-phase-diagnostic`
- 앞선 1~4단계(3,675건) 완료 뒤 이어서 돌린다. 5단계 배치 파일은 폐기한다.

---

## 1. 왜 설계를 바꾸는가

이 편광 연구의 목적은 **에어로졸 모델의 구분**이다. 그런데 단일산란알베도(흡수)와 입자 크기는 둘 다 세기
반사도의 분광 기울기를 움직인다. 그래서 세기만 보면 흡수와 크기를 갈라낼 수 없다. 이것이 축퇴다.

Ahmad 80종 라이브러리에는 이 축퇴가 그대로 들어 있다. 미세입자 비율을 고정하고 상대습도만 바꾸면
옹스트롬 지수는 거의 그대로인데 단일산란알베도만 크게 변한다.

| 미세입자 비율 | 옹스트롬 지수 폭 | 단일산란알베도(865) 폭 |
|---|---|---|
| 20 % | 0.37 | 0.0152 |
| 30 % | 0.27 | 0.0229 |
| **50 %** | **0.10** | **0.0364** |
| 80 % | 0.33 | 0.0524 |

미세입자 50 % 행이 가장 순수한 축퇴다. 옹스트롬 1.81~1.91 로 사실상 같은데 알베도는 0.947~0.984 다.

**문제는 지금까지 쓴 8종이 이 축퇴를 만들지 않는다는 점이다.** 알베도를 고르게 덮도록 골랐기 때문에
(크기, 흡수) 평면에서 대각선으로 늘어서 있다. 두 효과가 함께 움직여 서로 분리되지 않는다.

## 2. 새 격자

미세입자 비율 4수준 × 상대습도 4수준 = 16종이다. **행이 크기 축, 열이 흡수 축**이 되어 두 효과를 나눠 볼 수
있다.

| | 습도 30 | 습도 50 | 습도 80 | 습도 95 |
|---|---|---|---|---|
| **미세 20 %** | r30f20v01 | r50f20v01 | r80f20v01 | r95f20v01 |
| **미세 30 %** | r30f30v01 | r50f30v01 | r80f30v01 | r95f30v01 |
| **미세 50 %** | r30f50v01 | r50f50v01 | r80f50v01 | r95f50v01 |
| **미세 80 %** | r30f80v01 | r50f80v01 | r80f80v01 | r95f80v01 |

알베도(865)는 0.9288 에서 0.9890 까지, 옹스트롬 지수는 0.87 에서 2.39 까지 덮는다.

**기존 8종을 버리지 않는다.** 합쳐서 24종으로 설계했으므로 이미 돌린 실행이 모두 그대로 살아난다. 기존
8종은 알베도를 고르게 덮으므로 일반적인 분리력 측도에 계속 쓰인다.

## 3. 두 하위 단계

| 하위 단계 | 무엇을 하는가 | 새 실행 | 24코어 |
|---|---|---|---|
| `p6` | 새 16종을 해수 결합으로 실행 | 1,680 | 약 4.3시간 |
| `p6a` | 24종 전부를 해수 없는 대기로 실행 | 2,520 | 약 0.4시간 |

`p6a` 는 **해수가 없는 기준 대기**(흑색 Fresnel 해면)다. 같은 조건의 결합 실행에서 이 값을 빼면 해수가 더한
몫이 세기와 편광 모두에서 정확히 나온다. 이것 없이는 "해수 편광이 무시할 만한가"를 잴 수 없다. 대기만
계산이라 1건에 11초로 매우 싸다.

조건은 태양천정각 25·50·75도, 밴드 5개, 풍속 5 m/s, 광학두께 7단(0.02~0.60), 선글린트 분리, 의사구면근사
켜짐이다.

## 4. 실행 순서

`run_star.py` 를 새것으로 덮어쓴다. 축퇴 격자와 표면 분기가 그 안에 있다.

```
.\run_part6.bat --dry-run     "이미 끝난 것 3675건" 이 나와야 한다
.\run_part6.bat
```

하위 단계를 따로 돌릴 수도 있다.

```
.\run_part6.bat p6a           기준 대기만 먼저 (24분)
.\run_part6.bat p6
```

산출물 폴더는 반드시 이전과 같은 곳을 가리켜야 한다. 기본값은 `<패키지루트>\star_run` 이다.

## 5. 산출물

기준 대기 결과는 `A_` 로 시작하는 파일에 들어간다. 해수 결합은 `O_` 다.

```
A_b555_s50_mr30f50v01_a0p15.csv     해수 없는 기준 대기
O_b555_s50_chl1_tsm1p5_cdm0p035_mr30f50v01_a0p15.csv    해수 결합
```

기록 파일은 `index_p6.csv`, `index_p6a.csv` 로 따로 쌓인다. 기록에 `surface` 열이 새로 생겨 두 종류를
구분한다.

## 6. 축약

끝나면 축약 스크립트를 돌린다. 기준 대기를 함께 처리하도록 고쳐 두었다.

```
python reduce_grid.py --run-dir D:\temp\star_run --out D:\temp\reduced --workers 12
```

`atm_reference.csv.gz` 가 새로 생긴다. 해수 없는 대기의 전 각도 반사도이며 열은 `atm_I`, `atm_Q`,
`atm_U` 다. 결합 실행에서 이 값을 빼면 해수의 몫이 세기와 편광 모두에서 나온다.

축약 폴더만 올리면 분석을 이어갈 수 있다. 전체가 20 MB 안팎으로 예상된다.

## 7. 이 자료로 답할 질문

1. 같은 세기 분광 기울기를 가진 에어로졸을 편광이 갈라내는가. 개념 검증(해수 없는 대기, 태양천정각 50도)
   에서 미세입자 50 % 행의 편광 반사도가 세기 기울기보다 9배 크게 벌어졌다. 결합 실행으로 확정한다.
2. 해수가 그 분리를 얼마나 흐리는가. 기준 대기 차감으로 정확히 잰다.
3. 어느 밴드와 어느 관측 방향에서 가장 잘 갈라지는가.
4. 편광의 분광 기울기(예: 443 nm 대 865 nm 편광 비)가 단일 밴드 편광보다 나은가. 다섯 밴드가 모두 편광
   밴드이므로 실제로 쓸 수 있는 관측량이다.

## 8. 뒤로 미룬 것

풍속 3단과 태양천정각 6단 확장은 보류했다. 축퇴 해소가 관측 방향에 크게 의존하는 것으로 확인되면
태양천정각 확장을 다시 우선순위에 올리는 것이 맞다. 그 판단은 6단계 결과를 보고 한다.
