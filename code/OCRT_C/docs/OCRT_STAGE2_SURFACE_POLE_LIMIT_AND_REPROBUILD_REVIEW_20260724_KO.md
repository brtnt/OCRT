# OCRT–OSOAA 2단계 수면 극한 회전 및 재현 빌드 업데이트 검토 보고서

문서 기준일: 2026-07-24  
검토 대상: `files (41).zip`  
적용 전 기준선: `OCRT-v1.2-2026-07-24-KST-stage2-interface-fix123-dtp-uclosure`  
적용 후 기준선: `OCRT-v1.2-2026-07-24-KST-stage2-interface-fix123-dtp-uclosure-pole-limit`

## 1. 종합 판정

첨부 업데이트의 두 핵심 제안은 타당하다.

1. rough-interface Mueller kernel에서 입사 또는 출사 방향이 정확히 연직일 때 일반 자오면 회전식이 0/0이 되므로 해석 극한을 명시적으로 사용하여야 한다.
2. 릴리스 빌드에 `-march=native`를 사용하면 실행 호스트의 CPU에 따라 코드 생성과 부동소수 연산순서가 달라지므로 고정 ISA를 사용하여야 한다.

이에 따라 다음을 코드에 반영하였다.

- 네 rough-interface kernel의 exact-pole meridian limit
- 극 연속성과 연직 spin-2 공변성 회귀시험
- 기본 릴리스 ISA를 `cascadelake`로 고정
- 독립 두 번 빌드의 binary SHA-256 재현성 시험
- pandas를 통하지 않는 CSV binary64 비트 비교 도구
- 누적 작업지시서 개정 4

대기 물리모듈과 1단계 대기 하네스, public water RAA/U 출력, FIX1–FIX4는 변경하지 않았다.

## 2. 첨부 패키지 무결성

ZIP에 실제 포함된 다음 네 파일은 첨부 manifest의 SHA-256과 일치하였다.

```text
NOTE_double_pole_and_build_2026-07-24.md
PATCH_surface_POLE_2026-07-24.diff
surface.c
surface_double_pole_perturbed.c
```

반면 manifest에 기재된 다음 다섯 파일은 ZIP에 포함되지 않았다.

```text
ADDENDUM2_atmval_pole_2026-07-24.md
atmval_v4_2026-07-24.csv
atmval_v4_2026-07-24.png
polecont.c
probe_nadir_pole_constraint.py
```

따라서 첨부의 종단 정확도 수치 및 224조건 전수 교란시험은 원자료로 직접 재현하지 않았다. 대신 포함된 패치의 수학적 극한과 별도의 C 회귀시험을 독립 구성하였다.

## 3. 물리·수학 검토

일반 자오면 회전의 한쪽 각도는 다음 구조를 갖는다.

\[
\cos\sigma
=
\frac{\mu_{\mathrm{other}}-\mu_{\mathrm{self}}\cos\Theta}
{\sin\Theta\,\sin\theta_{\mathrm{self}}},
\qquad
\sin\sigma
=
\frac{\sin\theta_{\mathrm{other}}\sin\phi}{\sin\Theta}.
\]

`sin(theta_self) -> 0`에서 분자와 분모가 동시에 0으로 수렴한다. 분모에만 `1e-12` 바닥값을 적용하면 정확한 극에서 회전이 identity에 가까워지며 실제 일방 극한과 불연속이 된다.

방향벡터를 일차 전개하면 다음 극한을 얻는다.

\[
\cos\sigma\rightarrow -\operatorname{sign}(\mu_{\mathrm{self}})\cos\phi,
\qquad
\sin\sigma\rightarrow \sin\phi.
\]

첨부 패치는 이 식을 사용하며 물리적으로 타당하다. `OCRT_SURF_POLE_EPS=1e-8`은 대략 `sqrt(DBL_EPSILON)` 규모로, 원식의 상쇄오차와 극한식의 절단오차가 교차하는 합리적 문턱이다. 생산 격자에서는 연직 특수절점 또는 충분히 큰 비극 사인값만 나타나므로 문턱 미세조정에 민감하지 않다.

## 4. 코드 변경

### 4.1 `src/shared/surface.c`

다음 네 함수의 입사측·출사측 회전 블록에 극한 분기를 적용하였다.

```text
surface_R_coxmunk_trig
surface_T_coxmunk_trig
surface_R_ww_coxmunk_trig
surface_T_aw_coxmunk_trig
```

각 방향에서 `sin(theta_self) <= 1e-8`이면 해석 극한을 사용하고, 그 외에는 기존 일반식과 정규화를 그대로 사용한다.

공기→물 경로에서는 기존 FIX1+FIX2의 signed propagation cosine과 pi-shifted relative azimuth를 그대로 사용하였다. FIX3 incoming-U Fourier 열과 FIX4 diffuse-top incoming-U 열은 변경하지 않았다.

### 4.2 빌드 스크립트

다음 파일에서 `-march=native`를 제거하였다.

```text
scripts/build_release_v1.1.sh
scripts/build_release_v1.2.sh
scripts/regression_ext_bottom_mode_bounds.sh
```

기본값은 다음과 같다.

```bash
OCRT_MARCH=${OCRT_MARCH:-cascadelake}
```

다른 ISA가 필요한 경우 명시적으로 override할 수 있으나, 그 경우 승인 바이너리 해시를 재사용하지 않는다.

### 4.3 부수적인 회귀스크립트 수정

`regression_raa_convention.sh`가 `tests/test_raa_convention.c`를 링크할 때 `src/rt_fourier.c`를 포함하지 않아 현재 소스에서 링크 실패하던 문제를 수정하였다. 이는 물리코드 변경이 아니다.

### 4.4 CSV 비트 비교

`scripts/compare_csv_bitexact.py`를 추가하였다. numeric cell은 pandas가 아니라 Python 내장 `float()`로 파싱하여 IEEE-754 binary64 비트열을 비교한다. `--raw`를 사용하면 원시 파일 바이트를 비교한다.

## 5. 독립 모듈 회귀시험

### 5.1 정확한 극과 일방 극한의 연속성

수정 전에는 네 kernel의 입사/출사 극 8개 경우 모두에서 상대 불연속이 약 1.5였다.

```text
R_air 1.499990
T_wa  1.499951–1.499953
R_ww  1.499995
T_aw  1.500047–1.500049
```

수정 후 결과는 다음과 같다.

```text
R_air 9.43e-6
T_wa  4.67e-5–4.92e-5
R_ww  5.31e-6
T_aw  4.67e-5–4.92e-5
```

모든 경우가 회귀 기준 `1e-4` 이내이다.

### 5.2 연직 spin-2 공변성

수정 전에는 극에서 Q/U의 물리적 m=2 진폭이 거의 소실되고 Q에 큰 m=0 누출이 발생하였다. 수정 후에는 다음과 같다.

```text
kernel  spin-2 amplitude   m=0 leakage/amplitude  paired-mode mismatch
R_air   2.163445e-05       8.31e-17               3.13e-16
T_wa    7.897087e-04       7.88e-17               8.24e-16
R_ww    1.026800e-03       6.38e-17               0
T_aw    2.901195e-02       5.91e-17               1.56e-15
```

이는 기준코드의 종단값이 아니라 회전 연산자의 대칭성 자체를 검증한다.

### 5.3 기존 폐합 모듈의 비회귀

다음 회귀시험을 모두 통과하였다.

```text
FIX123 air-to-water signed rotation and U-column regression: PASS
FIX4 pure-Rayleigh source/operator closure: 2.776e-17
FIX4 Rayleigh–aerosol mixed closure: 5.551e-17
public water RAA/U regression: PASS
full-grid/single-geometry parity: PASS
wind=0 branch: PASS
external-bottom mode bounds: PASS
RAA glint anchors: unchanged
```

따라서 이미 폐합된 모듈을 건드리지 않았다는 조건이 충족된다.

## 6. 양쪽 극 축퇴

입사와 출사 방향이 동시에 연직이면 상대방위각 기준면은 본질적으로 정의되지 않는다. 생산 격자에서 해당 방향쌍은 삽입 special node이며 경계 적분에서 zero quadrature weight와 곱해진다.

첨부의 `surface_double_pole_perturbed.c`를 사용하여 양쪽 극의 Q/U 회전 블록을 identity로 크게 교란한 바이너리를 별도 빌드하였다. 본 세션에서는 다음 독립 범위를 실행하였다.

```text
pure water
SZA = 0 deg
bands = 412, 443, 490, 555, 660, 865 nm
24 full-grid rows per band
```

6개 출력 CSV가 정본과 모두 byte-identical하였다. 이는 양쪽 극의 회전 모호성이 해당 범위의 종단출력에 기여하지 않는다는 독립 근거이다.

첨부가 기재한 224조건 전수검사는 누락된 원자료 때문에 전부 반복하지 않았으며, 본 보고서는 6개 파장 독립시험까지만 주장한다.

## 7. 동일 아키텍처 수정 전·후 종단 영향

수정 전 FIX123+FIX4와 수정 후 pole-limit 바이너리를 모두 동일 compiler와 `-march=cascadelake`로 빌드하였다.

SZA 40°, 정확한 VZA 0°, 순수해수 6개 파장과 8개 RAA에서 다음을 확인하였다.

```text
TOA rho_I       : 정확히 불변
TOA rho_Q/rho_U : exact-pole 기준면 교정으로 변화
Rrs I/Q/U       : 정확히 불변
rrs I/Q/U       : 정확히 불변
```

443 nm에서 최대 절대변화는 다음과 같다.

```text
max |delta TOA rho_Q| = 4.48e-4
max |delta TOA rho_U| = 3.78e-4
```

412–865 nm 전체에서는 TOA Q/U 변화가 파장이 길어질수록 감소한다. 이 결과는 패치가 exact-nadir surface-coupled polarization을 교정하지만 기존 해양 Rrs/rrs 모듈을 종단오차 조정 목적으로 변경하지 않았음을 보여준다.

## 8. 72개 실행 및 Rrs/rrs 회귀 대조

기존 72개 manifest와 대기 설정을 그대로 사용하였다.

```text
72/72 physical runs successful
24 rows per run
water_converged = 1 for all runs
maximum water order = 29
600/600 matched regression cells
```

대기 optical depth, depolarization, layer/Fourier/SOS 옵션 및 gas/aerosol off 설정은 변경하지 않았다.

기존 OSOAA 600셀 표는 최신 수심 정정 이전 자료이므로 최종 물리 정확도 기준으로 사용하지 않고 코드 비회귀 대조로만 사용하였다. pole-limit 적용 전후의 Rrs/rrs 지표는 표시 자릿수에서 동일하다.

순수해수 150셀의 회귀 지표는 다음과 같다.

```text
Rrs I MAPE             2.276658%
Rrs RMS(delta Q / I)   3.340291%
Rrs RMS(delta U / I)   2.958578%

rrs I MAPE             2.137283%
rrs RMS(delta Q / I)   3.149931%
rrs RMS(delta U / I)   2.823816%
```

이 수치는 pole-limit 패치의 정확도 개선량이 아니라 동일 기준표에 대한 비회귀 상태를 나타낸다.

## 9. 빌드 재현성

호스트 의존 pre-patch 바이너리와 fixed-ISA pre-patch 바이너리의 SHA-256이 달랐다.

```text
old prebuilt/native provenance: 7e5b2938e3529cc068a70cd11514dfa22cbbc3563aa5ec8a910ce4a6aa258fd7
same source cascadelake build:   3ce8efee2983ae77774dfa04c42b4bec76cad07114bb128439eea3027a209e38
```

pole-limit 최종 소스를 동일 조건에서 두 번 독립 빌드한 결과는 모두 다음 해시였다.

```text
6c7f4b61d311323374233e0be9643541570d3b15289457ca9336d14d8f12ffde
```

따라서 고정 ISA는 현재 동일 toolchain에서 재현된다. 단, compiler/linker 버전까지 다른 환경의 cross-toolchain byte identity를 보장하는 것은 아니다.

`cascadelake` 명령어를 지원하지 않는 구형 CPU에서는 명시적으로 더 낮은 ISA를 사용해야 하며, 그 빌드는 별도 해시와 별도 검증기준을 가져야 한다.

## 10. 계산시간

443 nm 순수해수, SZA 40°, full-grid 24행, `OMP_NUM_THREADS=1` 조건에서 각 variant를 1회 예열한 뒤 5회 교차 측정하였다.

```text
before mean   = 4.831627 s
after mean    = 4.856947 s
mean change   = +0.5241%

before median = 4.819624 s
after median  = 4.873949 s
median change = +1.1272%
```

변화는 5% 병목조사 문턱보다 충분히 작다. 새 반복문, 적분, 메모리 할당 또는 phase-function 평가를 추가하지 않았으므로 별도 프로파일링은 수행하지 않았다.

## 11. 최종 코드 상태

최종 버전:

```text
OCRT-v1.2-2026-07-24-KST-stage2-interface-fix123-dtp-uclosure-pole-limit
```

승인 바이너리 SHA-256:

```text
6c7f4b61d311323374233e0be9643541570d3b15289457ca9336d14d8f12ffde
```

동결 모듈:

- 1단계 대기 모듈 및 대기 하네스
- public water RAA/U reconstruction
- FIX1+FIX2+FIX3 TAW kernel
- FIX4 diffuse-top incoming-U source
- 네 rough-interface kernel의 exact-pole meridian limit
- fixed-ISA release build contract

계속 미해결인 항목:

- corrected-depth OSOAA SZA 80° `rrs(0−)`
- corrected-depth OSOAA `Rrs(0+)` 두 extraction branch
- rough direct-solar boundary field와 in-water solver의 분리 폐합
- water-side `R_ww` 전 성분 독립 감사
- 동일 수심·층경계·n_mu·m_max·SOS tolerance의 최종 1% gate

이번 패치는 위 미해결 물리정합 문제를 해결했다고 주장하지 않는다.
