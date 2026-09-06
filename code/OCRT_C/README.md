# OCRT C — coupled atmosphere–ocean vector radiative transfer (reference implementation)

Current note (2026-09-05): this tree is OCRT C v1.11.1. It is platform-neutral: no operating
system is a required environment or an official performance gate. See
`README_FINAL_PLATFORM_NEUTRAL_2026-08-26_KO.md` for the platform policy, and
`validation/ipss_2026-09-05/OCRT_IPSS_REPLACEMENT_RECORD_2026-09-05.md` for the spherical-shell
correction introduced in v1.11.

## Validated release lineage

- Current: v1.11.1 (2026-09-05) — IPSS spherical-shell correction, surface-anchored viewing
  angle, phase-weighted κ.
- Base: `OCRT_C_FINAL_PLATFORM_NEUTRAL_PERF_2026-08-26`, derived from
  `OCRT-v1.2-2026-08-16-KST-mie-fr631-direct-truncation`.

### Mie FR631 direct-particle release (2026-08-16)

- The active aerosol catalog (176 models), the EAP catalog (17 models) and the AHN mineral
  catalog (4 models) use one fixed grid of 631 scattering angles (FR631) wherever source data
  exist or can be remapped without loss.
- The particle phase-matrix elements P11, P12 and P33 are evaluated with one common piecewise
  linear weight in scattering angle. Cubic interpolation and the L = 200 moment reconstruction are
  not used by the production kernel.
- The direct vector Fourier kernel is cached in worker-private memory. One case performs one
  native all-view solve and reconstructs every requested VZA × RAA cell; per-cell re-solves are
  not allowed in production.
- The hydrosol forward-peak truncation option uses `mu1 = 0.85`, `mu2 = 0.92` and
  `A_TRONCA >= 0.1`. The residual phase remains a direct theta-linear table, and the transport
  scattering coefficient is scaled as `b_eff = b (1 − A/2)`.
- When the truncation criterion is not met, the truncation step changes nothing: phase arrays,
  metadata, scalar IOPs and RT output stay identical.
- A no-atmosphere case (pressure zero) has its own one-water-solve all-view path; it does not fall
  back to per-cell replay.
- Batch rows are independent: case-local atmospheric all-view state is reset at every case
  boundary.

This source also includes:

- worker-private caching of parsed Mie models and of complete direct-vector Fourier kernels;
- mandatory native all-view execution (one case gives every requested VZA × RAA cell);
- independent batch-row state, including an explicit reset of the atmospheric S7/S7b caches;
- fixed-grid theta-linear direct-particle kernel reuse with cached Fourier trigonometric tables;
- exact skipping of zero-quadrature SOS columns, with a switch to disable it for byte regression;
- RAA-invariant Beer-factor hoisting in the native coupled angular LUT;
- the stage-2 water-output RAA convention fix: `Rrs(0+)` and `rrs(0-)` use the same public RAA and
  the same positive-sine Stokes-U reconstruction as the atmospheric output;
- explicit separation of the public output reconstruction from the local single-scatter `pi-RAA`
  propagation geometry;
- regression tests for the water RAA convention, full-grid/single-geometry consistency and the
  wind-zero branch;
- the canonical surface name `--surface black_fresnel_ocean` (rough Fresnel interface over a
  black ocean); `--surface coxmunk` remains as a deprecated alias;
- explicit slope-variance naming `--sigma-model ocrt-floor|nakajima-tanaka`, with the default
  OCRT low-wind floor law documented separately from the surface boundary condition;
- a water Fourier external-bottom-source mode-bound and shape fix;
- shape-preserving PCHIP wavelength interpolation of Mie P11/P12/P33;
- an RT-derived upward transmittance output with an explicit validity flag;
- a physically corrected Kd(0−) that uses the unscattered transmitted skylight at the surface
  and at the first water level;
- full-grid CSV output of constituent IOPs, including phytoplankton-only `a_chl`;
- 16 SnF, 80 Ahmad/AccuRT and 80 A2010ver aerosol Mie models in the full package.

### Stage-2 air-to-water interface correction (2026-07-24)

- Uses signed photon-propagation cosines and the matching pi-shifted relative azimuth in the
  rough air-to-water Stokes rotation.
- Stores the m > 0 incoming-U Fourier-column sign in the TAW operator, so the downstream
  contraction stays an ordinary 3 × 3 matrix–vector product.
- Adds independent polar-limit and Fourier-column regression tests.
- Closes the diffuse-top primary-source incoming-U column against the production SOS operator
  for every mode, node and basis direction tested. The three molecular U-input coefficients in
  `ocrt_add_diffuse_top_primary` carry the sign required by the existing U-output wrapper.
- Adds a reference-free source/operator contraction regression for pure-Rayleigh and mixed
  Rayleigh–aerosol media.
- The correction clearly reduces the high-SZA Q/U residuals. It also exposes separate low-SZA
  intensity residuals that are still under study. It is not an empirical fit to any reference
  model.

### Exact-pole surface rotation and reproducible build (2026-07-24)

- Applies the analytical meridian-rotation limit when the incident or outgoing direction is
  exactly vertical, in all four rough-interface Mueller kernels (`R_air`, `T_wa`, `R_ww`, `T_aw`).
- Keeps the validated non-pole equations, the FIX1–FIX4 conventions and the atmospheric harness.
- Adds exact-pole / one-sided-limit continuity and vertical spin-2 covariance regressions.
- Documents that the simultaneous double-pole ambiguity is multiplied by a zero quadrature weight
  in production boundary sums; a deliberately perturbed double-pole rotation leaves the audited
  full-grid outputs unchanged.
- Replaces host-dependent `-march=native` builds with a fixed default `-march=cascadelake`. Set
  `OCRT_MARCH=<target>` only when a different documented CPU baseline is required; changing it
  changes the binary fingerprint and may change last-bit floating-point results.

### Water-to-air coupling interpolation clamp correction (2026-07-25)

- Removes the silent 64-direction cap in the water-to-air coupling interpolation.
- Keeps the appended zero-weight view, solar and nadir directions used by the water solver; these
  directions must take part in the sorted interpolation grid.
- Uses a 320-entry stack workspace for the validated solver range and a fail-loud heap fallback
  for larger grids.
- Restores the exact nadir spin-2 constraint at `n_mu_water = 64` and the LUT default
  `n_mu_water = 96`; `rrs(0-)` is unchanged because the defect was downstream of the in-water
  solve.
- Adds regressions for `n_mu_water = 48/64/96`, forbidden nadir `m = 0` polarization,
  full-grid/single-geometry consistency and the no-trigger byte-identical path.
- Does not modify the atmospheric solver, the atmospheric harness, the FIX1–FIX4 interface
  kernels, the exact-pole equations or the in-water SOS physics.

## Build

```bash
./scripts/build_release_v1.2.sh build/ocrt
```

The default target is `cascadelake` for deterministic code generation on the validated toolchain.
A reproducibility claim requires the same source, compiler and linker versions, flags,
environment and `OCRT_MARCH`. For another explicit target:

```bash
OCRT_MARCH=x86-64-v3 ./scripts/build_release_v1.2.sh build/ocrt
```

For bit-level CSV comparisons use raw byte comparison or `scripts/compare_csv_bitexact.py`.
Values already parsed through `pandas` must not be the source of a bit-exact assertion.

## Verification

Run from the full package root:

```bash
./verify_package.sh
./verify_package.sh --full-rebuild
```

Spherical-correction gates (no Mie tables needed): `bash scripts/run_ipss_gates.sh`.

## Repository and distribution

Canonical repository: https://github.com/brtnt/OCRT

Session deliverables are provided directly in the session by default. GitHub upload is performed
only when the user explicitly requests it.

## License

Academic and non-commercial use only. See `LICENSE`. Commercial licensing inquiries:
brtnt@kiost.ac.kr.

## EAP 17-species phase generator status

The 17-species coated-sphere EAP generator and its FR631 P11/P12/P33 files are included and pass
the active-file physical-contract gate. Species-specific phytoplankton scattering stays disabled
in the constituent production model until the final end-to-end scientific acceptance gate for the
EAP closure is passed. It is never selected silently by a chlorophyll input.

- Chl without `--ocrt-phyto-group`: the frozen validated phytoplankton absorption table is used;
  phytoplankton `b = bb = 0`; detritus scattering remains.
- Chl with an explicit `--ocrt-phyto-group`: fail-loud, exit code 2.
- The 17 EAP files remain available for generator validation and for explicit
  fixed-bulk/direct-kernel diagnostics.
- The direct-kernel and forward-truncation implementation is validated on its own; enabling a
  species in the constituent closure needs a separate scientific-model decision.

## Direct particle kernel, truncation and spectral scope

The production particle-phase path is the direct theta-linear vector kernel. The legacy moment
representation is kept only for controlled diagnostics:

- default / direct: `OCRT_WATER_PARTICLE_KERNEL=direct`
- diagnostic legacy comparison: `OCRT_WATER_PARTICLE_KERNEL=moment`
- azimuth quadrature used while building a cold direct kernel: `OCRT_WATER_VALUE_NPHI=N`

Once a complete Fourier kernel is built, it is cached in worker-private memory and reused for the
case. The requested output RAA grid does not trigger repeated SOS solves.

The active AHN mineral files were remapped to FR631 while keeping every historical 0.5-degree node
at the common phase wavelengths. The temporary 330-nm mineral phase contract is a 350-nm endpoint
hold, because a canonical AHN microphysical generator was not present in the supplied package.
Source-informed regenerated candidates were kept under validation and were not installed.

The active legacy detritus file stays on its original variable angle grid and is evaluated by the
same theta-linear direct evaluator. The rejected 330–1100 nm detritus candidate stays out of
production.

At exactly VZA = 0 the polarized direct-kernel output extraction keeps the documented singularity.
Use VZA = 0.001 degree for polarized diagnostic output at nominal nadir; no silent geometry
substitution is applied.

The OCRT constituent spectral contract is 330–1100 nm, subject to the explicit data contract of
each selected component. Pure-water optical properties are tabulated over the wider Z09 interval;
the separate temperature-correction table must also cover the requested wavelength when the
temperature differs from 20 °C.

---

# OCRT C — 대기–해양 결합 벡터 복사전달 코드 (참조 구현)

현재 안내(2026-09-05): 이 트리는 OCRT C v1.11.1 이다. 플랫폼 중립이며, 특정 운영체제는 필수 실행환경도
공식 성능 게이트도 아니다. 플랫폼 정책은 `README_FINAL_PLATFORM_NEUTRAL_2026-08-26_KO.md`, v1.11 에서
도입한 구면 보정은 `validation/ipss_2026-09-05/OCRT_IPSS_REPLACEMENT_RECORD_2026-09-05.md` 를 참조한다.

## 검증된 릴리스 계보

- 현재: v1.11.1 (2026-09-05) — IPSS 구면 보정, 지표 기준 관측각, 위상함수 가중 κ.
- 기준: `OCRT_C_FINAL_PLATFORM_NEUTRAL_PERF_2026-08-26`, 그 이전은
  `OCRT-v1.2-2026-08-16-KST-mie-fr631-direct-truncation`.

### Mie FR631 직접 입자 커널 릴리스 (2026-08-16)

- 활성 에어로졸 목록(176 모델), EAP 목록(17 모델), AHN 광물 목록(4 모델)은 원자료가 있거나 손실 없이
  재배치할 수 있는 경우 모두 고정 산란각 631점 격자(FR631)를 쓴다.
- 입자 위상행렬 원소 P11, P12, P33 은 산란각에 대한 공통 구간선형 가중치 하나로 계산한다. 3차 보간과
  L = 200 모먼트 복원은 생산 커널에서 쓰지 않는다.
- 직접 벡터 푸리에 커널은 worker 전용 메모리에 캐시한다. 케이스 하나가 native all-view 해 하나를 풀고
  요청된 모든 VZA × RAA 셀을 복원한다. 셀별 재해석은 생산에서 금지한다.
- 하이드로졸 전방 피크 절단 옵션은 `mu1 = 0.85`, `mu2 = 0.92`, `A_TRONCA >= 0.1` 을 쓴다. 잔여 위상은
  직접 theta-linear 표로 유지하고, 수송 산란계수는 `b_eff = b (1 − A/2)` 로 조정한다.
- 절단 기준을 만족하지 않으면 절단 단계는 아무것도 바꾸지 않는다: 위상 배열, 메타데이터, 스칼라 IOP,
  RT 출력이 그대로다.
- 대기 없는 케이스(기압 0)는 별도의 one-water-solve all-view 경로를 가지며, 셀별 재실행으로 돌아가지 않는다.
- 배치 행은 서로 독립이다: 케이스 국소 대기 all-view 상태를 케이스 경계마다 초기화한다.

이 소스에는 다음도 포함된다.

- 파싱된 Mie 모델과 완전한 직접 벡터 푸리에 커널의 worker 전용 캐시;
- 필수 native all-view 실행(케이스 하나가 요청된 모든 VZA × RAA 셀을 산출);
- 독립 배치 행 상태(대기 S7/S7b 캐시의 명시적 초기화 포함);
- 고정 격자 theta-linear 직접 입자 커널 재사용과 푸리에 삼각함수 표 캐시;
- 가중치 0 인 SOS 구적 열의 정확한 건너뛰기(바이트 회귀용 비활성 스위치 포함);
- native 결합 각도 LUT 에서 RAA 에 무관한 Beer 인자 끌어올리기;
- stage-2 수중 출력 RAA 규약 수정: `Rrs(0+)` 와 `rrs(0-)` 가 대기 출력과 같은 공개 RAA 와 양의 사인
  Stokes-U 복원을 쓴다;
- 공개 출력 복원과 국소 단일산란 `pi-RAA` 전파 기하의 명시적 분리;
- 수중 RAA 규약, 전체격자/단일기하 일관성, 풍속 0 분기의 회귀 테스트;
- 표준 표면 이름 `--surface black_fresnel_ocean`(흑색 해양 위의 거친 Fresnel 계면). `--surface coxmunk` 는
  폐기 예정 별칭으로 남긴다;
- 경사분산 모델 이름 `--sigma-model ocrt-floor|nakajima-tanaka` 명시. 기본값인 OCRT 저풍속 하한식은
  표면 경계조건과 별도로 문서화한다;
- 수중 푸리에 외부 바닥원 모드 상한과 형상 수정;
- Mie P11/P12/P33 의 형상 보존 PCHIP 파장 보간;
- 유효성 플래그가 붙은 RT 기반 상향 투과율 출력;
- 표면과 첫 수중 준위의 비산란 투과 하늘광을 쓰는 물리적으로 수정된 Kd(0−);
- 식물플랑크톤 단독 `a_chl` 을 포함한 성분 IOP 의 전체격자 CSV 출력;
- 전체 패키지의 에어로졸 Mie 모델: SnF 16, Ahmad/AccuRT 80, A2010ver 80.

### Stage-2 공기→물 계면 보정 (2026-07-24)

- 거친 공기→물 Stokes 회전에서 부호 있는 광자 전파 코사인과 그에 맞는 π 이동 상대방위각을 쓴다.
- m > 0 입사-U 푸리에 열의 부호를 TAW 연산자에 저장해, 이후 축약이 보통의 3 × 3 행렬–벡터 곱으로 유지된다.
- 독립적인 극한(pole) 회귀와 푸리에 열 회귀를 추가한다.
- 확산 상단 1차 원천의 입사-U 열을 시험한 모든 모드·절점·기저 방향에서 생산 SOS 연산자와 닫는다.
  `ocrt_add_diffuse_top_primary` 의 분자 U 입력 계수 세 개는 기존 U 출력 래퍼가 요구하는 부호를 가진다.
- 순수 Rayleigh 와 Rayleigh–에어로졸 혼합 매질에 대해 기준값 없는 원천/연산자 축약 회귀를 추가한다.
- 이 보정은 고 SZA 의 Q/U 잔차를 크게 줄이고, 아직 조사 중인 저 SZA 세기 잔차를 별도로 드러낸다.
  어떤 참조 모델에도 경험적으로 맞춘 것이 아니다.

### 정확 극점 표면 회전과 재현 가능한 빌드 (2026-07-24)

- 입사 또는 출사 방향이 정확히 연직일 때 네 개의 거친 계면 Mueller 커널(`R_air`, `T_wa`, `R_ww`, `T_aw`)
  모두에서 해석적 자오선 회전 극한을 적용한다.
- 검증된 비극점 식, FIX1–FIX4 규약, 대기 하니스를 유지한다.
- 정확 극점 / 한쪽 극한 연속성과 연직 spin-2 공변성 회귀를 추가한다.
- 동시 이중 극점의 모호성은 생산 경계 합에서 가중치 0 과 곱해짐을 문서화한다. 이중 극점 회전을 일부러
  흔들어도 감사된 전체격자 출력은 변하지 않는다.
- 호스트 의존 `-march=native` 빌드를 고정 기본값 `-march=cascadelake` 로 바꾼다. 다른 CPU 기준선이
  필요할 때만 `OCRT_MARCH=<target>` 를 설정한다. 이를 바꾸면 바이너리 지문이 바뀌고 부동소수 마지막
  비트가 달라질 수 있다.

### 물→공기 결합 보간 상한 수정 (2026-07-25)

- 물→공기 결합 보간 경로의 숨은 64방향 상한을 제거한다.
- 수중 솔버가 쓰는 가중치 0 의 관측·태양·천저 방향을 유지한다. 이 방향들은 정렬된 보간 격자에 들어가야 한다.
- 검증된 솔버 범위에는 320개 스택 작업공간을, 더 큰 격자에는 fail-loud 힙 대체를 쓴다.
- `n_mu_water = 64` 와 LUT 기본값 `n_mu_water = 96` 에서 정확한 천저 spin-2 제약을 복원한다.
  결함이 수중 해의 하류에 있었으므로 `rrs(0-)` 는 변하지 않는다.
- `n_mu_water = 48/64/96`, 금지된 천저 `m = 0` 편광, 전체격자/단일기하 일관성, 비발동 바이트 동일 경로의
  회귀를 추가한다.
- 대기 솔버, 대기 하니스, FIX1–FIX4 계면 커널, 정확 극점 식, 수중 SOS 물리는 수정하지 않는다.

## 빌드

```bash
./scripts/build_release_v1.2.sh build/ocrt
```

기본 대상은 검증된 툴체인에서 결정적 코드 생성을 위한 `cascadelake` 다. 재현성을 주장하려면 같은 소스,
컴파일러·링커 버전, 플래그, 환경, `OCRT_MARCH` 가 필요하다. 다른 대상을 명시하려면:

```bash
OCRT_MARCH=x86-64-v3 ./scripts/build_release_v1.2.sh build/ocrt
```

비트 수준 CSV 비교에는 원시 바이트 비교 또는 `scripts/compare_csv_bitexact.py` 를 쓴다. `pandas` 로
이미 파싱한 값은 비트 동일 판정의 근거로 쓰지 않는다.

## 검증

전체 패키지 루트에서 실행한다.

```bash
./verify_package.sh
./verify_package.sh --full-rebuild
```

구면 보정 게이트(Mie 표 불필요): `bash scripts/run_ipss_gates.sh`.

## 저장소와 배포

정본 저장소: https://github.com/brtnt/OCRT

세션 산출물은 기본적으로 세션 안에서 직접 전달한다. GitHub 업로드는 사용자가 명시적으로 요청할 때만 한다.

## 라이선스

학술·비상업 이용만 허용한다. `LICENSE` 참조. 상업 라이선스 문의: brtnt@kiost.ac.kr.

## EAP 17종 위상 생성기 상태

17종 코팅 구형 EAP 생성기와 그 FR631 P11/P12/P33 파일이 포함되어 있으며 활성 파일 물리 계약 게이트를
통과한다. 종별 식물플랑크톤 산란은 EAP 폐합에 대한 최종 end-to-end 과학 승인 게이트를 통과할 때까지
성분 생산 모델에서 비활성으로 둔다. 엽록소 입력만으로 조용히 선택되는 일은 없다.

- `--ocrt-phyto-group` 없이 Chl 을 주면: 동결된 검증 식물플랑크톤 흡수표를 쓰고, 식물플랑크톤
  `b = bb = 0`, 쇄설물 산란은 유지한다.
- `--ocrt-phyto-group` 을 명시하면: fail-loud, 종료 코드 2.
- 17개 EAP 파일은 생성기 검증과 명시적 fixed-bulk/직접 커널 진단용으로 남긴다.
- 직접 커널과 전방 절단 구현은 독립적으로 검증되었다. 성분 폐합에 종을 넣으려면 별도의 과학 모델 결정이
  필요하다.

## 직접 입자 커널, 절단, 분광 범위

생산 입자 위상 경로는 직접 theta-linear 벡터 커널이다. 옛 모먼트 표현은 통제된 진단용으로만 남긴다.

- 기본 / 직접: `OCRT_WATER_PARTICLE_KERNEL=direct`
- 진단용 옛 방식 비교: `OCRT_WATER_PARTICLE_KERNEL=moment`
- 직접 커널을 처음 만들 때의 방위각 구적: `OCRT_WATER_VALUE_NPHI=N`

완전한 푸리에 커널이 한 번 만들어지면 worker 전용 메모리에 캐시되어 케이스 안에서 재사용된다. 요청된
출력 RAA 격자가 SOS 재해석을 유발하지 않는다.

활성 AHN 광물 파일은 공통 위상 파장의 과거 0.5도 절점을 모두 유지한 채 FR631 로 재배치했다. 330 nm 의
임시 광물 위상 계약은 350 nm 끝값 유지다. 공급된 패키지에 표준 AHN 미세물리 생성기가 없었기 때문이다.
원자료 기반 재생성 후보는 검증 폴더에 보존했고 설치하지 않았다.

활성 옛 쇄설물 파일은 원래의 가변 각도 격자에 그대로 있으며 같은 theta-linear 직접 평가기로 계산한다.
기각된 330–1100 nm 쇄설물 후보는 생산에서 제외한다.

정확히 VZA = 0 에서 편광 직접 커널 출력 추출에는 문서화된 특이점이 남는다. 명목 천저의 편광 진단 출력에는
VZA = 0.001 도를 쓴다. 기하를 조용히 바꾸지 않는다.

OCRT 성분 분광 계약은 330–1100 nm 이며, 선택한 각 성분의 명시적 자료 계약을 따른다. 순수해수 광학
특성은 더 넓은 Z09 구간에 표로 있고, 온도가 20 °C 와 다르면 별도의 온도 보정표도 요청 파장을 덮어야 한다.
