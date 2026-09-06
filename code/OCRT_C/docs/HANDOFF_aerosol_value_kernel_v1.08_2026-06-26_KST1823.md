# OCRT v1.08 — aerosol angle-space value-kernel 작업 handoff / 컨텍스트 백업
버전: v1.08_2026-06-26_KST1823 (KST 18:23)

## 한 줄 요약
atmospheric aerosol phase의 backscatter glory 결함(moment kernel의 Legendre Gibbs)을 고치기 위한 **angle-space value kernel을 코드에 완전히 wiring**했고, **빌드·발화까지 확인**했으나 **결과가 틀렸다**(mid/high VZA 과대예측). 다음 작업 = value kernel이 채우는 `phase_fourier_m`의 **정규화/규약을 atmospheric SOS의 moment kernel(`rt_kernel_phase_fourier`)에 맞추는 것**.

## 안전성 (중요)
`--aer-phase-kernel` 기본값 = `moment` (= 기존 동작). value path만 opt-in으로 깨져 있음. 즉 **기존 LUT driver/검증은 영향 없음**. value path만 다음 세션에서 수정 대상.

---

## 1. 목표 & 확정된 root cause (memory #19, 이번 세션 이전 확정)
- 목표: aerosol-only nadir이 OSOAA로 매끄럽게 수렴하도록 aerosol phase를 angle-space로 주입.
- 확정 결함: v1.08 atmospheric aerosol은 **moment kernel 전용**(rt_aerosol.c가 betal/gammal/alphal Legendre moment 계산, rt_solver.c가 moment 소비). raw .mie phase를 Legendre로 재구성하면 backscatter glory(M80C@412: P(140°)=0.155→P(160°)=0.245→P(180°)=0.492)를 못 잡고 Gibbs 진동(L40 P(140°)=−0.096 음수, L80 P(180°)=1.155 과대, L200 P(180°)=0.261 과소; 비수렴).
- 증상: **aerosol-only @412/aod0.2/sza30/raa90 (`--pressure 0`) nadir −20%** (M80C: OCRT 2.0033e-2 vs OSOAA 2.5061e-2). mid-VZA(20~80°)는 moment에서 −0.3~2.5%로 정상.
- 결론(확정): OCRT(moment kernel) 결함, OSOAA(angle-grid raw phase 직접 사용) 옳음.
- 해법(확정): in-water의 value kernel(`fixed_bulk_direct_phase_fourier`)과 동일 원리로 raw P11 **값**을 azimuth Fourier 직접 적분 → moment 우회 → Gibbs 원천 차단.

## 2. 현재 코드 상태 (이번 세션에서 한 일)
value kernel을 **완전히 wiring**:
1. `src/rt_solver.h`: `rt_aerosol_input_t` 구조체에 phase-grid 필드 추가 — `int use_value_kernel; int n_ang_phase; const double *theta_phase,*P11_phase,*P12_phase,*P33_phase;`
2. `src/rt_aerosol_runtime.h`: `rt_aerosol_runtime_options_t`에 `int use_value_kernel;` 추가
3. `src/rt_aerosol_runtime.c`:
   - `rt_aerosol_runtime_prepare`: moment 계산 후 capped phase(P11_w/P12_w/P33_w)를 **free하지 않고** aer에 보존. theta_phase=mie->angles 복사, P11/P12/P33_phase=P11_w/P12_w/P33_w 소유권 이전. `aer->use_value_kernel = ropts->use_value_kernel`.
   - `rt_aerosol_runtime_free`: theta_phase/P11/P12/P33_phase 추가 해제.
4. `src/rt_types.h`: `rt_atm_t`에 borrowed phase-grid 필드 추가 — `int aer_use_value_kernel; int aer_n_ang_phase; const double *aer_theta_phase,*aer_P11_phase,*aer_P12_phase,*aer_P33_phase;` (rt_atm_free가 **free하면 안 됨** — borrowed).
5. `src/rt_solver.c`:
   - includes 직후 `#ifndef M_PI ... #endif` + 두 static 함수: `aer_phase_p11_at_costheta(atm, cth)` (capped P11을 acos(cth)°에서 binary-search 선형보간) + `aerosol_value_phase_fourier(ws, atm, m, nphi)` (kernel; in-water line 1215 구조 그대로: j=0..n_mu, k=−n_mu..n_mu, q=0..nphi, cth=mu_j*mu_k+sj*sk*cos(phi), acc+=P11(cth)*cos(m*phi), phase_fourier_m[j][k]=acc/nphi).
   - **공유 internal 함수 `rt_solve_case_pol_impl`** 안에서: aerosol build(`rt_atm_build_aerosol_rayleigh`) 직후 atm.aer_* 필드를 aer로부터 설정. scalar phase kernel 호출부에서 분기: `if (atm.aerosol_active && atm.aer_use_value_kernel && atm.aer_P11_phase) aerosol_value_phase_fourier(&ws,&atm,m,720); else rt_kernel_phase_fourier(&ws,m,betal_aer);`
   - **주의**: 모든 public wrapper(rt_solve_case_pol / _aerosol / _lut / _for_ocean)가 `rt_solve_case_pol_impl`을 호출하므로 이 한 곳 wiring으로 LUT 포함 전 경로 커버됨.
   - vector kernel(gr/arr/att = P12/P33 moment)은 **아직 moment 그대로** (Q/U 채널 미수정).
6. `src/main.c`: 로컬 `int aer_use_value_kernel = 0;` + `--aer-phase-kernel value|moment` 파싱 + ropts 4곳 중 **init 2(단일-case, line~1498)와 init 4(ropts2, LUT path, line~1674)에 `.use_value_kernel = aer_use_value_kernel` 추가**. init 1(batch struct, mie_path 있음)·init 3(ocean path)은 미수정(의도적).
   - **핵심 발견**: LUT(`--lut --output-full-grid --vector --surface coxmunk`)는 `rt_solve_case_pol_lut`(main.c line~1688) → ropts2(init 4) → aer2를 씀. 처음에 init 2만 고쳐서 발화 안 했음. init 4 고친 뒤 발화 확인.

빌드 OK: `cd /home/user/work/OCRT_AC_LUT_bundle && gcc -std=c11 -O2 -fopenmp -Isrc $(find src -name '*.c') -o build/v2_solver_vk -lm`
(주의: `src/*.c shared/*.c`는 실패 — shared .c는 `src/shared/`에 있음. 반드시 `$(find src -name '*.c')` 사용.)

## 3. 틀린 증상 (정확한 수치)
약어: VZA=뷰 천정각(°), OSOAA=참조 벡터 RT의 REFL 컬럼(=π·L/Ed), moment=기존 Legendre kernel, value=신규 angle-space kernel, d%=100·(OCRT/OSOAA−1).
대상: M80C aerosol-only, 412nm, aod=0.2, sza=30, raa=90, wind=3, `--pressure 0`, `--aer-h-km 2.0 --n-layers 400 --trunc-aer-loglin --aer-l-max 80 --m-max 16`.

| VZA | OSOAA | moment | d% | value | d% |
|----:|------:|-------:|---:|------:|---:|
| 0 | 2.5061e-2 | 2.0033e-2 | −20.1 | 2.1442e-2 | −14.4 |
| 10 | 2.2669e-2 | 1.9681e-2 | −13.2 | 2.1752e-2 | −4.0 |
| 20 | 1.8888e-2 | 1.8157e-2 | −3.9 | 2.2730e-2 | +20.3 |
| 30 | 1.4819e-2 | 1.4780e-2 | −0.3 | 2.4535e-2 | +65.6 |
| 40 | 1.3054e-2 | 1.2956e-2 | −0.8 | 2.7520e-2 | +110.8 |
| 50 | 1.3889e-2 | 1.3740e-2 | −1.1 | 3.2544e-2 | +134.3 |
| 60 | 1.8181e-2 | 1.8015e-2 | −0.9 | 4.1570e-2 | +128.6 |
| 70 | 2.8851e-2 | 2.8203e-2 | −2.2 | 5.9148e-2 | +105.0 |
| 80 | 5.1645e-2 | 5.0331e-2 | −2.5 | 9.6740e-2 | +87.3 |

MAPE: moment 4.99% → value **74.50%** (악화). nadir만 올바른 방향(−20.1→−14.4), 나머지 폭발.

## 4. 진단 진행 상황 & 다음 단계
**핵심 관찰**: value kernel 공식(acc/nphi)은 in-water `fixed_bulk_direct_phase_fourier`(rt_water_rt.c:1215)와 **완전히 동일**. 그런데 결과가 틀림. → **in-water SOS와 atmospheric SOS의 phase_fourier_m 규약이 다르다**. value kernel을 in-water가 아니라 **atmospheric moment kernel(`rt_kernel_phase_fourier`, rt_kernel.c)의 규약에 맞춰야** 한다.

**다음 세션 첫 작업**: `grep -rn rt_kernel_phase_fourier src/rt_kernel.c`로 정의 찾아서 phase_fourier_m을 **정확히 어떻게** 채우는지 확인 (정규화 인자, (2−δ_0m) 포함 여부, P11 정규화 규약).

**PRIMARY 가설 (가장 유력)**: value kernel은 raw `P11_w`(저장된 capped phase)를 그대로 씀. 그러나 moment(β_0)는 β_0=1 규약(=½∫P11 dμ=1). 만약 `rt_aerosol_compute_vector_legendre`가 **내부에서 P11을 재정규화**해서 β_0=1을 만든다면, raw P11_w의 스케일은 β_0=1이 **아닐 수** 있음 → value kernel이 잘못된 스케일 사용. **확인법**: ½∫P11_w dμ를 계산 → 1이 아니면 그 인자로 P11_w를 나눠서 value kernel에 주입.
- 단, 이 가설은 "전 VZA 일정 비율 과대"를 예측하는데 실제는 "VZA 증가할수록 과대"(nadir는 오히려 낮음)라 **순수 스케일 인자만으로는 설명 안 됨**. 스케일 + 구조(인덱싱/대칭/m-mode) 문제 복합 가능성.

**확정 디버깅 방법 (권장)**: moment kernel과 value kernel이 채운 `phase_fourier_m[j][k]`를 **같은 케이스(m=0, 몇 개 j,k)에서 stderr 덤프해 직접 비교**.
- 일정 배수 차이 → 정규화 인자.
- 구조적 차이(특정 j 또는 k에서만, 또는 m별로) → 공식/인덱싱/대칭 문제.
- 후보: (a) P11 정규화(위 PRIMARY), (b) (2−δ_0m) 인자, (c) ½ 인자, (d) j 범위(0..n_mu vs −n_mu..n_mu), (e) atmospheric kernel이 layer별 ssa/tau 가중을 phase_fourier_m에 접지 않는지(in-water와 다를 수 있음).

**VZA 증가=과대 패턴 단서**: 과대가 VZA(grazing)로 갈수록 커짐 → side/backward 산란 기여 또는 고차 m-mode가 과대. moment 대비 value의 m>0 모드가 너무 클 가능성(예: (2−δ_0m) 이중적용 또는 누락). nadir(μ_v=1)은 m=0 지배라 영향 작음 → 일관됨.

## 5. 검증 환경 (matched config, 재유도 금지 — memory #18/#19)
- OCRT 명령(value 테스트): `./build/v2_solver_vk --surface coxmunk --wind-speed 3 --sza 30 --vza 0 --raa 90 --wavelength 412 --pressure 0 --mie inputs/M80C.mie --aod 0.2 --aer-l-max 80 --m-max 16 --aer-h-km 2.0 --n-layers 400 --trunc-aer-loglin --aer-phase-kernel value --vector --lut --lut-vza-step 10 --lut-vza-max 80 --lut-raa-step 30 --output-full-grid /tmp/vk_value.csv`
- moment 비교는 `--aer-phase-kernel moment`로 동일 실행.
- aerosol-only는 `--pressure 0` 필수(Rayleigh OFF). coupled로 보면 412서 Rayleigh가 nadir 지배해 결함이 가려짐(coupled@412 nadir +0.37%로 맞아버림).
- OSOAA 기준: `/home/user/work/OSOAA_OCRT_purewater_default_source` 내 AERONLY M80C 412 sza30 raa90 run의 `LUM_vsVZA.txt`. REFL 컬럼(=π·L/Ed) 사용. OCRT rho_I(=π·L/Ed)와 직접 비교. VZA<=0(upward) 행만, raa=90.
- OSOAA 48-Gauss VZA 격자 → OCRT vza 격자로 linear 보간(nearest 금지).
- OSOAA 재컴파일됨: `inc/OSOAA.h`의 `CTE_MAXNB_ANG_EXT` 200→400.

## 6. tool 전송 깨짐 (이번 세션 반복 문제)
긴 preamble + 복잡한 shell(특히 grep `\|` escape, `${PIPESTATUS}`, 다중 heredoc)이 붙으면 tool 호출이 중간에 잘려 텍스트로 렌더됨. **완화**: tool 호출 전 preamble 짧게, 한 줄 단순 명령, grep alternation은 `\|` 대신 단일 패턴 여러 번 또는 `-E`. create_file은 안정적.

## 7. 이어가기 체크리스트 (다음 세션)
1. 이 파일 + memory #18/#19 먼저 읽기.
2. `grep -rn rt_kernel_phase_fourier src/rt_kernel.c` → moment kernel의 phase_fourier_m 규약 정독.
3. moment vs value phase_fourier_m 덤프 비교(env-gate stderr 덤프 추가, default-off).
4. PRIMARY 가설(P11_w β_0=1 정규화) 먼저 검증: ½∫P11_w dμ 계산.
5. 규약 일치시킨 뒤 aerosol-only nadir이 OSOAA로 수렴 + mid-VZA 무손상 확인(scatter plot, linear axes).
6. 성공 후: vector kernel(gr=P12, arr/att=P33)도 value로 확장 → `--aer-phase-kernel` 기본을 value로 전환 검토 → LUT driver `ocrt_ac_lut.py`에 `--aer-h-km 2.0 --n-layers 400 --trunc-aer-loglin` + value 반영.
7. aerosol phase consistency 하네스 구축(raw vs OCRT-effective backscatter + aerosol-only nadir regression guard).
8. Windows로 migration(static .exe cross-compile) 후 Jae가 올바른 LUT 재생성(현재/과거 LUT는 an23+40layer 잘못된 기본값이라 폐기 대상).
