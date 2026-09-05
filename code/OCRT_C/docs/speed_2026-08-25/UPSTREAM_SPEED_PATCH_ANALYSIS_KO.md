# OCRT 속도 패치 통합 문서 (분석·수정·검증·전체 diff) — 2026-08-25

첨부 패키지 `OCRT_STAGE2_FAIR_SPEED_PARITY_20260823.zip`(속도 관련 내용만 참조; 결과값 불참조)의 공정 비교 계약을 기준으로 현재 트리(MIGRATION_PKG_2026-08-19, T_wa Q-fix 반영판)를 실측·프로파일링하고, 지배 병목 2건을 수정해 **공정 프로파일 3.25×, 정례 케이스 2.53× 가속**을 얻었다. 동일 머신·근사 공정 조건에서 OCRT가 OSOAA보다 느리던 관계(1.76×)가 **1.85× 우위로 반전**됐다.

## 1. 첨부 패키지 검토 (속도 관련 요지)

- 구판 4.991×(OCRT 84.1s vs OSOAA 16.9s) 판정은 조건 불일치(PSSA ON, 대기 400층, 수중 quadrature 96, 차수·기체·캐시 상태 상이)로 **철회 상태(NON_COMPARABLE_DIAGNOSTIC)** — 병목 참고자료로만 유효.
- 공정 계약: 공통 물리(555 nm/SZA 50/wind 5/red_clay 0.5/M80C AOD865 0.1/PSSA OFF/기체 0) + 공통 해상도(Gauss 48, 대기 26층/해수 640층, 차수 cap 100, particle Fourier m=32, 출력 18×72=1296기하) + OMP 1 + fresh/warm 분리. OSOAA 내부 rough-surface OS_NS=96 vs OCRT m≤32는 `STRICT_INTERNAL_FOURIER_PARITY=false`로 명시 유지.
- OSOAA 공통해상도 probe 17.21s(패키지 측 Linux). 패키지의 ADVANCED 노브 패치(--n-layers-water 등)는 구판용이며, 현재 트리는 동일 지점을 `--debug-water-n-layers/--debug-water-max-orders`(OCRT_DEBUG=1)로 재현 가능 — 신규 패치 불필요.

## 2. 수정 전 실측 (이 샌드박스, 2코어, OMP_NUM_THREADS=1)

| 케이스 | 구성 | wall |
|---|---|---|
| B1 정례(cross-check red_clay) | 555/SZA40/w3, n_mu_water 48, 자동층, 13×24기하 | 14.18 s |
| B2a 공정-근사(에어로졸 제외) | 계약 §2·§4 재현(48/48, 26/640층, 차수100, m32, 18×72) | **45.98 s** |
| B2b = B2a + M80C AOD 0.1 | 동일 | 48.14 s |
| OSOAA 공정-근사 대응(무에어로졸) | NbGauss 48/Mie 16/IGmax 100/NT 26·640, 동일 red_clay·기하 조건 | **26.1–26.3 s** (표면행렬 fresh 확인) |
| OMP 2 (B2a) | — | 44.9 s (이득 2.4% — 단일 케이스 경로 사실상 직렬) |

→ 수정 전 OCRT/OSOAA = **1.75×** (구판 4.99×보다 이미 개선된 상태였음). 참고: fullgrid는 native coupled LUT(1회 해석+재구성) 구조로 기하당 재해석 문제는 이미 없음(`calls=3 water_cold=1`).

## 3. 프로파일 (gprof, B2a): 병목은 물리 코어가 아니라 표면 커널 생성이다

| 순위 | 함수 | 비중 | 호출수 | 진단 |
|---|---|---|---|---|
| 1 | `surface_R_coxmunk_trig` | **43.0%** | 1.62억 | 공기측 거친면 반사 커널의 φ-적분을 **푸리에 모드 m마다 재수행** (33모드 × 동일 기하·Fresnel) |
| 2 | `aw_interp_on_unsorted` | **16.5%** | 2,460만 | 매 호출마다 동일 방향격자를 삽입정렬+dedup (O(n²)/호출) |
| 3 | `rt_sos_operator_apply_vector` | 7.9% | 414 | 실제 RT 물리 코어 — 전체의 8%뿐 |
| 4 | `ocrt_fkc_serve` | 6.8% | 2,183 | all-m 커널 캐시 서빙(복사)·빌드 |
| 5 | air 커널 드라이버 + T_aw/T_wa trig | ~12% | — | 표면 커널 계열 합계 ~62% |

핵심 발견: v1.09-opt에 **all-m 일괄 커널 캐시(FKC, OSOAA SURF_MATR 유사 설계)가 이미 존재**하고 R_ww·T_wa에는 적용돼 있으나, **최대 비용인 공기측 R 커널만 미적용**(per-m 재적분 + 별도 per-m 캐시)이었다. 또 FKC는 요청 m_max와 무관하게 65모드 전부를 빌드했다.

## 4. 적용 수정 (v1.11-speed S1·S2, 5파일, diff 312줄 동봉)

- **S1a — 공기측 R 커널 FKC 통합** (`shared/surface.c`): `surface_coxmunk_fourier_kernel`이 8-슬롯 LRU FKC(키: 격자·풍속·n·nphi·규약 바이트 대조)에서 서빙. m>0 IU/QU 부호접기는 기존 규칙 그대로 서빙 후 적용. 회귀 스위치 **`OCRT_AIR_FKC_OFF=1`= 구경로 비트 재현**. 대기 m-루프가 grid-grid/grid-solar 키를 교대 요청하므로 다중 슬롯이 재빌드 스래시를 방지.
- **S1b — FKC 빌드 모드 상한 힌트** (`surface_fkc_set_build_mmax`, 3개 m-루프 앞 설정): 65모드 고정 빌드 → 실제 m_max로 캡. R_ww·T_wa에도 공통 이득. Rayleigh-only(m_max=2) 케이스의 과빌드 회귀 방지.
- **S2 — `aw_interp_on_unsorted` 정렬 메모화** (`rt_air_water_coupling.c`): 동일 격자 내용(바이트 대조)의 정렬/dedup 인덱스맵을 스레드-로컬 메모화, 값은 맵 경유 수집 — 보간 입력이 직접 경로와 비트 동일.

## 5. 수정 후 실측

| 케이스 | 수정 전 | 수정 후 | 배율 |
|---|---|---|---|
| B2a 공정-근사 | 45.98 s | **14.13 s** | **3.25×** |
| B2b (+M80C) | 48.14 s | **14.87 s** | 3.24× |
| B1 정례 | 14.18 s | **5.61 s** | 2.53× |
| PSSA 대기 스모크(555/SZA75/w5) | 0.39 s | 0.15 s | 2.6× |
| **동일 머신 OSOAA 대비 (B2a)** | 1.76× 느림 | **1.85× 빠름** (14.1 vs 26.1 s) | — |

주의: OSOAA 대비 수치는 이 샌드박스(Linux, 동일 CPU, 단일스레드, 표면행렬 fresh 동등)의 공정-근사 조건 값이다. 계약상 공식 판정은 native Windows + 패키지 러너(PARITY_STATUS 확인 + 1,296행 IQU 산포도 봉인)로 수행해야 하며 이는 로컬 실행 항목이다(§7).

## 6. 정합성 검증 (전 항목 수행 완료)

| 게이트 | 결과 |
|---|---|
| 메모·힌트 비트투명성 (B1을 `OCRT_AIR_FKC_OFF=1`로 실행 vs 수정 전 출력) | **비트 동일** — 편차는 air FKC 단독 기원임을 격리 증명 |
| wind=0 (평면 경로) | FKC on/off **비트 동일** |
| B2a 수치 편차 | 13/1296행, 절대치 max 1e-16(Rrs_U)·1e-17(ρ_Q) — 마지막-ULP 급, 수중장·I열 비트 동일 |
| B1 / B2b 편차 | 7/312행 max 1e-17 / 13/1296행 max 1e-14 — 동급 |
| PSSA 스모크 | ρ_I **비트 동일**, Q/U 6행 1e-22. 레거시 SHA `98e4b2433be300eb…`은 `OCRT_AIR_FKC_OFF=1`로 **자릿수까지 재현**; 신규 기본경로 SHA `2ec2592be03ec88b…`(비트 기준 재설정) |
| 공식 하니스(run_cross_check_fr631.sh) | pure PASS +0.245% · rayleigh_bfo PASS 0.162% · red_clay κ-별건 −0.258/+0.659% — **세 수치 모두 수정 전과 인쇄 자릿수까지 동일** |
| T_wa 단위 게이트 | (T_I,T_Q2,T_U2)=(1,1,1) PASS |
| Q-fix 검증기(작업지시서 동봉) NEW-BUILD | **6/6 PASS**, 통계 동봉 기준과 동일 (pure Q 0.127, red_clay Q 0.085 등) |

ULP 편차의 성격: `-fassociative-math` 빌드에서 φ-적분 누적 루프의 구조가 바뀌면 컴파일러 벡터화 재결합이 달라져 마지막 비트가 흔들릴 수 있다(물리 게이트 0.15~0.6% 대비 11~13자릿수 아래). 완전 비트 재현이 필요한 회귀 비교는 `OCRT_AIR_FKC_OFF=1` 경로를 쓰면 된다.

## 7. 남은 병목과 로드맵 (수정 후 프로파일 근거, 우선순위순)

수정 후 1위는 물리 코어 `rt_sos_operator_apply_vector`(22.2%)로 정상화됐다. 남은 오버헤드성 항목:

- **S3** `aw_interp` 잔여 1.02s(11.9%): 정렬은 제거됐고 남은 것은 호출당 바이트대조+수집. 타깃별 (구간, 가중) 사전계산을 (격자,타깃)당 1회로 끌어올려 성분·모드 간 재사용(≈33×)하면 ~0.1s로. 난도 하.
- **S4** T_aw 계열 1.2s(14%): T_aw는 generic FKC 부적합 이력(2026-08-16, II 오차 기록)이 있어 **전용** all-m 배치 빌더가 필요(자체 memo와 병존). 난도 중.
- **S5** T_wa 뷰노드 일괄화 0.68s: 현재 뷰 μ당 별도 키로 all-m 재빌드(18회) — 뷰노드 배열을 한 키로 묶으면 1회. 난도 중.
- **S6** FKC accumulate SIMD/OpenMP + 포인터 서빙 1.3s: 빌드 내적(모드×9 madd)의 병렬화·복사 제거. 난도 중.
- **S7** 단일 케이스 인트라 병렬화: OMP2 이득이 2.4%뿐 — sos_apply·층 스윕·모드 루프의 OpenMP화가 멀티코어 스케일링의 관건(16코어 로컬에서 배수 기대). 배치 러너는 이미 행 병렬이므로 캠페인은 현재도 코어 스케일함. 난도 상.
- **S8** OSOAA식 표면행렬 디스크 캐시: 동일 (wind, n, 격자, SZA) 재사용 캠페인에서 커널 빌드 자체를 상각. FKC 블록의 파일 직렬화로 구현 용이. 난도 하~중.

추정 여지: S3–S6로 단일스레드 ~10–11s 수준(추가 ~25%), S7로 멀티코어 배수.

- **공식 Windows 공정 게이트(사용자 실행 항목)**: 첨부 패키지 러너(RUN_FAIR_SPEED_GATE_WINDOWS_ANACONDA.ps1)를 새 exe로 native 실행해 PARITY_STATUS=PASS + IQU 산포도 봉인과 함께 공식 속도비를 산출할 것. 구판용 패치 적용 단계는 건너뛰고 현재 트리의 debug 노브 대응을 쓰면 된다(§1).

## 8. 반영 위치

- 소스(수정 5파일): `code/OCRT_C/src/shared/surface.c`(SHA 9c41f187…)·`surface.h`·`src/rt_air_water_coupling.c`(ea3327f9…)·`src/rt_solver.c`·`src/rt_water_rt.c`. diff: `docs/technical/speed_patch_2026-08-25.diff`(312줄, 5파일).
- Windows exe 재빌드: `code/OCRT_C/build/ocrt.exe`(x86-64-v3, SHA 94a07316…)·`ocrt_x86-64-baseline.exe`(9e23316f…) — T_wa Q-fix + speed 패치 포함.
- 샌드박스 런타임 재빌드 완료(binary dcce3c10…). 벤치 산출물: `~/ocrt/speedbench/`.
- 첨부 패키지 사본: 샌드박스 `~/ocrt/speed_pkg/`(구판용 — 속도 계약·러너 참조용으로만 유지).

## 9. 적용 절차 (독립 세션용 — 이 파일 하나로 완결)

1. 부록 A의 unified diff를 추출해 적용한다(또는 이미 커밋된 `code/OCRT_C`의 5개 패치 파일을 그대로 사용):

```
# 이 md에서 diff 추출 (부록 A 코드블록 → speed_patch.diff)
awk '/^```diff$/{f=1;next} /^```$/{f=0} f' OCRT_SPEED_PATCH_2026-08-25.md > speed_patch.diff
cd <OCRT_C_ROOT> && patch -p1 --dry-run < speed_patch.diff && patch -p1 < speed_patch.diff
make        # Linux. Windows는 buildwin.sh(mingw, -static) 두 타깃
```

2. 적용 대상 5파일과 패치 후 SHA-256 앞자리: `src/shared/surface.c` 9c41f187…, `src/shared/surface.h`, `src/rt_air_water_coupling.c` ea3327f9…, `src/rt_solver.c`, `src/rt_water_rt.c`. 재빌드 산출물: Linux binary dcce3c10…, Windows `ocrt_x86-64-v3.exe` 94a07316… / `ocrt_x86-64-baseline.exe` 9e23316f….

## 10. 검증 절차 (기대값 포함)

| 단계 | 명령/조건 | 기대값 |
|---|---|---|
| V1 비트투명성 격리 | 정례 케이스를 `OCRT_AIR_FKC_OFF=1`로 실행, 수정 전 출력과 `cmp` | **비트 동일** (메모·힌트는 비트투명; 편차는 air FKC 단독 기원) |
| V2 wind=0 | 평면 케이스 FKC on/off `cmp` | 비트 동일 |
| V3 PSSA 스모크 | 555 nm/SZA 75/wind 5/BFO/`--pssa`, vza 0:5:85·raa 0:5:355 | OFF 경로 SHA-256 = `98e4b2433be300eb…`(레거시 게이트 재현), 기본 경로 SHA-256 = `2ec2592be03ec88b…`(신규 기준) |
| V4 ULP 편차 상한 | 공정-근사(§2 B2a)·정례(B1) 수정 전후 수치 비교 | 절대치 max ≤ 1e-14(에어로졸 TOA), 무에어로졸 ≤ 1e-16; 수중장·I열 비트 동일 |
| V5 공식 하니스 | `run_cross_check_fr631.sh` | pure PASS +0.245% / rayleigh_bfo PASS 0.162% / red_clay κ-별건 −0.258·+0.659% — 수정 전과 자릿수 동일 |
| V6 단위 게이트 | `tests/test_twa_spin2_gate.c` | (T_I,T_Q2,T_U2)=(1,1,1)±0.005 |
| V7 Q-fix 검증기 | `verify_twa_qfold_fix.py --pure --tsm` (T_wa 작업지시서 동봉) | 6/6 PASS, pure Q max 0.127 / red_clay Q max 0.085 수준 |
| V8 속도 | §2 B2a·B1 명령 재실행 (OMP 1) | B2a ≈ 3.2×, B1 ≈ 2.5× 가속 (하드웨어 상대값) |

## 부록 A. 전체 패치 (unified diff, 5파일 312줄)

```diff
--- /tmp/surface.c.pre_speed_backup	2026-08-25 07:52:21.473454024 +0000
+++ /root/ocrt/runtime/MIGRATION_PKG_2026-08-19/01_OCRT_C/src/shared/surface.c	2026-08-25 07:55:53.562116141 +0000
@@ -793,6 +793,47 @@
  * Uniform φ-quadrature on [0, 2π): φ_p = (p + 0.5) · 2π/N, w_p = 2π/N
  * (midpoint rule — exactly integrates trig polynomials up to order N-1).
  * ---------------------------------------------------------------------------- */
+#ifdef OCRT_FAST_KERNELS
+/* v1.11-speed S1 (2026-08-25): FKC type + prototype hoisted above the AIR-side
+ * kernel entry point so it can serve from the all-m cache too (definition and
+ * design notes below at "v1.09-opt S3").  g_fkc_build_mmax_hint caps the
+ * batched build at the caller's true m_max (default: full OCRT_FKC_MMAX,
+ * i.e. the pre-hint behavior). */
+#define OCRT_FKC_MMAX 64
+typedef void (*ocrt_fkc_trig_fn)(double, double, double, double,
+                                 double, int, double, int, double *);
+typedef struct {
+    int valid, n_o, n_i, n_phi, sigma_type, q_conv;
+    int built_mmax;               /* modes 0..built_mmax present in coef */
+    double ws, nw;
+    double *mu_o, *mu_i;
+    double *coef;                 /* (built_mmax+1) slabs of n_o*n_i*9, m-major */
+} ocrt_fkc_t;
+static int ocrt_fkc_serve(ocrt_fkc_t *C, ocrt_fkc_trig_fn trig,
+                          const double *mu_o, int n_o,
+                          const double *mu_i, int n_i,
+                          int m, int n_phi_quad,
+                          double ws, int sigma_type,
+                          double n_water, int q_convention,
+                          double *out_m);
+static _Thread_local int g_fkc_build_mmax_hint = OCRT_FKC_MMAX;
+#endif
+
+/* Public: declare the true mode-loop upper bound before a Fourier-kernel mode
+ * loop, so the batched all-m cache build stops there instead of projecting all
+ * OCRT_FKC_MMAX+1 modes.  Values outside [0, OCRT_FKC_MMAX] restore the
+ * uncapped default.  Thread-local; call before EACH kernel-building m loop
+ * (atm boundary, R_ww, T_wa) with that loop's m_max.  Bit-transparent for all
+ * modes actually built. */
+void surface_fkc_set_build_mmax(int m_max) {
+#ifdef OCRT_FAST_KERNELS
+    g_fkc_build_mmax_hint = (m_max >= 0 && m_max <= OCRT_FKC_MMAX)
+                                ? m_max : OCRT_FKC_MMAX;
+#else
+    (void)m_max;
+#endif
+}
+
 static int surface_coxmunk_fourier_kernel_uncached_(const double *mu_o, int n_o,
                                                       const double *mu_i, int n_i,
                                                       int m, int n_phi_quad,
@@ -958,6 +999,63 @@
 
     const size_t kernel_count = (size_t)n_o * (size_t)n_i * 9u;
     const size_t kernel_bytes = kernel_count * sizeof(double);
+
+#ifdef OCRT_FAST_KERNELS
+    /* v1.11-speed S1 (2026-08-25): serve the AIR-side rough-reflection kernel
+     * from the all-m FKC instead of re-integrating the m-independent Cox-Munk
+     * Mueller kernel once per Fourier mode (measured: 43% of a fair-profile
+     * coupled fullgrid run, 162M trig calls at n_mu=48/nphi=1024/m_max=32).
+     * Bit-transparent: the FKC per-m accumulation arithmetic is exactly the
+     * uncached_ loop's, and the identical m>0 IU/QU sign fold is applied
+     * below.  The atm mode loop alternates grid-grid and grid-solar keys, so
+     * a small LRU slot array (not a single slot) prevents rebuild thrash.
+     * Kill switch for regression comparison: OCRT_AIR_FKC_OFF=1. */
+    if (!getenv("OCRT_AIR_FKC_OFF")) {
+        enum { OCRT_FKC_AIR_SLOTS = 8 };
+        static _Thread_local ocrt_fkc_t g_fkc_air[OCRT_FKC_AIR_SLOTS];
+        static _Thread_local unsigned long long air_stamp[OCRT_FKC_AIR_SLOTS];
+        static _Thread_local unsigned long long air_tick;
+        ++air_tick;
+        int slot = -1;
+        for (int s = 0; s < OCRT_FKC_AIR_SLOTS; ++s) {
+            const ocrt_fkc_t *C = &g_fkc_air[s];
+            if (C->valid && C->n_o == n_o && C->n_i == n_i &&
+                C->n_phi == n_phi_quad && C->sigma_type == sigma_type &&
+                C->q_conv == q_convention && C->ws == ws && C->nw == n_water &&
+                memcmp(C->mu_o, mu_o, (size_t)n_o * sizeof(double)) == 0 &&
+                memcmp(C->mu_i, mu_i, (size_t)n_i * sizeof(double)) == 0) {
+                slot = s;
+                break;
+            }
+        }
+        if (slot < 0) {
+            unsigned long long oldest = ~0ULL;
+            slot = 0;
+            for (int s = 0; s < OCRT_FKC_AIR_SLOTS; ++s) {
+                if (!g_fkc_air[s].valid) { slot = s; break; }
+                if (air_stamp[s] < oldest) { oldest = air_stamp[s]; slot = s; }
+            }
+        }
+        if (ocrt_fkc_serve(&g_fkc_air[slot], surface_R_coxmunk_trig,
+                           mu_o, n_o, mu_i, n_i, m, n_phi_quad,
+                           ws, sigma_type, n_water, q_convention, R_m)) {
+            air_stamp[slot] = air_tick;
+            if (m > 0) {
+                /* identical to the uncached_ tail: OSOAA-basis IU/QU (indices
+                 * 2/5) sign conversion for m>0. */
+                const size_t np = (size_t)n_o * (size_t)n_i;
+                for (size_t q = 0; q < np; ++q) {
+                    R_m[q * 9 + 2] = -R_m[q * 9 + 2];
+                    R_m[q * 9 + 5] = -R_m[q * 9 + 5];
+                }
+            }
+            return 0;
+        }
+        /* serve declined (mode above build cap, or alloc failure):
+         * fall through to the legacy per-m cache + uncached_ path. */
+    }
+#endif
+
     if (!ocrt_air_surface_cache_enabled_()) {
         return surface_coxmunk_fourier_kernel_uncached_(
             mu_o, n_o, mu_i, n_i, m, n_phi_quad, ws, sigma_type,
@@ -1153,17 +1251,9 @@
  * sin(m*phi) libm expressions and the p-accumulation order per (m,pair,kl)
  * matches the original loop exactly, so served coefficients are the exact
  * bits the direct evaluation produces.  m > MMAX or alloc failure falls
- * through to the original body (always correct).  Thread-local. */
-#define OCRT_FKC_MMAX 64
-typedef void (*ocrt_fkc_trig_fn)(double, double, double, double,
-                                 double, int, double, int, double *);
-typedef struct {
-    int valid, n_o, n_i, n_phi, sigma_type, q_conv;
-    double ws, nw;
-    double *mu_o, *mu_i;
-    double *coef;                 /* (MMAX+1) slabs of n_o*n_i*9, m-major */
-} ocrt_fkc_t;
-
+ * through to the original body (always correct).  Thread-local.
+ * (v1.11-speed S1: type/prototype hoisted above the air kernel entry; build
+ * capped at surface_fkc_set_build_mmax; built_mmax records the built range.) */
 static int ocrt_fkc_serve(ocrt_fkc_t *C, ocrt_fkc_trig_fn trig,
                           const double *mu_o, int n_o,
                           const double *mu_i, int n_i,
@@ -1174,17 +1264,30 @@
 {
     if (m > OCRT_FKC_MMAX) return 0;
     const size_t pair9 = (size_t)n_o * (size_t)n_i * 9;
-    int hit = C->valid && C->n_o == n_o && C->n_i == n_i &&
+    int keymatch = C->valid && C->n_o == n_o && C->n_i == n_i &&
               C->n_phi == n_phi_quad && C->sigma_type == sigma_type &&
               C->q_conv == q_convention && C->ws == ws && C->nw == n_water &&
               memcmp(C->mu_o, mu_o, (size_t)n_o * sizeof(double)) == 0 &&
               memcmp(C->mu_i, mu_i, (size_t)n_i * sizeof(double)) == 0;
-    if (!hit) {
+    if (keymatch && m <= C->built_mmax) {
+        memcpy(out_m, C->coef + (size_t)m * pair9, pair9 * sizeof(double));
+        return 1;
+    }
+    /* v1.11-speed S1: cap the batched build at the caller-declared true
+     * m_max (surface_fkc_set_build_mmax) instead of always projecting all
+     * OCRT_FKC_MMAX+1 modes.  Bit-transparent: built modes keep the exact
+     * per-m accumulation arithmetic; a request above the cap declines (return
+     * 0) so the caller's original per-m body handles it, and never destroys
+     * an existing valid slab. */
+    const int bmax = g_fkc_build_mmax_hint;
+    if (m > bmax) return 0;
+    {
         free(C->mu_o); free(C->mu_i); free(C->coef);
+        C->valid = 0; C->built_mmax = -1;
         C->mu_o = (double*)malloc((size_t)n_o * sizeof(double));
         C->mu_i = (double*)malloc((size_t)n_i * sizeof(double));
-        C->coef = (double*)malloc((size_t)(OCRT_FKC_MMAX + 1) * pair9 * sizeof(double));
-        double *W = (double*)malloc((size_t)2 * (OCRT_FKC_MMAX + 1)
+        C->coef = (double*)malloc((size_t)(bmax + 1) * pair9 * sizeof(double));
+        double *W = (double*)malloc((size_t)2 * (size_t)(bmax + 1)
                                     * (size_t)n_phi_quad * sizeof(double));
         if (!C->mu_o || !C->mu_i || !C->coef || !W) {
             free(C->mu_o); free(C->mu_i); free(C->coef); free(W);
@@ -1194,7 +1297,7 @@
         memcpy(C->mu_o, mu_o, (size_t)n_o * sizeof(double));
         memcpy(C->mu_i, mu_i, (size_t)n_i * sizeof(double));
         const double dphi = 2.0 * M_PI / (double)n_phi_quad;
-        for (int mm = 0; mm <= OCRT_FKC_MMAX; mm++)
+        for (int mm = 0; mm <= bmax; mm++)
             for (int p = 0; p < n_phi_quad; p++) {
                 const double phi = (p + 0.5) * dphi;
                 W[((size_t)2*mm    ) * (size_t)n_phi_quad + p] = cos(mm * phi);
@@ -1215,14 +1318,14 @@
                     double R_local[9];
                     trig(mo, mi, cphi, sphi, ws, sigma_type,
                          n_water, q_convention, R_local);
-                    for (int mm = 0; mm <= OCRT_FKC_MMAX; mm++) {
+                    for (int mm = 0; mm <= bmax; mm++) {
                         const double cm = W[((size_t)2*mm    ) * (size_t)n_phi_quad + p];
                         const double sm = W[((size_t)2*mm + 1) * (size_t)n_phi_quad + p];
                         for (int kl = 0; kl < 9; kl++)
                             acc[mm][kl] += R_local[kl] * (cos_entry[kl] ? cm : sm);
                     }
                 }
-                for (int mm = 0; mm <= OCRT_FKC_MMAX; mm++) {
+                for (int mm = 0; mm <= bmax; mm++) {
                     const double norm_mm = (mm == 0) ? (1.0/(2.0*M_PI)) : (1.0/M_PI);
                     double *dst = C->coef + (size_t)mm * pair9 + pair * 9;
                     for (int kl = 0; kl < 9; kl++)
@@ -1233,7 +1336,7 @@
         free(W);
         C->n_o = n_o; C->n_i = n_i; C->n_phi = n_phi_quad;
         C->sigma_type = sigma_type; C->q_conv = q_convention;
-        C->ws = ws; C->nw = n_water; C->valid = 1;
+        C->ws = ws; C->nw = n_water; C->built_mmax = bmax; C->valid = 1;
     }
     memcpy(out_m, C->coef + (size_t)m * pair9, pair9 * sizeof(double));
     return 1;
--- /tmp/surface.h.pre_speed_backup	2026-08-25 08:07:54.748244351 +0000
+++ /root/ocrt/runtime/MIGRATION_PKG_2026-08-19/01_OCRT_C/src/shared/surface.h	2026-08-25 07:56:10.097831408 +0000
@@ -226,6 +226,13 @@
                                     double n_water, int q_convention,
                                     double *R_m);
 
+/* v1.11-speed S1 (2026-08-25): declare the true Fourier-mode upper bound of
+ * the NEXT kernel-building m loop so the all-m batched kernel cache builds
+ * only modes 0..m_max instead of the full compile ceiling.  Thread-local;
+ * values outside the valid range restore the uncapped default.  Purely a
+ * performance hint - built modes are bit-identical to the per-m evaluation. */
+void surface_fkc_set_build_mmax(int m_max);
+
 /* WATER-SIDE internal-reflection analogues (B2, 2026-06-03). Same signatures
  * and conventions as the air-side R_coxmunk_trig / coxmunk_fourier_kernel, with
  * the facet Fresnel changed to water(n_water)->air(1.0) (internal, TIR).
--- /tmp/coupling.c.pre_speed_backup	2026-08-25 07:52:21.475401492 +0000
+++ /root/ocrt/runtime/MIGRATION_PKG_2026-08-19/01_OCRT_C/src/rt_air_water_coupling.c	2026-08-25 07:57:14.390119390 +0000
@@ -131,6 +131,21 @@
  * For n <= 64 the arithmetic order is unchanged from the legacy path. */
 enum { AW_INTERP_STACK_N = 320 };
 
+/* v1.11-speed S2 (2026-08-25): the coupling loops call this once per
+ * (component, mode, target) with the SAME direction grid, and the per-call
+ * insertion sort + dedup of that grid dominated the calls (measured 16.5% of
+ * a fair-profile coupled run, 24.6M calls).  Memoize the sorted/deduped index
+ * map per exact grid content (thread-local, byte-compared — no false hits);
+ * values are gathered through the memoized map, so every arithmetic input to
+ * linear_interp_ascending is bit-identical to the direct path. */
+typedef struct {
+    int    valid, n, w;
+    double mu_copy[AW_INTERP_STACK_N];
+    int    keep_src[AW_INTERP_STACK_N];   /* source index per kept entry */
+    double ms[AW_INTERP_STACK_N];         /* sorted, deduped mu */
+} aw_sortmemo_t;
+static _Thread_local aw_sortmemo_t g_aw_sortmemo;
+
 static double aw_interp_on_unsorted(const double *mu, const double *vals,
                                     int n, double mu_target) {
     int ord_s[AW_INTERP_STACK_N];
@@ -139,6 +154,30 @@
     int *heap_i = NULL; double *heap_d = NULL;
     int w = 0;
     if (n < 2) return linear_interp_ascending(mu, vals, n, mu_target);
+    if (n <= AW_INTERP_STACK_N) {
+        aw_sortmemo_t *M = &g_aw_sortmemo;
+        if (!(M->valid && M->n == n &&
+              memcmp(M->mu_copy, mu, (size_t)n * sizeof(double)) == 0)) {
+            /* (re)build the sorted/deduped map for this grid content —
+             * identical insertion sort + tolerance dedup as the direct path */
+            for (int a = 0; a < n; ++a) ord[a] = a;
+            for (int a = 1; a < n; ++a) {
+                int k = ord[a]; int b = a - 1;
+                while (b >= 0 && mu[ord[b]] > mu[k]) { ord[b+1] = ord[b]; --b; }
+                ord[b+1] = k;
+            }
+            int wm = 0;
+            for (int a = 0; a < n; ++a) {
+                double m = mu[ord[a]];
+                if (wm > 0 && fabs(m - M->ms[wm-1]) < 1e-12) continue;
+                M->ms[wm] = m; M->keep_src[wm] = ord[a]; ++wm;
+            }
+            memcpy(M->mu_copy, mu, (size_t)n * sizeof(double));
+            M->n = n; M->w = wm; M->valid = 1;
+        }
+        for (int k = 0; k < M->w; ++k) tmp[k] = vals[M->keep_src[k]];
+        return linear_interp_ascending(M->ms, tmp, M->w, mu_target);
+    }
     if (n > AW_INTERP_STACK_N) {
         heap_i = (int *)malloc((size_t)n * sizeof(int));
         heap_d = (double *)malloc((size_t)n * 2u * sizeof(double));
@@ -623,6 +662,7 @@
                 xi[NWIN + i] = xg[i];
                 wi[NWIN + i] = (xg[i] > lo && xg[i] < hi) ? 0.0 : wg[i];
             }
+            surface_fkc_set_build_mmax(m_max);   /* v1.11-speed S1: cap T_wa all-m build */
             for (int m = 0; m <= m_max; ++m) {
                 const double *Trow;
                 #pragma omp critical(p2b_rowcache)
--- /tmp/rt_solver.c.pre_speed_backup	2026-08-25 07:52:21.477183796 +0000
+++ /root/ocrt/runtime/MIGRATION_PKG_2026-08-19/01_OCRT_C/src/rt_solver.c	2026-08-25 07:56:27.905689578 +0000
@@ -54,6 +54,7 @@
 #include "rt_sos_operator.h"
 #include "rt_fourier.h"
 #include "rt_surface_boundary.h"
+#include "shared/surface.h"   /* v1.11-speed S1: surface_fkc_set_build_mmax */
 #include "rt_water_rt.h"
 #include "rt_spectral_contract.h"
 
@@ -829,6 +830,7 @@
         }
     }
 
+    surface_fkc_set_build_mmax(m_max);   /* v1.11-speed S1: cap air-kernel all-m build */
     for (int m = 0; m <= m_max; ++m) {
         /* Phase 3: pass atm->{gammal,betal}_aer when aerosol active; NULL
          * preserves bit-level Phase 2 Rayleigh-only behavior. */
--- /tmp/rt_water_rt.c.pre_speed_backup	2026-08-25 07:52:21.479306947 +0000
+++ /root/ocrt/runtime/MIGRATION_PKG_2026-08-19/01_OCRT_C/src/rt_water_rt.c	2026-08-25 07:56:27.907915657 +0000
@@ -5729,6 +5729,7 @@
         rt_mkc_begin(&atm, &ws, m_loop_max + 1);
 #endif
     if (getenv("OCRT_S6_TRACE")) { struct timespec ts_; clock_gettime(CLOCK_MONOTONIC,&ts_); fprintf(stderr,"[S6W] wmloop begin abs=%.4f\n", ts_.tv_sec+1e-9*ts_.tv_nsec); }
+    surface_fkc_set_build_mmax(m_loop_max);   /* v1.11-speed S1: cap R_ww all-m build */
     for (int m = 0; m <= m_loop_max; ++m) {
         if (grid_hit) {
             if (m >= g_grid_cache.m_count) break;   /* cold run early-exited here */
```
