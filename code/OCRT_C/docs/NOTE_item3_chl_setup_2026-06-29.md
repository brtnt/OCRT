# item#3 chl(phytoplankton) fixed-bulk: 셋업 완료 + rrs −34.8% gap 진단 확정 (2026-06-29)

## 약어
chl = chlorophyll/phytoplankton 농도[mg/m³], fixed-bulk = OCRT 검증모드(외부 a/b/bb/phase 주입), PM_PHYTO = OSOAA phyto 위상 Legendre 모멘트 출력, a/b/bb = 흡수/산란/후방산란[m⁻¹], A = OSOAA truncation coefficient, b* = 절단 산란계수, ω* = 절단 single-scattering albedo, Bp = backscatter fraction(후방산란 비율), moment/value kernel = Legendre/angle-space phase 처리.

## 확립·검증된 것

### OSOAA chl 작동 인자셋 (확정)
```
-HYD.Model 1 -HYD.DirMie <mie_cache_dir>
-PHYTO.Chl <X> -PHYTO.ProfilType 1
-PHYTO.JD.slope 4.0 -PHYTO.JD.rmin 0.01 -PHYTO.JD.rmax 200. -PHYTO.JD.MRwa 1.05 -PHYTO.JD.MIwa -0.000 -PHYTO.JD.rate 1.0
-SED.Csed 0.0 -YS.Abs440 0.0 -DET.Abs440 0.0
-OSOAA.ResFile.Adv.Up RESLUM_Adv_UP.txt   (rrs(0-)용)
```
필수 누락 에러: HYD.Model→pure water fallback, HYD.DirMie→ERROR_2610, PHYTO.ProfilType→ERROR_22400, JD.rate→ERROR_2702. 참조 `exe/run_OSOAA_demo.ksh`.

### IOP 추출 (PROFILE_SEA.txt) — 검증됨
c = TAU_EXT(level)/depth(level), b_w = MOL_PC·c, b_phy = PHY_PC·c, a_total = (1−ΣPC)·c.
검증: chl=1@555 b_w=0.00186 = Rayleigh b_w(555) **정확 일치**. PROFILE_SEA는 "sans ajustement à la troncature" = **절단 전 original b**.
chl=1@555: c=0.36986, b_phy=0.29744, a_total=0.07056, ω=0.8092.

### truncation convention — 일치 확인 (CRUX 해결)
OSOAA(OSOAA_HYDROSOLS.F:1172) `PIZTR = PIZ·(1−A/2)/(1−PIZ·A/2)`, 즉 **b* = b·(1−A/2)** = OCRT(rt_water_rt.c:2038) 동일.
검산: 절단 g' = (g−f)/(1−f), f=A/2=0.7176, g=0.965 → g'=0.876 = PM_PHYTO 모멘트 a₁/3=0.875 **일치**.
→ original b 주입 → OCRT가 OSOAA와 동일 절단. **IOP/ω* parity 확인: OCRT ω*=0.545 = OSOAA ω*(0.809→0.545)**.

## ★ rrs −34.8% gap 진단 확정 (OCRT 2.616e-3 vs OSOAA 4.012e-3)

IOP·ω*·truncation 모두 일치하는데 rrs −34.8%. 원인 격리:
- **bb 무감각**(주입 0.001/0.003/0.006 모두 rrs 동일) → OCRT는 주입 bb 안 쓰고 phase에서 backscatter 유도
- **lmax 무감각**(L 80→200: 2.616→2.582e-3) → 모멘트 개수 아님

**결정적 판정 — 모멘트 closed-form backscatter:**
PM_PHYTO BETA11=a_k(절단 Legendre 계수). Bp = 0.5·Σ a_k·J_k (J_k=∫_{-1}^0 P_k dμ, closed-form, Gibbs 없음) = **0.02521**.
- OSOAA rrs(4.012e-3) 역산 Bp = 0.0246 → **모멘트 진짜값(0.0252)과 일치 ✓ → OSOAA 정확**
- OCRT rrs(2.616e-3) 역산 Bp = 0.0122 → **진짜값의 절반 ✗ → OCRT 틀림**
- 모멘트로 예측 rrs(Gordon) = 4.07e-3 = OSOAA, OCRT 아님

**버그 정체: OCRT의 OSOAA-PM phase 경로(rt_water_rt.c:1765)가 backscatter를 2배 과소평가.**
메커니즘: PM_PHYTO 모멘트는 **이미 절단된** phase인데, OCRT가 각도 재구성 → **OSOAA forward cap 또 적용 → 재정규화**(double-cap) → backscatter 절반.
- value 커널(fixed_bulk_direct_phase_fourier, line 1212)은 **CCRR FF analytic phase 전용**(angular VALUES 소비). OSOAA-PM 모멘트 경로는 이걸 안 씀.
- water_phase_kernel은 "no public CLI override"(line 700) → CLI로 경로 전환 불가.
- memory가 flag한 "moment kernel vs value kernel / Gibbs" 이슈와 정확히 일치.

## fix 경로 (code 수정 필요, focused 작업)
1. **권장**: OSOAA-PM 모멘트는 이미 절단됨 → 재구성·재cap 없이 그대로 moment 커널에 쓰되, 저해상도 backscatter Gibbs를 OS_NB(=m_max·n_mu) 충분히 키워 억제. 또는
2. value 커널 경로로 보내되 **untruncated** angular phase 필요(OSOAA cap 1회 적용). untruncated Mie phase는 miepython(미설치) 또는 OSOAA 텍스트 출력 부재 → Mie 직접 계산 필요.
3. double-cap 제거: rt_water_rt.c:1765 경로에서 OSOAA-PM 입력에 대해 cap·재정규화 skip(이미 절단이므로).

## 성능 이슈
OCRT fixed-bulk가 고산란(chl=1, b~0.3)서 n-mu-water 16+에서 timeout. chl{0.03,0.3,3,30}×6band 위해 n-mu-water-sky + 튜닝 필요.

## 다음 단계
1. **OSOAA-PM phase backscatter 버그 fix**(위 경로1 or 3) → OCRT ω*·rrs가 OSOAA와 일치하는지 재검증
2. 성능 튜닝
3. chl{0.03,0.3,3,30} TOA IQU + Rrs 비교 (고흡수 잔차도 red band 동반 예상)
