# 준비 구간 배치화(prep-batching) 완료 노트 — 2026-07-20

## 목표(Jae 지적 반영)
동기식 단일 스트림에서는 케이스별 CPU 처리가 GPU 처리 중간에 하나라도 있으면
그동안 GPU가 대기한다(오버헤드 크기 무관). 따라서 구조를 다음으로 재편했다.
  [host 셋업 전부(스케줄링): 링·스칼라·깊이정책·Greek·에어로졸 Greek]
   → [GPU 배치 전부: 대기·수중 매질 빌드, 표면 커널, 결합, SOS, 조립 — 중간에 host 없음]
   → [결과 readback·인덱스 정렬]
링 구성 같은 순수 host 작업은 GPU가 없는 셋업 단계라 유휴를 만들지 않는다.

## 이번 구간 완성 항목 (전부 CPU 비트 일치)
자료생산 3종(R1=ρ_TOA, R2=ρ_R, R3=ρ_R+A) 전부에서 per-case GPU 커널 런치를 제거했다.

1. 대기 매질 빌드 배치: atmos_batch.build_atm_aerosol_batch + apply_gas_absorption_batch.
   레일리 매질 + 에어로졸 수직프로파일 bisection((B,nt) 벡터화) + 기체흡수(층 그리드
   평탄화/reshape)를 1회 배치 GPU 연산으로. 필드별 전부 최대 절대차 0.00(z_km/xdel/ydel/gas).
   핵심: h 계산 순서를 per-case와 정확히 맞춰야 함((k·tau)/nt). 1 ULP가 얇은 대기 bisection
   민감도로 z_km에서 1e-14로 증폭됨.
2. 레일리 매질 빌드 배치(R2용): atmos_batch.build_atm_rayleigh_batch.
   US62 직접역산(us62_altitude_from_grid_fraction)으로 per-case build_atm_rayleigh와
   비트 일치(에어로졸 bisection 경로와 1e-14 어긋남 회피). 필드별 최대 절대차 0.00.
3. 대기 표면 커널 배치: surface_batch.fourier_kernel_coxmunk_batch(air-side Cox-Munk 반사
   Fourier 커널, B축 + 케이스별 풍속). _coxmunk_surf_kernels의 케이스 루프를 배치 2회
   호출(R_m·R_solar)로 대체하고 surf_seed도 벡터화. 전 m R_m·R_solar 최대 절대차 0.00.
4. forward 결합 T_aw 배치: surface_batch.fourier_kernel_T_aw_batch(air→water 투과 Fourier
   커널, B축). couple_atm_to_water_batch의 케이스 루프를 배치 커널로 대체. 전 m 0.00.
5. R3 값커널 배치: aerosol.aerosol_value_pfm_batch + _aer_p11_at_costheta_batch.
   모든 OPAC 에어로졸이 동일 각도격자(361각도)를 써서 th_asc 공통, P11만 케이스별.
   searchsorted를 공통격자에 nD 질의 + take_along_axis 케이스별 gather. 전 m 0.00.
   solve_atm_batch 값커널 경로의 케이스 루프를 배치판으로 대체(P11 스택 1회).
6. solve_atm_batch 인터페이스 확장: atms 인자가 리스트(내부 스택) 또는 사전빌드 BatchAtm
   (직접 사용)을 모두 받음. R1/R2/R3 전부 배치 BatchAtm 경로 사용.
7. solve_r2_grid/solve_r3_grid 재구성: prep 루프는 링·스칼라·Greek만(host), 그룹별로
   매질 배치 빌드 + 기체흡수 배치 후 solve_atm_batch(bAa)에 전달.

## 검증 (프로덕션 복사본 /tmp/pyocrt_prod_2026-07-20)
- R1 앵커(wl555 sza30vza30raa90 wind3 C50 chl1.2): rho_TOA 2.2579111759e-01,
  Rrs 6.9346186098e-02 정확 일치.
- R1 B=2(다른 에어로졸·각도·조성, wl490): 최대 상대차 2.94e-16(1 ULP). 대기 매질·표면
  커널·forward 결합 배치화 각 단계마다 값 불변.
- R2 배치 vs per-case 리스트: 상대차 0.00(2케이스). rho_R=6.2837331469e-02.
- R3 앵커(wl490 C50): rho_R+A 8.0602185738e-02(값커널·매질 배치화 후 불변, SOS 차수
  잡음 1.78e-12은 게이트 무관).

## 구조 상태
prep 루프는 순수 host 스케줄링(GPU 매질 빌드 없음). GPU 단계(R1 _r1_group, R2/R3 그룹
루프)에 per-case GPU 커널 런치 없음. 남은 케이스별 Python 루프는 전부 host 스케줄링
(리스트 스택, 링 구성, 스칼라 계산, (B,) 배열 조립)이다.
- atmos_batch 832행 값커널: 이제 배치판(케이스 루프 없음).
- atmos_batch 1277·1285행: host 리스트→(B,) 배열 변환(GPU 런치 아님).
- 층 루프(integrate_up/dn)·링노드 루프(diffuse-top): 순차 불가피(케이스 아님).

## GPU 게이트 (Jae 로컬)
샌드박스 GPU 없음 → numpy 비트 기준까지 확인. Jae 로컬 CUDA(RTX 5090)에서 gpu_gate.py로
numpy vs cupy 7유효자리 게이트 통과 후 1200케이스(100 × 12밴드) 자료생산에 사용.
실행: run_produce_gpu.bat <grid> <out> (set OCRT_PY_GPU=1).

## 남은 것
- GPU 게이트 실측(Jae 로컬 CUDA): gpu_gate.py numpy vs cupy 7유효자리.
- CCRR 구성모델(Jae 로컬 준비), wind=0 pass-2, PSSA(기본 OFF), V2/V3 12자리 provenance.

---

# GPU 유휴 진단·층 루프 스캔 벡터화 — 2026-07-20 (추가)

## 증상
Jae 로컬(RTX 5090, WDDM)에서 gpu_gate.py 실행 시 GPU-Util 12%, 71W로 사실상 논다.
Ctrl+C 트레이스백이 integrate_bcs_batch의 층 루프(for k in range(nt))에서 멈춤.

## 원인
- 멈춘 곳은 케이스 루프가 아니라 **층 루프**다. 케이스별 GPU 런치는 이전 구간에 없앴지만,
  SOS는 (m모드 × 차수 × 층 × 패스) 만큼 작은 커널을 순차로 던진다. 이 런치 "횟수"는 B와
  무관하다. WDDM은 커널 실행 지연이 커서(수십 µs~ms), 수만~수백만 개의 작은 순차 런치
  사이에 GPU가 논다. B=6이라 각 연산이 작아 더 심하다.
- 게이트가 전체 자료생산 파라미터(nt_atm=400, m_max=16, max_it_water=500)로 돌아 층·차수·m
  조합이 특히 많다.

## 조치: 층 루프 스캔 벡터화 (integrate_bcs_batch)
층 점화식 I[k]=c[k]·I[k±1]+contrib[k]를 나눗셈 없는 **아핀사상 Hillis-Steele 병렬 스캔**으로
바꿨다. 아핀 합성 (a,b)#(a',b')=(a·a', a'·b+b'). a는 투과율(≤1) 곱이라 [0,1], b는 유계
복사량 누적이라 언더플로·오버플로 없음(큰 광학두께 수중에도 안전). nt개 순차 런치가
log2(nt)개 전배열 단계로 감소(약 10배 런치 감소).
- 검증(스캔 vs 순차, CPU): 대기 τ~0.5(nt100/400) 1.1e-11, 수중 τ~30(nt600) 1.65e-11,
  수중 τ~100(nt1300) 1.89e-11. NaN 없음. 전부 7유효자리 이내.
- R1 전체: 스캔 vs 순차 rho 1.23e-16(사실상 비트 일치), 앵커 대비 1.05e-11.
- R2 5.68e-12, R3 1.78e-12(앵커 대비). 전부 7유효자리.
- 검증용 비트 경로 보존: 환경변수 OCRT_PY_SEQ_BCS=1이면 순차 층 루프(비트 일치) 사용.
- 주의: 스캔은 부동소수 재결합이라 순차 folding과 비트 일치가 아님(~1e-11). 게이트가 이미
  7유효자리 기준이라 일관됨.

## 남은 유휴 원인(구조적)
- 차수 루프(Neumann 급수, 717·1004행)는 각 차수가 앞 차수에 의존해 벡터화 불가. 다만 이제
  각 차수의 층 연산이 스캔이라 차수당 런치가 nt→log2(nt)로 줄었다. 수렴 판정은 이미
  8차수마다만 GPU→host 동기화(수렴 케이스는 0 기여 마스킹)라 per-차수 동기화 없음.
- **B가 작으면(예: 6) 각 배치 연산이 작아 WDDM 런치 지연에 계속 묶인다.** GPU를 채우려면
  자료생산에서 chunk(B)를 크게 잡아야 각 연산이 커진다.

## 권장 사용법
- 정확성 게이트(빠른 확인): 축소 파라미터로 실행.
  python gpu_gate.py --grid <소규모>.csv --data data --bands 555 --chunk 6 \
      --nt-atm 60 --m-max 4 --max-it-water 150
  (numpy 대 cupy 7유효자리 대조가 목적. 전체 해상도 불필요.)
- 처리량(GPU 채우기): chunk를 크게(메모리 한도 내 200~500 권장). 각 층/스캔 연산이
  (B, nt, n_mu)로 커져 GPU가 바빠진다. 24GB에서 수중 nt~1300이면 B~200-500이 현실적.
- WDDM에서 nvidia-smi는 프로세스를 안 보여주므로(정상), Power/Memory/Util로 판단한다.

---

# 속도 개선 1차: 스토크스 묶음 + 프로파일 진단 — 2026-07-20 (추가)

## 게이트 결과(Jae 로컬, 축소 파라미터 nt-atm60 m-max4 max-it-water150)
6케이스 1밴드: GPU 35.8s, CPU 418.4s, 가속 11.7배. 정확성 7유효자리 전 항목 통과
(rho_TOA 2.49e-14, rho_R 1.60e-14, rho_RpA 5.66e-14, Rrs 2.76e-15).

## 프로파일 진단(B=6 CPU, 런치 횟수 = GPU 런치 수 대응)
- _bcs_scan_fwd 2658회, mixdot 15588회, integrate_bcs 1329회가 최다 런치.
- 표면 커널(_mat3_mul 등)은 20~50회로 적고 각 연산이 커서 GPU에선 빠름(문제 아님).
- 병목은 "작은 런치가 수천~수만 번"이고, B=6이라 각 커널이 작아 WDDM 런치 지연에 묶임.

## 조치: 스토크스(I·Q·U) 묶음
integrate_bcs가 I·Q·U를 따로 3번 불렸다. 층 투과율(c)이 세 성분에 공통이므로, 스토크스
축(S)으로 묶어 한 번에 처리하도록 integrate_bcs_batch·_bcs_scan_fwd를 S축 지원으로
확장하고 integrate_bcs_iqu 헬퍼를 추가했다. 성분 간 독립이라 비트 일치(재결합 없음).
- 효과: _bcs_scan_fwd 2658→886회, integrate_bcs 1329→443회(각 3배 감소). GPU 런치 수
  3배 감소(CPU 총시간은 연산량 동일이라 비슷).
- 검증: R1 앵커 대비 rho 1.05e-11·Rrs 1.00e-12, R2 5.68e-12, R3 1.78e-12. 전부 7유효자리.

## 핵심 권장: 큰 B로 측정
게이트는 B=6이라 각 커널이 작아 GPU가 여전히 런치에 묶인다. **자료생산에서 chunk(B)를
크게 잡으면(200~500) 같은 런치 수라도 커널이 33~80배 커져 GPU가 찬다.** 11.7배는 B=6의
런치 묶임 상태 수치이며, 큰 B에서 훨씬 커질 것으로 예상된다. 다음 최적화의 대상(런치가
여전한지, 메모리 대역인지)은 큰 B 측정 결과에 따라 정해야 하므로, 큰 B 측정을 먼저 권한다.

## 다음 최적화 후보(측정 후 결정)
1. mixdot(소스 빌드) 런치 감소: 현재 소스 빌드당 36회 행렬곱. F(6종)별로 a[e]를 쌓아
   einsum 1회로 묶으면 런치 약 3배 감소 가능(재구현 위험 있어 측정 후 진행).
2. m모드 묶음: m루프(4~17회)를 배치축에 접어 넣으면 차수 루프의 모든 런치가 m배 감소.
   가장 큰 구조적 이득이나 per-m 커널·수렴·재구성 처리로 복잡. 큰 B로 부족할 때 진행.

---

# 출력 형식 확장: 입력값 + td_s/td_v + rho_wn — 2026-07-20 (추가)

## 추가된 컬럼 (총 21열)
- 입력값(입력 CSV 그대로): sza, vza, raa, wind, aerosol, aod865, chl, tsm, acdom440.
- 대기 투과율: td_s(하향, 태양→표층), td_v(상향, 표층→센서).
- 정규화 water-leaving 반사도: rho_wn = (rho_RC - rho_A_RA)/(td_s·td_v) = t_rho_w/(td_s·td_v).
- 기존 컬럼(rho_TOA, rho_R, rho_RpA, rho_RC, rho_A_RA, t_rho_w, Rrs)은 유지.

## td_s·td_v 정의 (C rt_solver.c와 동일)
- td_s = T_dir_dn + T_diff_dn_hemi = exp(-tau_atm/mu_sun) + T_diff_dn_hemi.
- td_v = T_total_up_view = exp(-tau_atm/mu_view) + T_diff_dn_hemi (상반성: 상향 확산=하향 확산 반구).
- tau_atm = tau_R + tau_a_eff + tau_abs_col. T_diff_dn_hemi는 대기 1차 패스(res1)에서 산출됨.
- 이 총투과율 정의에서 항등식 t_rho_w = td_s·td_v·rho_w 가 성립하므로, rho_wn 은 표층 위
  water-leaving 반사도 rho_w 의 추정값이다.

## 검증
- td_s·td_v는 0~1 범위(490nm 예: 0.87, 0.90)로 물리적으로 타당.
- rho_wn/(pi·Rrs) = 0.9994(저에어로졸)~0.9375(고에어로졸·고산란). pi·Rrs는 전체 RT로 구한
  표층 rho_w이므로, 두 값이 근접함은 td_s·td_v 정의가 일관됨을 확인한다. 잔차(수%)는 단순식
  투과율이 실제 RT 두방향 투과율의 근사인 데서 오는 대기보정 오차이며, rho_wn(대기보정 추정)과
  Rrs(전체 RT 참값)를 함께 두면 대기보정 정확도 평가에 쓸 수 있다.

## 코드 변경
- batch_driver.solve_r1_grid 반환이 (rho, rrs) -> (rho, rrs, td_s, td_v)로 변경.
  _r1_group 시그니처에 td_s_out, td_v_out 추가, 최종 조립부에서 산출·채움.
- produce_grid.py: OUT_FIELDS 21열로 확장, 행 조립에서 입력값·td_s·td_v·rho_wn 기록.
- gpu_gate.py _run: solve_r1_grid 4-tuple 언팩(정확성 게이트는 rho/rrs만 사용).

## 주의: 재개(resume)
기존 out_100.csv는 옛 9열 형식이다. 새 21열 형식과 한 파일에서 섞이면 열 수가 불일치하므로,
새 형식 자료는 **새 출력 파일**로 생성한다(옛 파일에 append 금지).

---

# td_v 정정: 해석식 -> RT-정확 비율 — 2026-07-20 (추가)

## 문제
상향 투과율 td_v를 처음에 C 해석식(exp(-tau/mu_view)+T_diff_dn_hemi, 상반성 근사)으로
구현했으나, 이는 근사다. 고산란·저태양 케이스에서 이 근사가 8%까지 틀려(td_v 0.955 대
실제 0.883) rho_wn/(pi·Rrs)가 0.9375로 6% 어긋났다.

## 정정 (Jae 지시)
상향 투과율은 별도 루틴 없이 RT 결과의 비율로 얻는다:
  td_v = Lu_TOA / Lu_0plus = (TOA 수중발광 라디언스)/(표층 수중발광 라디언스, 센서방향).
두 양은 R1에서 이미 산출(TOA_wl=res2['I_TOA'], Lu_0plus=resW['Lu_0plus']).
하향은 td_s = Ed(0+)/(F_sun·mu_sun)로 명시화(= 기존 exp(-tau/mu_sun)+T_diff_dn_hemi와 동일).

## 검증
비율법 적용 후 rho_wn/(pi·Rrs) = 0.99998, 0.99993 (해석식법 0.9375 -> 개선). 남은 잔차
0.007%는 차감식(R1 해양 대기경로 대 R3 흑색 coxmunk 대기경로)의 결합 오차이며, 이는 실제
대기보정 분해 오차다(투과율 근사 아티팩트가 아님). rho_wn과 Rrs를 함께 두면 이 차감식
오차를 정량화할 수 있다.

---

# C v1.2 RT-derived transmittance fix와 파이썬 일치화 — 2026-07-21 (추가)

## C fix 요약(첨부 zip)
상향 투과율의 상반성 근사(T_diff_up_hemi=T_diff_dn_hemi, T_total_up_view=exp(-tau/mu_view)
+T_diff_dn_hemi)를 제거하고, 실제 bottom-source pass 결과로 T_total_up_view=(I_TOA_total
-I_atm_path)/Lu0plus로 변경. T_total_dn_hemi(직달+확산), T_up_rt_valid 플래그, 성분 분리
필드를 추가.

## 이미 일치했던 것
- td_v(=TOA_wl/Lu_0plus) = C의 T_total_up_view. 지난 세션에 이미 비율법으로 교체됨.
- td_s(=Ed_0plus/(F_sun·mu_sun)) = C의 T_total_dn_hemi.
- 검증: rho_wn/(pi·Rrs)=0.99998. 정의가 C와 동일.
- C가 언급한 "ocean full-grid driver가 aerosol을 solver에 전달 안 함" 버그: 파이썬에는
  없음. produce_grid의 R1은 aerosol을 대기 패스에 전달하며(앞서 C 해양 실행과 rho_TOA
  5.6e-12 일치), rho_TOA가 에어로졸 종류에 따라 정확히 달라짐.
- C의 flux-ratio 진단 컬럼 개명(rho_glint*mu_sun 등): 파이썬은 해당 컬럼을 출력하지 않아 무관.

## 이번에 반영한 것(C 표준 명명·정의로 CSV 확장)
produce_grid 출력에 C stdout과 동일한 투과율 필드를 추가(td_s/td_v 컬럼을 대체):
- T_dir_dn = exp(-tau_atm/mu_sun)
- T_diff_dn_hemi = 대기 SOS BOA 하향확산 반구적분
- T_total_dn_hemi = T_dir_dn + T_diff_dn_hemi (구 td_s)
- T_dir_up_view = exp(-tau_atm/mu_view)
- T_diff_up_view = T_total_up_view - T_dir_up_view (각도재분배로 음수 가능)
- T_total_up_view = TOA_water_signal_I / Lu0plus (구 td_v, RT 유효 상향투과율)
- TOA_water_signal_I = I_TOA_total - I_atm_path (= TOA 수광 radiance)
- T_up_rt_valid = 1(수광 source 존재; Lu0plus=0이면 상향 NaN·플래그 0)
- rho_wn = t_rho_w/(T_total_dn_hemi·T_total_up_view) 유지
검증: 항등식 T_total_dn_hemi=T_dir_dn+T_diff_dn_hemi(1e-11), T_total_up_view=T_dir_up_view
+T_diff_up_view(0), rho_wn/(pi·Rrs)=0.99998. CSV 27열.

## 미반영(가용성 문제, 필요 시 후속)
- T_diff_dn_dir(방향성 하향확산): solve_atm_batch가 현재 반환하지 않음. 추가하려면 대기 SOS
  BOA 하향장에서 관측방향 성분을 뽑아 반환하도록 solver 수정 필요. 요청 시 반영.

## 주의
CSV 컬럼이 바뀌었으므로(td_s/td_v -> 8개 투과율 필드) 자료는 새 출력 파일로 생성한다.

---

## 2026-07-21 추가: 구성성분 IOP 15개 컬럼 + Kd(0⁻)

### 추가된 CSV 컬럼 (Rrs 뒤, 총 42컬럼)
각 행에 다음 15개를 추가했다. 값은 파장·수중 IOP에만 의존하고 (vza,raa)에는
의존하지 않는다(행 = 한 조건).

- 총량: a_total, b_total, bb_total
- 순수수: a_w, b_w, bb_w
- 식물플랑크톤 흡수 단독: a_chl (= a_phyto)
- 입자 결합(식물플랑크톤+detritus): a_phyto_detritus, b_phyto_detritus, bb_phyto_detritus
- CDOM 흡수: a_dom
- 무기물: a_min, b_min, bb_min
- 수면 직하 감쇠: Kd0minus

### 명칭·정의 (C 대응)
- a_phyto_detritus = a_phyto + a_det = C의 a_pig (rt_water_rt.c:3059,
  OCRT는 comp.pigment=0). b/bb도 동일하게 결합.
- a_chl = a_phyto 단독 = C의 comp.eap_phyto.a. detritus 흡수 0이면(기본)
  a_chl = a_phyto_detritus.
- a_dom = C의 a_cdom_used.
- Kd0minus = C의 Kd_0minus. 정의 Kd = -ln(Ed(z1)/Ed(0⁻))/z1,
  z1 = tau_lvl1/ext, ext=a+b. m=0 하향조도 2준위 로그미분.

### 구현 위치
- ocrt_py/atmos_batch.py: solve_water_batch_r1에 ext_list 파라미터 추가,
  m=0 장(ti_m0) 포착, out['Kd0minus']/out['Ed_0minus'] 산출.
- ocrt_py/batch_driver.py: IOP_KEYS(15개) 정의, _r1_group에 iop_out 인자,
  prep[i]['r']에서 성분 채움, solve_r1_grid가 4-튜플(rho,rrs,trans,iop) 반환.
- produce_grid.py: OUT_FIELDS에 15컬럼 추가, 4-튜플 언팩·행 기록.
- gpu_gate.py: 4-튜플 언팩.

### 검증 (C 단일실행 대조)
조건: --surface ocean --water-model ocrt --ocrt-chl 1.2 --ocrt-tsm 3.5
--ocrt-adom440 0.04 --ocrt-phyto-group micro --ocrt-tsm-species red_clay
--sza 30 --vza 30 --raa 90 --wind-speed 3 --wavelength 490 --pressure 1013.25
--mie inputs/C50.mie --aod-865 0.2 --n-mu-water 16 --n-layers 100 --m-max 3.
- 성분 IOP 12개: C와 8자리 일치(a_w=1.500007e-2, a_pig=3.067200e-2,
  a_min=1.904000e-1, a_total=2.559355e-1 등).
- Kd0minus: 파이썬 -0.85427 = C -0.854268 (~1e-5).
- 맑은 물(chl=0.1 tsm=0.1): 파이썬 -0.0465 = C -0.046498.

### Kd0minus 부호 (하늘광 미산란 투과분 수정 완료)
Kd0minus는 수면 직하 하향확산 감쇠계수로 물리적으로 반드시 양수(~a+bb)다.
초기 이식본은 근표층 하향조도 Ed에서 하늘광(확산-top)의 미산란 투과분을
빠뜨려 음수가 나왔다(C도 같은 버그). 미산란 투과분
Ed_sky(k)=2π·f_scale·Σ_c dtI(μ_c)·exp(-h_k/μ_c)·μ_c·w_c 를 level 0·1 Ed에
더해 Kd를 재계산하도록 수정했다. 결과:
- 탁수(chl1.2): Kd0minus=+0.342 (a+bb=0.325 근접), Ed(0⁻)=0.837=C 0.839 일치.
- 맑음(chl0.1): Kd0minus=+0.028 (a+bb=0.032 근접).
에너지 보존(Gershun: E=Ed−Eu 단조감소)으로 물리성 확증. 수정 위치는
ocrt_py/atmos_batch.py solve_water_batch_r1(ext_w_list 파라미터 추가),
batch_driver.py 호출부(w_wg 전달). C는 같은 버그가 남아 있어(출력 Kd가
장 기반 음수) C_WORKORDER 작업 3-4에서 수정을 지시한다.

### C 작업지시서
같은 항목을 C 해양 전체격자 CSV에도 추가하는 지시서를 첨부한다:
C_WORKORDER_IOP_KD_FULLGRID_CSV.md. 요지: C는 성분 IOP·Kd0minus를 단일실행
stdout에는 이미 찍으나(rt_io.c) 전체격자 CSV(main.c:1068)에는 총량만 있다.
(1) 결과 구조체에 a_chl_used(=a_phyto) 노출, (2) 전체격자 CSV 헤더·행에 성분
IOP·a_chl·Kd0minus 추가가 필요하다. Kd 계산은 이미 있으므로 CSV 출력만 늘린다.

## 주의 (갱신)
CSV 컬럼이 다시 바뀌었으므로(투과율 8개 + IOP·Kd 15개) 자료는 새 출력 파일로
생성한다. 기존 파일과 혼용 금지.
