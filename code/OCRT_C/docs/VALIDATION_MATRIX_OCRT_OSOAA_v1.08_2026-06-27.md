# OCRT ↔ OSOAA 검증 matrix

> ⛔ **세션 시작 시 `COMPARISON_PITFALLS_AND_PROCEDURE_v1.08.md`를 맨 먼저 읽고 §2 체크리스트를 준수할 것.** OSOAA는 sandbox서 빌드·실행 가능(그 문서 §4 recipe). 기억/prose 참조값과 비교 금지 — 그게 −20% 같은 가짜 버그의 원인이었다. (done/todo + 재현 프로토콜)
버전 v1.08 · 2026-06-27 · 하네스 문서(`HARNESS_OCRT_OSOAA_v1.08`) 동반. 규약·명령은 하네스 문서 참조.

> 이 문서는 **무엇을 검증했고(done) 무엇이 남았는지(todo)**를 추적하고, **매 세션 시작 시 done 항목의 byte-단위 재현 테스트**를 정의한다. 모든 항목은 공통 기하 grid를 쓴다.

---

## 0. 공통 기하 grid (모든 항목 공통)

sza ∈ {0, 40, 80}°, vza ∈ {0, 30, 60}°, raa ∈ {0, 45, 90, 135, 180}° → nominal **45**.

**물리적 축퇴 (불필요 run 방지):**
- **sza=0**: 태양 천정 → 방위대칭 → 주어진 vza에서 5개 raa 모두 동일, U=0. → vza당 1개. (sza0 행: 15 nominal → **3 distinct**)
- **vza=0 (nadir)**: 시선 on-axis → 5개 raa 모두 동일. (sza≠0 × vza0: 10 nominal → **2 distinct**)
- **sza≠0 × vza≠0**: 2 sza × 2 vza × 5 raa = **20 distinct**.
- → **구별되는 기하 = 3 + 2 + 20 = 25** (45 중). 검증·재현은 25개 기준.

**raa 규약 주의(하네스 §1.0)**: OCRT glint는 raa=180. **RAA=90만 convention-invariant** — raa{0,45,135,180}은 OSOAA와 비교 전 매핑/산란각 Θ 일치 확인 필수(레벨별 매핑 상이).

**스펙트럼 밴드 (확정)**: {412, 443, 490, 555, 660, 865} nm. aot는 865nm 기준.

---

## 1. 검증 항목 (6개)

약어: IQU=Stokes 반사도 성분, Rrs=표층 remote-sensing reflectance(0−/0+), aot865=865nm aerosol optical depth, aDOM=yellow substance(CDOM) 흡수, TSM=total suspended matter, MOT=molecular optical thickness.

| # | 항목 | 레벨/물리량 | surface | 변수 | 상태 |
|---|---|---|---|---|---|
| 1 | Black Fresnel ocean, Rayleigh / aerosol / Ray+aerosol | TOA IQU 반사도 | black Fresnel | aot865 {0.05, 0.3, 1.0} | **부분 PASS** (아래) |
| 2 | Pure ocean | TOA Rayleigh IQU + 표층 Rrs | Fresnel ocean | — | **PARTIAL** (skylight Ed/Lu 버그 FIX-SKY-EDLU 해결; sza40 OSOAA 일치, TOA 0.63%; 잔차: 412/443 분자위상 +4.3/+1.8%, 660/865 고흡수 −1.5%, sza80 grazing +2%) |
| 3 | 다양한 chl ocean | TOA Ray IQU + Rrs | Fresnel ocean | chl {0.03, 0.3, 3, 30} mg/m³ | **TODO** (native delta-M 버그) |
| 4 | 다양한 aDOM ocean | TOA Ray IQU + Rrs | Fresnel ocean | aDOM(440) {0.01, 0.1, 1} m⁻¹ | **PARTIAL** (sza40: TOA 0.52%, Rrs 흡수구동 잔차 -0.5~-1.9%; sza0/sza80 + 고흡수 잔차 규명 남음) |
| 5 | 다양한 TSM ocean | TOA Ray IQU + Rrs | Fresnel ocean | TSM {0.1, 1, 10, 50} g/m³ | **TODO** |
| 6 | TSM+aDOM+chl, Ray+aerosol | TOA IQU + Rrs | Fresnel ocean | 복합 | **TODO** (최종 통합 — 통과 시 정합성 종료) |

각 항목 = (밴드) × (변수) × (25 기하). 통과 기준은 §3.

---

## 2. 항목별 상세 상태 (done/todo 근거)

**#1 Black Fresnel ocean TOA (Rayleigh / aerosol / Ray+aerosol)** — **OSOAA sandbox 직접 대조(2026-06-28)**
- **Rayleigh ρ_I**: sza≤40 전 밴드 **MAPE ~0.4% PASS**. sza=80 −2~−7%(λ↑에 증가, **grazing surface BRDF**로 추정, registry P7 open).
- **Rayleigh Q/U**: 크기 일치, **U 부호 handedness 반전**(OCRT_U=−OSOAA_U), nadir Q 부호 모호 — 규약 적용 시 일치.
- **aerosol IQU 전 grid**(M80C, 6밴드×aot{0.05,0.3,1.0}×sza, raa90): **sza≤40 I MAPE 0.48%, Q RMS/I 0.28%, U RMS/I 0.20% PASS**(U 부호규약). sza=80은 P7(surface, aot↑에 6.2%→0.5%). **★구 −20%는 폐기** — 틀린 참조값(2.5061e-2)이 원인, 정상 OSOAA(2.0116e-2)는 OCRT(2.0033e-2)와 −0.4%. 이전 "surface/MS 분해/open nadir gap" 진단 전부 무효.
- **ray_aer IQU 전 grid**(raa90): **sza≤40 I MAPE 0.41%, Q RMS/I 0.27%, U RMS/I 0.41% PASS**. sza=80도 0.73%(결합이 P7 희석). 
- **#1 raa90 IQU 3종(Rayleigh/aerosol/ray_aer) 전부 PASS(sza≤40).** 남은 것: raa{0,45,135,180} Φ매핑 확장, P7(sza=80 surface) 규명.

**#2 Pure ocean (TOA IQU + Rrs(0−) 동시검증, 2026-06-28)** — OCRT in-water RT ON, 한 OSOAA run 북킵핑(Adv_UP 전레벨+Flux Ed). water IOP는 OSOAA와 exact match(재구성 LUT).
- **TOA IQU**: sza≤40 I MAPE **0.79%**, Q RMS/I 0.57%, U RMS/I 0.93%(U 부호규약). black water(#1, 0.4%)보다 약간 큼 = water-leaving이 in-water 차이를 운반. sza=80은 P7(grazing, 3.10%).
- **Rrs(0−)**: 단일버그 아님 — **파장 2중구조**. (a) **skylight Ed/Lu 버그**: 저-omega(490+) −2%, grazing(sza80) −5%까지(OCRT 낮음, atm-ON 시 발생). (b) **고-omega 효과**: 412 **+3%**(OCRT 높음) → 분자 위상/depolarization 차이 추정(coupled Rayleigh +2% mid-VZA@412와 동일 근원 가능).
- **남은 것**: skylight 버그 수정(direction B), 412 고-omega 위상차 규명, Rrs off-nadir Snell 정밀화.

**#3 chl ocean** — chl=1 phyto 일부 spot 검증. **native simple-chl delta-M 절단 버그**(412nm rrs +48% vs MC truth, 하네스 §8). ⇒ 버그 수정 후 검증. TODO.

**#4 aDOM ocean** — 체계적 검증 없음. TODO. (chl=0, TSM=0, aDOM만 변화.)

**#5 TSM ocean** — TSM_mineral 일부 작업 이력. 체계적 OSOAA 비교 없음. TODO.

**#6 복합(TSM+aDOM+chl, Ray+aerosol)** — 최종 통합. TODO.

⇒ **결론: full grid 기준 모두 미완.** #2가 부분 anchor, #1-aerosol이 active open. #1 grid 기준 정합성 테스트의 시작점.

---

## 3. 통과 기준 (vs OSOAA)
- **반사도(IQU TOA)**: 글린트-free 기하에서 I MAPE 목표 ≤1%, |bias| ≤0.5%(v1.07 baseline I 0.66% 수준). Q/U는 RMS/I ≤2-3%(편광 관례차 허용).
- **Rrs(0−/0+)**: MAPE ≤1%.
- 비교 규약은 하네스 §1 frozen 적용. 불일치 시 양쪽 물리 검증(하네스 §1.0G — OSOAA도 틀릴 수 있음).
- 통과 = 해당 항목 OCRT 출력을 **golden CSV로 freeze** → done 전환 → 재현 set 등록.

---

## 4. 세션 시작 재현 테스트 (done 항목, byte-단위)

**목적**: 코드 변경이 이미 검증된 OCRT 출력을 깨뜨리지 않았는지 byte-단위 확인(OSOAA 불필요 — OCRT 자기 golden 대조).

**프로토콜**:
1. done 항목별 **frozen golden CSV**에서 **대표 3-5 cell**만 추출(전 케이스 아님). cell = (밴드, 기하, 변수) 한 점. 대표 추출 규칙: 기하 {sza40/vza30/raa90}(일반), {sza0/vza0}(축퇴), {sza80/vza60/raa45}(극단) + 항목 변수 양끝.
2. 동일 명령으로 OCRT 재실행 → golden과 **rho_I/Q/U 전 자리 byte-비교**(상대차 0, 또는 부동소수 last-ULP).
3. 한 cell이라도 불일치 = **재현 FAIL** → 코드 회귀. 마지막 gate-pass 체크포인트와 diff.
4. 전 cell 일치 = PASS → 그 세션 작업 진행.

**현 재현 set**: 비어 있음(full-grid golden 미생성). 유일한 기존 golden은 #2 pure-water 555(sandbox 부재). ⇒ **Jae가 golden 올리면 즉시 재현 set 등록**, 또는 항목이 done 되는 대로 golden freeze하며 채운다.

**예외(이미 sandbox 재현 확인됨, 2026-06-27)**: #1 aerosol M80C@412 moment 6 VZA가 rev3 §3과 byte 일치(재현 성공). 단 이 케이스는 **vs OSOAA 실패(−20%)**라 validation golden 아님 — 수정 후 변경되므로 regression anchor로도 보류.

---

## 5. 확정 결정 (2026-06-27 Jae)
1. **밴드**: {412,443,490,555,660,865} nm (aot 865 기준). ✅
2. **IOP 범위**: chl {0.03,0.3,3,30} mg/m³, aDOM(440) {0.01,0.1,1} m⁻¹, TSM {0.1,1,10,50} g/m³. ✅
3. **#1 black water**: **Fresnel surface 유지 + 흡수 water**(OSOAA `YS.Abs440 1000` + coxmunk; OCRT coxmunk + 흡수 water). ✅
4. **golden 생성 순서**: #1 → #2 → #3(delta-M 버그 후) → #4 → #5 → #6. ✅
5. **재현 cell**: 항목당 **5개**. ✅

## 부록: 25 distinct 기하 목록
```
sza=0:  (0,0,*) (0,30,*) (0,60,*)                          # raa 무관, 3개
sza=40: (40,0,*)                                            # nadir, raa 무관, 1개
        (40,30,0) (40,30,45) (40,30,90) (40,30,135) (40,30,180)
        (40,60,0) (40,60,45) (40,60,90) (40,60,135) (40,60,180)
sza=80: (80,0,*)                                            # nadir, raa 무관, 1개
        (80,30,0) (80,30,45) (80,30,90) (80,30,135) (80,30,180)
        (80,60,0) (80,60,45) (80,60,90) (80,60,135) (80,60,180)
```
(sza,vza,raa). `*`=raa 축퇴(대표 1개만 실행). 총 25.
