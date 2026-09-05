# OCRT 2단계 EAP 17종 generator 통합 및 phytoplankton scattering 비활성화 검토보고서

**기준일:** 2026-07-26  
**기준 소스:** `OCRT-v1.2-2026-07-25-KST-stage2-coupling-clamp-fix`  
**최종 버전:** `OCRT-v1.2-2026-07-26-KST-stage2-eap-generator-phyto-scattering-disabled`

## 1. 판정

첨부 작업지시서의 핵심 수정방향을 현재 누적 OCRT 기준선에 선택적으로 통합하였다.
기존에 폐합된 대기, water RAA, FIX1+FIX2+FIX3, FIX4, exact-pole 및 coupling-clamp
모듈은 변경하지 않았다.

적용 결과는 다음과 같다.

1. 17종 EAP coated-sphere phase generator와 단일 공개 API를 통합하였다.
2. 17종 대표 `.mie`를 생성하여 배포 입력자료에 포함하였다.
3. constituent water model의 species-specific EAP scattering을 비활성화하였다.
4. Chl>0에서 명시적 `--ocrt-phyto-group`은 조용히 대체하지 않고 종료코드 2로 실패한다.
5. species 미지정 Chl 실행은 phytoplankton absorption을 유지하고 `b_phyto=bb_phyto=0`으로 둔다.
6. Chl 연계 입자산란은 detritus 하나만 사용한다.
7. constituent truncation 관련 플래그는 미구현 상태를 조용히 무시하지 않고 fail-loud 처리한다.
8. fixed-bulk IOP truncation 경로는 그대로 유지한다.
9. 비-Chl CDOM/TSM 18조건은 구버전과 stdout/stderr가 byte-identical하다.

## 2. 주요 코드 변경

### 2.1 EAP generator

신규 공개 경계는 다음과 같다.

```c
species ID + wavelength array + theta array
    -> normalized P11/P12/P33
```

생성기는 Chl, a/b/bb, bb/b, RT moments 또는 파일경로를 공개 인자로 받지 않는다.
미세물리 catalog와 수치 Mie 설정은 내부 구현으로 고정한다.

### 2.2 EAP scattering 비활성화

`organic_phyto_scattering` gate의 생산 기본값은 0이다. 내부 diagnostic re-enable
필드는 남겼지만 production CLI에서는 활성화할 수 없다.

Chl>0에서 species를 명시하면 다음 사유를 포함하여 실패한다.

- current in-water phase representation: L≤200
- strongly forward-peaked species의 mid-angle reconstructed P11 음수 가능성
- validated component-level truncation 또는 higher-order representation이 아직 없음

species를 명시하지 않으면 선택된 기본 micro table의 absorption spectrum만 사용한다.
phase ratio와 phytoplankton phase moments는 계산하지 않는다. detritus phase만 준비한다.

### 2.3 필요한 자료만 초기화

기존 compatibility initializer는 catalog-wide 진단을 위해 20개 phyto table을 모두
읽는다. 생산 단일 Chl 실행에는 `rt_iop_organic_init_selected()`를 사용하여 선택된
absorption table과 detritus만 읽는다. 이는 비활성화된 17종 phase cache의 불필요한
런타임 I/O를 방지한다.

### 2.4 미구현 truncation fail-loud

다음 플래그는 constituent model에서 작동하지 않았으므로 오류로 전환하였다.

```text
--ocrt-mie-truncation
--ocrt-mie-ss-mode
```

사용자는 검증된 fixed-bulk 경로를 사용해야 한다.

```text
--water-model iop
--iop-mie-phase FILE
--iop-mie-truncation
```

## 3. EAP phase API 검증

species 11(Prochlorococcus), 412/443 nm, 0–180° 0.5° 격자에서 다음을 확인하였다.

- phase normalization
- `P11>=0`
- `|P12|<=P11`
- `|P33|<=P11`
- 반복호출 bit identity
- invalid species/range error codes

결과: PASS.

17종 대표 `.mie`가 모두 생성되었다. Deff 24 µm 두 종은 각각 약 195초가 소요되었고,
나머지는 약 0.17–28초였다. species 11을 독립 2회 재생성한 파일과 배포 파일의
SHA-256은 모두 동일하였다.

## 4. constituent model gate 검증

다음을 자동시험하였다.

- canonical 17종 이름, legacy pico/nano/micro를 Chl>0과 명시하면 rc=2
- species를 생략한 Chl=0.3에서 phytoplankton `a>0`, `b=bb=0`
- detritus `b>0`, `bb>0`
- phase component는 detritus 하나만 존재
- constituent truncation flags는 fail-loud
- fixed-bulk truncation path는 성공

결과: PASS.

## 5. 비회귀

### 5.1 비-Chl byte regression

CDOM 12조건과 TSM 6조건, 총 18조건에서 구 기준 바이너리와 최종 바이너리의
stdout 및 stderr가 모두 byte-identical하였다.

### 5.2 동결 모듈

다음 집중 회귀를 통과하였다.

- public water RAA
- FIX1+FIX2+FIX3
- FIX4 pure-Rayleigh closure: `2.776e-17`
- FIX4 Rayleigh+aerosol closure: `5.551e-17`
- exact-pole continuity 및 spin-2 covariance
- coupling clamp n_mu_water=48/64/96
- single/full-grid parity
- direct-glint RAA anchors

대기 소스 및 대기 하네스는 수정하지 않았다.

## 6. 제공 수중 정합자료 검증

11,154행 원자료에서 50행 summary를 다시 계산하였다. 모든 수치열의 최대 절대차는
0이었다. Rrs/rrs I/Q/U 18개 산포도도 독립 재생성하였다.

수면 아래 rrs I 평균오차 범위는 제공자료에서 다음과 같다.

```text
CDOM  0.0346–0.3409%
TSM   0.0446–0.4943%
Chl   0.1321–0.2483%
```

단, Chl 결과는 production constituent path의 직접 종단 결과가 아니다. 전달
`ff_chl2.py`는 scalar IOP를 추출한 뒤 detritus-only truncation, pure-water Rayleigh
재혼합 및 b-rescaling을 수행하는 fixed-bulk 우회경로를 사용한다. 따라서 위 Chl
수치는 fixed-bulk validation evidence로 분류한다.

수정된 OSOAA 실행파일과 완전한 provenance가 포함되지 않아 OSOAA 원시 실행은 이번
통합에서 새로 수행하지 않았다. 전달 CSV의 내부 통계와 코드 경로는 검증하였으나,
OSOAA 외부 독립 재현은 별도 provenance 확보 후 수행해야 한다.

## 7. 성능

동일 물리인 무-Chl CDOM 단일기하를 구/신 바이너리로 각각 5회 교차 측정하였다.

```text
old mean   0.228257 s
new mean   0.226025 s
mean       -0.978%
old median 0.227057 s
new median 0.226742 s
median     -0.139%
```

유의한 계산시간 증가는 없다. 5% 병목조사 기준에 해당하지 않는다.

Chl=0.3 absorption-only 단일기하의 최종 실행시간은 6.11초였다. 이 값은 detritus
phase projection을 포함한 절대 실행시간이며, 구버전은 EAP scattering까지 포함하여
물리가 다르므로 성능회귀 비율로 직접 사용하지 않는다.

## 8. provenance 제한

동결 generated catalog는 포함되었으나 다음 원자료는 누락되었다.

```text
EAP_invivo_means.csv
501nm_extended_e1701000.mat
```

catalog 주석에 원자료 SHA-256은 남아 있다. 따라서 catalog→phase/.mie는 재현되지만
raw source→catalog 단계는 현재 패키지 하나만으로 재현되지 않는다.

## 9. 향후 복원 조건

species-specific EAP scattering은 다음 중 하나가 독립 모듈 검증을 통과한 뒤에만
복원한다.

1. component-level forward-peak truncation
   - removed mass
   - b rescaling
   - truncated bb/b 재계산
   - truncation-aware cache key
   - pure-water Rayleigh 비절단
2. 충분히 높은 phase expansion order의 정식 옵션화와 수렴검증

종단 OCRT–OSOAA 오차가 일시적으로 감소한다는 이유만으로 gate를 해제하지 않는다.
