# OCRT 330-1100 nm 교체자료 적용·통합 가이드

**문서 목적:** 이 세션에서 확정된 분광자료를 최신 OCRT C/Python paired 구현에 안전하게 설치하고, 코드의 범위·보간 계약과 테스트를 일관되게 갱신하기 위한 실행 지침이다.

**인계 상태:** 자료 생성·구조검증 완료, 최신 본 OCRT tree에는 아직 미적용.

## 1. 통합 범위와 비범위

### 1.1 완료 대상으로 동결한 범위

- 파장: 330-1100 nm
- elastic atmosphere-ocean coupled RT
- C와 Python의 paired data 및 matched-input regression
- 기본 phytoplankton absorption, AHN mineral scalar/vector, atmospheric aerosol, six-gas absorption, pure-water temperature absorption correction

### 1.2 이번 완료조건에서 제외한 기능

- EAP 17종 species-specific phytoplankton phase: 계속 비활성
- Raman: excitation auxiliary domain 확장 미포함
- legacy CCRR Chl branch: 기본모델이 아니며 미확장
- O2-O2/O4 collision-induced absorption: 현재 6-gas 계약 외
- OSOAA 자체의 330 nm 확장: 수행하지 않음

## 2. 통합 전에 반드시 확인할 과학적 상태

| 구성요소 | 상태 | 통합 판정 |
|---|---|---|
| PLOPS `a_ph*` | 330-1100 nm 명시자료 | 기본 설치 |
| AHN mineral `a*`, `b*` 8개 | 330-399 nm OLS, >=400 보존 | 기본 설치 |
| AHN mineral Mie 4개 | source-informed seawater forward Mie | 설치 후 필수 RT regression |
| Detritus Mie | 330-1100 후보, legacy 잔차 존재 | 조건부, 자동설치 제외 |
| Aerosol Mie 176개 | 330·1100 explicit nodes | 기본 설치 |
| Gas xsec 6개 | 330-1100 nm, AFGL 40 level | 기본 설치 |
| `psi_T` | 330-440 exact zero, 441-1100 WOPP family | 기본 설치 |

> **중요:** PLOPS 기본값은 Diatom 단일종이 아니다. 각 분류군 내 중앙값을 계산한 뒤 6개 분류군 중앙값의 중앙값을 취한 class-balanced 대표 스펙트럼이다. 합성값과 가장 가까운 실제 시료가 Diatom CCMP1316이었을 뿐이다.

## 3. 패키지 무결성 검사

패키지 최상위에서 다음을 실행한다.

```bash
python 06_TOOLS/verify_handoff_payload.py
```

PASS 조건:

- 설치맵의 모든 C/Python payload가 존재
- 각 paired 파일의 SHA-256이 설치맵과 일치
- aerosol 176개, AHN Mie 4개, gas 6개, AHN scalar 8개 수량 확인
- C/Python 대응 파일이 byte-identical

무결성 실패 시 파일을 임의 복구하지 말고 원 component archive에서 다시 추출한다.

## 4. 최신 본 OCRT tree 준비

### 4.1 기준 tree 고정

통합 시작 전에 다음을 기록한다.

```text
OCRT C commit/version:
OCRT Python commit/version:
OSOAA reference package:
compiler / Python / OS:
baseline regression archive:
```

C와 Python은 같은 물리버전이어야 한다. 한 구현만 먼저 배포하지 않는다.

### 4.2 원 tree 백업

Git working tree가 clean한지 확인한 뒤 별도 integration branch를 만든다. 설치도구도 교체 전 파일을 패키지 내부 `local_backups/<UTC>/`에 복사하지만, Git commit/tag 백업을 대체하지 않는다.

## 5. 자료 설치 절차

### 5.1 Dry run

```bash
python 06_TOOLS/apply_data_patch.py \
  --c-root <OCRT_C_ROOT> \
  --python-root <OCRT_PY_ROOT>
```

출력의 `COPY`, `SAME`, `SKIP-CONDITIONAL`을 검토한다.

### 5.2 승인된 기본자료 적용

```bash
python 06_TOOLS/apply_data_patch.py \
  --c-root <OCRT_C_ROOT> \
  --python-root <OCRT_PY_ROOT> \
  --apply
```

기본 적용은 detritus를 제외한 196개 runtime 파일을 설치한다.

### 5.3 조건부 detritus 적용

아래 세 조건이 충족된 경우에만 실행한다.

1. 본 개발 세션이 source-recipe candidate를 새 canonical로 채택한다고 명시
2. `ORGANIC_WAVELENGTH_MIN_NM`를 330 nm로 내리는 코드수정 준비
3. C/Python matched-input I/Q/U regression 계획 확정

```bash
python 06_TOOLS/apply_data_patch.py \
  --c-root <OCRT_C_ROOT> \
  --python-root <OCRT_PY_ROOT> \
  --apply --include-detritus
```

## 6. 정확한 목적지

### 6.1 Water IOP

| 파일 | C | Python |
|---|---|---|
| `phyto_absorption_default.csv` | `inputs/water_iop/` | `data/water_iop/` |
| `psi_T_rottgers2014_OCRT.txt` | `inputs/water_iop/` | `data/water_iop/` |
| `Detritus_Stramski2001.mie` | `inputs/water_iop/` | `data/water_iop/` |

### 6.2 AHN mineral TSM

4개 `astarmin_*`, 4개 `bstarmin_*`, 4개 `*_AHN.mie`를 각각 다음에 배치한다.

```text
C      inputs/tsm_ahn/
Python data/tsm_ahn/
```

### 6.3 Atmospheric aerosol

- SnF/OPAC 16개: C `inputs/`, Python `data/`
- Ahmad paper 80개: `inputs|data/aerosol_ahmad2010_paper_mie/`
- Ahmad AccuRT 80개: `inputs|data/aerosol_ahmad2010_accurt_mie/`

전체 개별 경로와 해시는 `05_VALIDATION_AND_MANIFESTS/INSTALLATION_MAP_C_PYTHON.csv`를 단일 기준으로 사용한다.

### 6.4 Gas xsec

```text
C      inputs/xsec/xsec_{h2o,o2,co2,ch4,o3,no2}.dat
Python data/xsec/xsec_{h2o,o2,co2,ch4,o3,no2}.dat
```

## 7. 필수 코드 검토·수정

자료 복사만으로 완료되지 않는다. 최신 tree에서 이름과 위치가 달라졌을 수 있으므로 아래 검색 anchor를 기준으로 수정한다.

### 7.1 Chl-linked organic phase 하한

Detritus를 채택한 경우 C의 현재 하한을:

```c
#define ORGANIC_WAVELENGTH_MIN_NM 350.0
```

에서:

```c
#define ORGANIC_WAVELENGTH_MIN_NM 330.0
```

으로 변경한다. 관련 오류문구의 `350-1100`도 `330-1100`으로 갱신한다. Detritus를 미채택하면 이 하한을 내리면 안 된다.

### 7.2 정확한 endpoint 상태 플래그

`rt_water_iop.c`의 generic helper가 exact minimum/maximum 조회를 extrapolation으로 표시하지 않도록:

```c
xq <= x[0]      -> xq < x[0]
xq >= x[n-1]   -> xq > x[n-1]
```

로 변경한다. 값 자체의 경계보간이 아니라 상태 플래그 정정이다.

### 7.3 Silent endpoint clamp 제거 또는 상위 fail-loud 보장

다음 계열을 검색한다.

```text
rt_iop_organic.c
rt_iop_ahn_mineral.c
rt_aerosol.c
shared/numerics.c
Python interpolation / phase_batch helpers
```

helper 내부가 endpoint 값을 반환하더라도, constituent loader/evaluator는 요청 파장이 선언된 table/phase 범위 밖이면 명확히 실패해야 한다. 이번 패치의 330-1100 지원은 자료로 구현하며, 범위 밖 clamp를 과학모델로 사용하지 않는다.

### 7.4 Gas reader

`xsec_*.dat`는 40 x 771이다. 파장 수 751 또는 350 nm 시작을 하드코딩한 배열·검사·문서가 없는지 확인한다. loader가 header/첫 행에서 동적으로 파장 수를 읽고 330, 1100 nm exact endpoint를 허용해야 한다.

### 7.5 EAP 비활성 유지

EAP generator/catalog의 350-850 nm 제한은 이번 production 범위의 blocker가 아니다. 그러나 EAP가 사용자 옵션 또는 default 경로로 우회 활성화되지 않도록 C/Python 모두 smoke test를 추가한다.

### 7.6 문서·버전

- `SPECTRAL_SUPPORT_330_1100.json`을 최종 계약과 일치하도록 갱신
- README의 지원범위와 예외목록 갱신
- version string / release notes 갱신
- 모든 파일의 source DOI, project policy, caveat 유지

## 8. 구성요소별 필수 단위시험

### 8.1 PLOPS `a_ph*`

- 330, 350, 443, 555, 750, 775, 799, 800, 1100 nm
- 800-1100 exact zero
- `TChl-a=1 mg m-3`에서 수치적 `a_ph=a_ph*`
- C/Python exact parity

### 8.2 AHN scalar/vector

- 4종 x 330/350/400/750/1100 nm
- scalar 399-400 join과 >=400 preservation
- phase `P11>0`, `abs(P12)<=P11`, `abs(P33)<=P11`
- half-integral normalization, `g`, `bb/b`

### 8.3 Aerosol 176종

- 176/176 load
- 330·1100 node 직접조회
- 350·860·1240 legacy node 비교
- 대표 16+80+80 family smoke

### 8.4 Gas xsec

- 6/6 shape = 40 x 771
- finite/nonnegative
- 350-1100 canonical preservation
- 330/340/349/350 boundary
- O2/CO2/CH4 330-349 explicit zero

### 8.5 `psi_T`

- 330-440 exact zero
- 442-1100 even source nodes exact preservation
- `T=20 C`에서 base `a_w` exact
- 0, 5, 10, 20, 25, 30 C에서 `a_w>0`

## 9. RT acceptance regression

### 9.1 공통 파장

```text
330, 340, 349, 350, 400, 412, 443, 490, 555,
620, 680, 709, 745, 865, 900, 940, 1000, 1100 nm
```

### 9.2 구성요소 분리시험

1. Pure water, T=20 C, gas OFF
2. Pure water, T=0/10/30 C
3. Chl-only: 0.1, 1, 3 mg m-3
4. AHN mineral-only: 4종, low/high TSM
5. Detritus candidate on/off (채택 시)
6. Rayleigh + each aerosol representative
7. Gas-only: O3, NO2, H2O and all gases
8. Coupled atmosphere-ocean

출력은 최소 `Rrs/rrs I,Q,U`, TOA I/Q/U, Ed/Eu, Kd를 포함한다.

### 9.3 기존범위 regression 해석

- Gas 350-1100은 exact preservation 대상이다.
- Aerosol SnF/OPAC·AccuRT 기존 node는 exact; Ahmad paper family의 작은 residual은 사용자 승인범위이다.
- PLOPS와 AHN Mie는 과학모델 자체가 교체되므로 old RT exact match 대상이 아니다. 대신 새 C/Python parity, 물리성, reference trend, 변경영향 기록이 acceptance 기준이다.
- pure-water, CDOM, surface, solver 등 비변경 경로는 기존 regression을 유지해야 한다.

## 10. 권장 통과기준

| 검증 | 통과기준 |
|---|---|
| 파일구조·해시 | 100% PASS |
| C/Python data | byte-identical |
| C/Python matched RT | 부동소수 허용범위 내 일치, 차이 원인 문서화 |
| phase physicality | `P11>0`; `abs(P12)<=P11`; `abs(P33)<=P11` |
| phase normalization | 파일별 기존 validator 기준 통과 |
| gas finite/nonnegative | 전 층·전 파장 PASS |
| 20 C pure-water | 기존 base exact |
| out-of-range | silent clamp 없이 명확한 오류 |
| EAP/Raman | 비활성·비지원 상태 명시 |

## 11. Rollback

1. integration commit 이전 tag/branch로 복귀하거나 installer backup을 복원한다.
2. C/Python을 반드시 함께 되돌린다.
3. partial rollback을 금지한다. 예: scalar만 구버전, Mie만 신버전인 혼합상태 금지.
4. regression artifact와 실패조건을 보존한다.

## 12. 최종 release 체크리스트

- [ ] 196개 기본자료 paired 설치
- [ ] Detritus 채택 여부 명시
- [ ] 필요 시 detritus 포함 197개 설치
- [ ] code range / endpoint / fail-loud 정리
- [ ] C/Python build·unit·RT tests PASS
- [ ] 330/340/349/350 boundary plots
- [ ] 350-1100 intentional/non-intentional 변화 분리 보고
- [ ] source DOI와 project policy를 data header에 유지
- [ ] release note와 migration doc 갱신
- [ ] 최종 package SHA-256 생성
