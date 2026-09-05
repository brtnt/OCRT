# OCRT 330–1100 nm 분광자료 통합 명세

## 1. 기준선과 적용 대상

| 항목 | 기준 |
|---|---|
| C 기준선 | `OCRT-v1.2-2026-08-09-KST-chl-zero-extension-legacy-mie-cleanup` |
| Python 기준선 | `pyOCRT-v1.2-2026-08-09-chl-zero-extension-legacy-mie-cleanup` |
| C 최종본 | `OCRT-v1.2-2026-08-14-KST-spectral-data-330-1100-component-gated` |
| Python 최종본 | `pyOCRT-v1.2-2026-08-14-spectral-data-330-1100-component-gated` |
| OSOAA | 현행 2026-07-28 reference 유지; 이번 작업에서 미수정 |
| 공식 elastic envelope | 330–1100 nm |
| Raman | 제외 |
| EAP species scattering | 비활성 유지 |

## 2. 설치자료 수량

handoff paired mapping은 197개였다.

| 구성요소 | handoff 수량 | production 채택 |
|---|---:|---:|
| SnF/OPAC aerosol Mie | 16 | 16 |
| Ahmad2010 AccuRT aerosol Mie | 80 | 80 |
| Ahmad2010 paper aerosol Mie | 80 | 80 |
| AHN mineral Mie | 4 | 4 |
| AHN mineral scalar `a*`,`b*` | 8 | 8 |
| PLOPS-derived `a_ph*` | 1 | 1 |
| pure-water `psi_T` | 1 | 1 |
| gas xsec | 6 | 6 |
| detritus Mie | 1 conditional | 0 active replacement; candidate archived |
| **합계** | **197** | **196 accepted + 1 rejected candidate archive** |

모든 accepted 파일은 C/Python에서 byte-identical하고 handoff payload와 SHA-256이 동일하다.

## 3. 목적지

### C

```text
inputs/
inputs/aerosol_ahmad2010_paper_mie/
inputs/aerosol_ahmad2010_accurt_mie/
inputs/tsm_ahn/
inputs/xsec/
inputs/water_iop/
inputs/water_iop/candidates/
```

### Python

```text
data/
data/aerosol_ahmad2010_paper_mie/
data/aerosol_ahmad2010_accurt_mie/
data/tsm_ahn/
data/xsec/
data/water_iop/
data/water_iop/candidates/
```

## 4. 데이터별 최종 상태

### 4.1 Pure water

- 파일: `water_coef_z09_1nm.txt`
- 범위: 200–2449 nm, 1 nm
- 열: wavelength, `a_w`, total `b_w`
- OCRT: `bb_w=0.5*b_w`
- handoff/첨부/C/Python SHA-256:
  `f173e1a4514b653973def109a5415678afbe4f39564a7b19f94a71ef55ea0e46`
- 판정: 이미 canonical과 동일, 수치 교체 없음

### 4.2 Pure-water temperature correction

- 파일: `psi_T_rottgers2014_OCRT.txt`
- 범위: 330–1100 nm, 1 nm
- 330–440 nm: exact zero
- 441–1100 nm: handoff의 Röttgers/WOPP-family PCHIP 계열
- 범위 밖 endpoint hold 금지

### 4.3 Phytoplankton absorption

- 파일: `phyto_absorption_default.csv`
- 범위: 330–1100 nm, 1 nm
- 330–750 nm: PLOPS class-balanced representative spectrum
- 750–800 nm: 800 nm exact zero로 선형 taper
- 800–1100 nm: explicit zero rows
- runtime `if (lambda > cutoff) a=0` 제거
- phytoplankton production scattering: 계속 0

### 4.4 Chl-linked detritus

- active: 기존 `Detritus_Stramski2001.mie`, 350–850 nm
- rejected candidate: `candidates/Detritus_Stramski2001_330_1100_CANDIDATE_REJECTED_L200.mie`
- candidate는 handoff 파일과 byte-identical하게 보존
- active range 밖 Chl>0 요청은 fail-loud

### 4.5 AHN mineral TSM

- scalar 8개와 Mie 4개 적용
- target: 330–1100 nm
- 330–399 nm scalar: approved 400–420 nm 구간 OLS
- 400 nm 이상 scalar: 보존
- vector Mie: source-informed seawater forward-Mie restoration
- exact historical reproduction을 주장하지 않음

### 4.6 Atmospheric aerosol

- 총 176 `.mie`
- explicit 330 and 1100 nodes
- 176개 모두 C/Python paired exact
- active Mie physicality/range audit 통과

### 4.7 Gas absorption

- H2O, O2, CO2, CH4, O3, NO2
- 330–1100 nm, 1 nm, 771 points
- AFGL 40 levels
- 350–1100 canonical overlap 보존
- O2-O2/O4 CIA 제외

## 5. 코드변경

### 공통 계약

- official envelope: 330–1100 nm
- table interpolation: in-range only
- out-of-range: fail-loud
- silent endpoint clamp를 과학적 확장으로 사용 금지
- exact 330/1100 endpoints는 in-range

### C

- 신규 `src/rt_spectral_contract.h`
- `main.c`, `rt_solver.c`, `rt_water_rt.c`에 global envelope 적용
- `rt_water_iop.c`: pure water/psi table range check와 exact endpoint 처리
- `rt_iop_organic.c/.h`: data-driven PLOPS, detritus native-range API 및 fail-loud
- `rt_iop_ahn_mineral.c`: scalar/Mie range validation
- `shared/mie_io.c/.h`, `rt_aerosol_runtime.c`: target/reference/phase table coverage 검증
- `rt_absorption.c`: 771-point gas axis 검증, out-of-range 금지

### Python

- 신규 `ocrt_py/spectral_contract.py`
- `constituent.py`: pure water, psi, PLOPS, AHN, detritus query range 검증
- `aerosol.py`: bulk/phase target/reference range 검증
- `absorption.py`: gas table coverage와 query 검증
- `driver.py`, `atmos.py`: official envelope 적용
- EAP production scattering gate 유지

## 6. 삭제·격리

- 과거 legacy root-level narrow Mie는 이전 2026-08-09 기준선에서 이미 제거된 상태를 유지
- Python root legacy `r*.mie` 중복 80개 제거
- rejected detritus candidate는 active filename에서 제거하고 candidates 디렉터리에 격리
- 이전 350–1100 spectral manifest는 `docs/historical_spectral_contracts/`로 이동

## 7. machine-readable 기준

- `water_iop/SPECTRAL_SUPPORT_330_1100.json`
- `validation/INSTALLATION_DECISION_MAP_C_PYTHON.csv`
- `validation/ACCEPTED_DATA_SHA256SUMS.txt`
- `validation/validation_summary.json`
