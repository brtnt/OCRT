# OCRT 330–1100 nm FR631 분광자료 적용·통합 가이드

**문서 버전:** 2026-08-20  
**대상:** paired OCRT C/Python 본 개발 트리  
**전제:** elastic OCRT, EAP 기본 비활성, Raman 제외

## 1. 목적과 완료 범위

본 인계 세트는 2026-08-14 분광자료 인계 ZIP의 용량 축소/누락 문제를 대체한다. 이번 세트의 핵심은 기존 361각·희소 phase-wavelength 표현을 폐기하고, 모든 canonical Mie 파일을 아래 단일 계약으로 재생성한 것이다.

- Mie 파일 수: 198개
- bulk spectral table: 330–1100 nm, 1 nm, 771행
- phase wavelength columns: 330–1100 nm, 1 nm, 771열
- phase elements: P11, P12, P33
- angle grid: FR631 exact 631점, 180°→0° 내림차순
- 각도 계산: FR631 각도에서 직접 계산; 성긴 각도 결과의 후처리 각도보간 금지
- C/Python 소비계약: 파장·각도 선형 보간

198개 구성은 SNF/OPAC 16, Ahmad2010 paper 80, Ahmad2010 AccuRT 80, AHN 광물 4, detritus 1, EAP 17이다. `Phytoplankton.mie`, `TSM_mineral.mie`, `Calcareous_sand_from_OSOAA_phase_443_550.mie`는 legacy alias이므로 canonical 대상에서 제외한다.

## 2. 필수 패키지

PART00–PART10을 모두 같은 상위 폴더에 풀어 하나의 bundle root로 병합한다. PART11은 생성 재현 자료이며 runtime 설치에는 필수가 아니다.

추출 후 다음 파일이 있어야 한다.

```text
<BUNDLE>/README_FIRST_KO.md
<BUNDLE>/02_RUNTIME_MIE/aerosol/...
<BUNDLE>/02_RUNTIME_MIE/water/...
<BUNDLE>/02_RUNTIME_NON_MIE/...
<BUNDLE>/03_INTEGRATION_PATCH/...
<BUNDLE>/06_MANIFESTS/MIE_CANONICAL_INVENTORY_198.csv
```

전수 검증:

```bash
python <BUNDLE>/05_TOOLS/verify_extracted_bundle.py
```

예상 결과:

```text
expected=198
failure_count=0
status=PASS
```

## 3. 데이터 설치 전에 반드시 적용할 코드 패치

### 3.1 C reader line buffer

771개의 phase wavelength 값과 angle 1개를 한 줄에 기록하면 한 행이 기존 8192-byte buffer를 초과한다. `src/shared/mie_io.c`와 `mie_io.h`의 line buffer를 131072 byte 이상으로 확장해야 한다. 이 패치 없이 대용량 phase row가 잘리거나 parser가 phase block을 찾지 못한다.

### 3.2 각도 표현

FR631은 전방에서 매우 조밀한 piecewise-linear representation이다. 따라서 C/Python 모두 다음을 지켜야 한다.

- phase angle interpolation: linear
- Legendre/P2 moment용 densification: linear
- bb/b 적분: native FR631 grid에서 직접 trapezoid
- PCHIP/cubic angular interpolation 사용 금지

### 3.3 파장 표현

파일에 정수 파장 330–1100 nm가 모두 명시되어 있다. 정수 파장은 exact column lookup과 동일해야 하며, fractional wavelength만 adjacent columns 사이 선형 보간한다.

### 3.4 유기입자 범위

C `ORGANIC_WAVELENGTH_MIN_NM/MAX_NM` 및 Python 대응 검사를 330/1100 nm로 변경한다. 범위 밖 silent clamp는 제거하거나 fail-loud로 변경한다.

### 3.5 Python cache key

Python `bb_b_ratio()` cache는 `id(mie)` 대신 `(source_path, mtime_ns, size)` 기반 안정 key를 사용한다. 객체 ID 재사용에 의한 stale cache 가능성을 제거한다.

패치 위치:

```text
03_INTEGRATION_PATCH/OCRT_FR631_C_PYTHON_INTEGRATION.patch
03_INTEGRATION_PATCH/C/...
03_INTEGRATION_PATCH/Python/...
```

## 4. 설치 절차

### 4.1 dry-run

```bash
python <BUNDLE>/05_TOOLS/apply_data_patch.py \
  --c-root <OCRT_C_ROOT> \
  --python-root <OCRT_PY_ROOT>
```

### 4.2 실제 설치

```bash
python <BUNDLE>/05_TOOLS/apply_data_patch.py \
  --c-root <OCRT_C_ROOT> \
  --python-root <OCRT_PY_ROOT> \
  --apply
```

설치기는 canonical source를 C/Python에 복사한 뒤 SHA-256을 다시 비교한다. Python AccuRT 80종은 canonical subdirectory와 기존 flat lookup 위치에 같은 bytes를 함께 배치한다.

### 4.3 목적지 요약

| 자료 | C | Python |
|---|---|---|
| SNF/OPAC | `inputs/*.mie` | `data/*.mie` |
| Ahmad paper | `inputs/aerosol_ahmad2010_paper_mie/` | `data/aerosol_ahmad2010_paper_mie/` |
| Ahmad AccuRT | `inputs/aerosol_ahmad2010_accurt_mie/` | `data/aerosol_ahmad2010_accurt_mie/` + flat alias |
| AHN 광물 | `inputs/tsm_ahn/` | `data/tsm_ahn/` |
| detritus | `inputs/water_iop/` | `data/water_iop/` |
| EAP | `inputs/water_iop/eap/` | `data/water_iop/eap/` |
| gas xsec | `inputs/xsec/` | `data/xsec/` |
| PLOPS/psi_T | `inputs/water_iop/` | `data/water_iop/` |

정확한 215개 logical source mapping은 `06_MANIFESTS/INSTALLATION_MAP_C_PYTHON.csv`를 따른다.

## 5. 성능·메모리 주의

Mie 파일 한 개는 약 25 MB이며 198개 총 uncompressed size는 약 4.94 GB이다. 모든 파일을 동시에 메모리에 올리면 안 된다.

- active aerosol model 1개, active water constituent만 load/cache
- process worker 간 mutable cache 공유 금지
- I/O는 case hot loop 밖에서 한 번만 수행
- 같은 파일의 771 wavelength × 631 angle 배열은 read-only로 유지
- multi-geometry 계산은 한 SOS run에서 유지

## 6. 검증 순서

### Gate A — 구조

198개 모두 다음을 만족해야 한다.

- first declaration = 631
- angle text exact match
- bulk rows = 771
- phase blocks = 3
- phase rows/block = 631
- phase wavelength columns = 771
- range = 330–1100 nm, step 1 nm

### Gate B — 물리성

- P11 > 0
- |P12| ≤ P11
- |P33| ≤ P11
- 0 ≤ SSA ≤ 1
- Scatter = Extinction × SSA
- 1/2 ∫P11 dμ = 1
- bulk g = 1/2 ∫μP11 dμ

### Gate C — loader/build

- full C release build warning 0
- C reader representative 6 families PASS
- Python reader representative 6 families PASS
- C/Python same file and wavelength return exact/near-machine parity

### Gate D — RT smoke

우선 파장:

```text
330, 340, 350, 400, 443, 550, 865, 1100 nm
```

우선 모델:

```text
M50C
Ahmad paper r80f20
Ahmad AccuRT r80f20
Brown earth AHN
Detritus_Stramski2001
EAP_00 (advanced option only; default remains disabled)
```

검사량:

- aerosol-only TOA I/Q/U
- pure-water Rrs/rrs I/Q/U
- Chl-only
- each AHN TSM-only
- coupled atmosphere-ocean

### Gate E — regression

기존 350–1100 nm와 비교할 때 bulk IOP가 보존되는 계열은 exact 또는 machine-level을 요구한다. Phase representation은 FR631 및 1 nm column으로 변경되므로, 결과 차이는 새 angular/spectral representation에 기인하는지 분리한다.

## 7. 승인 기준

- Mie structural audit: 198/198 PASS
- phase normalization max error ≤ 2×10⁻⁹
- g phase-vs-bulk max error ≤ 1×10⁻⁹
- Mueller bound violation = 0
- C/Python matched-input parity 통과
- 20°C pure-water baseline exact 보존
- EAP default disabled 유지
- Raman 제외를 문서에 명시

## 8. 롤백

적용 전 C/Python tree 전체를 immutable archive 또는 Git tag로 보존한다. 문제가 발생하면 데이터와 코드 패치를 함께 롤백해야 한다. FR631 파일을 기존 cubic/PCHIP consumer로 읽는 혼합상태는 허용하지 않는다.
