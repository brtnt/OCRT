# EAP 17종 coated-sphere phase generator provenance — 2026-07-26

## 공개 API

```c
int rt_eap_mie_phase_compute(
    rt_eap_species_id_t species_id,
    const double *wavelength_nm,
    size_t n_wavelength,
    const double *theta_deg,
    size_t n_theta,
    double *p11,
    double *p12,
    double *p33);
```

- species ID: 0–16
- 운영 파장: 350–850 nm
- 출력: wavelength-major P11/P12/P33
- 정규화: `0.5*integral(P11 sin(theta) dtheta)=1`
- caller-owned memory, 파일 I/O 없음, re-entrant

## 구현 파일

- `src/rt_eap_mie_phase.{c,h}`
- `src/internal/rt_eap_coated_mie.{c,h}`
- `src/internal/rt_eap_species_catalog.{c,h}`
- `src/generated/rt_eap_species_data.inc`
- `tools/eap_mie_writer.c`
- `tools/eap_generate_all.sh`
- `tools/generate_eap_species_catalog.py`
- `tools/eap_catalog_freeze.py`

## 원자료 상태

동결 catalog에는 원자료 해시가 기록되어 있으나 원자료 파일 자체는 이번 전달
패키지에 없다.

```text
EAP_invivo_means.csv
92051b77cec90744294536dbfe611b65350532b98d7d0f8c32545bb5e4f273a8

501nm_extended_e1701000.mat
50638c22d4596c38da5f34e10d846188faa225d6fb262409f858fce6e9725c02
```

따라서 `rt_eap_species_data.inc`에서 phase를 계산하는 단계는 재현되지만, 원자료에서
catalog를 다시 생성하는 단계는 원자료가 확보될 때까지 provenance gap으로 남는다.

## 생성 결과

`tools/eap_generate_all.sh representative inputs/water_iop`의 17개 대표 `.mie`를
생성하였다. species 11/Deff 0.5 µm는 독립 2회 생성 결과와 배포 파일이 모두
byte-identical하였다. 전체 파일 해시는
`validation/eap_generation/EAP_17_SPECIES_SHA256_20260726.txt`에 기록하였다.

## 생산 RT 적용상태

17종 generator/API와 `.mie` 산출물은 통합되었으나, OCRT constituent water model의
species-specific phytoplankton scattering은 비활성화되어 있다. 현재 L=200 moment
representation에서 강한 forward peak를 가진 대부분 종의 reconstructed P11이
중간각에서 음수가 될 수 있기 때문이다.

- `--ocrt-phyto-group` + Chl>0: 종료코드 2
- species 미지정 + Chl>0: phytoplankton absorption만 적용
- `b_phyto=bb_phyto=0`
- Chl 연계 입자 phase: detritus only

복원 조건은 component-level truncation과 b-rescaling을 구현하여 모듈 폐합을
확보하거나, 충분히 높은 phase expansion order를 정식 지원하는 것이다.
