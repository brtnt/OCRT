# Ahn four-species mineral (TSM) optical data

This directory is the canonical OCRT mineral / TSM data set introduced in v1.15.

Species and vector-phase files:

- `red_clay` → `Red_clay_AHN.mie` (default)
- `brown_earth` → `Brown_earth_AHN.mie`
- `yellow_clay` → `Yellow_clay_AHN.mie`
- `calcareous_sand` → `Calcareous_sand_AHN.mie`

The four `.mie` files are the files supplied for this integration. The `astarmin_*` and
`bstarmin_*` audit tables are derived from their bulk spectral blocks without resampling:

- `a* = Extinct_Co − Scatter_Co`
- `b* = Scatter_Co`

OCRT uses these mass-specific coefficients with the dry-weight TSM concentration in g m⁻³. The
backscatter ratio is integrated from the matching P11 phase function. Set `OCRT_TSM_DIR` to
override this directory.

---

# Ahn 4종 광물(TSM) 광학 자료

이 디렉터리는 v1.15 에서 도입한 OCRT 표준 광물/TSM 자료다.

종과 벡터 위상 파일:

- `red_clay` → `Red_clay_AHN.mie` (기본값)
- `brown_earth` → `Brown_earth_AHN.mie`
- `yellow_clay` → `Yellow_clay_AHN.mie`
- `calcareous_sand` → `Calcareous_sand_AHN.mie`

`.mie` 파일 네 개는 이번 통합에 공급된 파일 그대로다. `astarmin_*` 과 `bstarmin_*` 감사표는 벌크 분광
블록에서 재표본 없이 유도했다.

- `a* = Extinct_Co − Scatter_Co`
- `b* = Scatter_Co`

OCRT 는 이 질량비 계수를 건조중량 TSM 농도(g m⁻³)와 함께 쓴다. 후방산란 비는 대응하는 P11 위상함수에서
적분한다. `OCRT_TSM_DIR` 로 이 디렉터리를 바꿀 수 있다.
