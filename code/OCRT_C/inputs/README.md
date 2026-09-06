# OCRT runtime inputs

Run OCRT from the `ocrt/` directory or through `../bin/ocrt`, so that these relative paths
resolve. The official data catalog and provenance documents are at the package-root `docs/`.

- `afgl_atm/`: atmospheric profiles
- `xsec/`: gas cross-section tables
- `water_iop/`: pure water, OCRT organic and CCRR Chl data
- `tsm_ahn/`: Ahn four-species mineral data
- top-level `.mie`: aerosol phase tables (and retained legacy tables)
- `deprecated/`: excluded diagnostic data

The 181 Mie tables (`*.mie`) are not stored in git; install them with `scripts/fetch_data.sh` or
`scripts/fetch_data.ps1` from the repository root (GitHub Release `data-v1`).

---

# OCRT 실행 입력 자료

상대 경로가 맞게 풀리도록 OCRT 는 `ocrt/` 디렉터리에서 실행하거나 `../bin/ocrt` 로 실행한다. 공식 자료
목록과 출처 문서는 패키지 루트의 `docs/` 에 있다.

- `afgl_atm/`: 대기 프로파일
- `xsec/`: 기체 흡수 단면적 표
- `water_iop/`: 순수해수, OCRT 유기물, CCRR Chl 자료
- `tsm_ahn/`: Ahn 4종 광물 자료
- 최상위 `.mie`: 에어로졸 위상 표(그리고 보존된 옛 표)
- `deprecated/`: 제외된 진단 자료

Mie 표 181개(`*.mie`)는 git 에 저장하지 않는다. 저장소 루트에서 `scripts/fetch_data.sh` 또는
`scripts/fetch_data.ps1` 로 설치한다(GitHub Release `data-v1`).
