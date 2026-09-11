# OCRT — Ocean Colour vector Radiative Transfer

Repository: https://github.com/brtnt/OCRT (public; academic and non-commercial license)

**Documentation:** [directory guide](docs/README.md) · [Korean PDF](docs/manual/OCRT_Manual_KO.pdf) · [English PDF](docs/manual/OCRT_Manual_EN.pdf) · [LaTeX sources](docs/manual/latex/README.md)

OCRT is a coupled atmosphere–ocean vector radiative transfer code for ocean-colour
atmospheric correction (GOCI-II / GOCI-III). It solves the polarized radiative transfer
equation with the successive-orders-of-scattering (SOS) method. The package has two parts:
a C reference implementation and a Python package that reproduces the same physics for
batch data production. The spectral range is 330–1100 nm.

Current version: **OCRT C v1.11.1 (2026-09-05)**. The spherical-shell correction now uses the
IPSS scheme of Zhai & Hu (2022) only; the older PSSA scheme was removed. The viewing zenith
angle is anchored at the surface pixel, and the correction factor κ is weighted by the
constituent phase functions in mixed atmospheres. Details:
`OCRT_MIGRATION_STATUS_2026-08-23.md` §14–§14.1 and
`code/OCRT_C/validation/ipss_2026-09-05/OCRT_IPSS_REPLACEMENT_RECORD_2026-09-05.md`.

## Repository layout

| Path | Contents |
|---|---|
| `code/OCRT_C/` | C reference implementation (`src/`), gate and regression tests (`tests/`, `scripts/`), diagnostic tools (`tools/`), technical documents (`docs/`), validation records (`validation/`), change logs (`CHANGES_*.md`, `RELEASE_NOTES_*.md`), file manifest (`SHA256SUMS.txt`) |
| `code/OCRT_Python/` | Python batch package (`ocrt_py/`, `produce_grid.py`, `ocrt_solve.py`, GPU gate) and small input data |
| `code/campaign_runner/`, `code/run_tools/`, `code/validation_05/` | Campaign runner, run tools, validation harness and bit baselines |
| `docs/technical/`, `docs/validation/` | Technical reports and validation records |
| `docs/reports/` | Session reports (validation documents and figures) |
| `docs/manual/` | User manual: installation, all CLI options, RAA, Rayleigh LUT and simulation examples, output formats |
| `scripts/` | Data installation scripts (`fetch_data.sh` / `.ps1`) and SHA-256 manifests |
| `OCRT_MIGRATION_STATUS_2026-08-23.md` | Cumulative development status |
| `MIGRATION_README_2026-09-05_KO.md` | **How to start a new work session; inventory of what is where** |

## Data installation (Mie tables, required for aerosol and hydrosol runs)

The 181 Mie tables (FR631 format, 4.3 GB) are not stored in git. They are distributed as four
zip files (1.3 GB) in the GitHub Release
[`data-v1`](https://github.com/brtnt/OCRT/releases/tag/data-v1).
Run one command after cloning; it installs the tables into `code/OCRT_C/inputs/` and
`code/OCRT_Python/data/`.

```
bash scripts/fetch_data.sh                                          # Linux / macOS
powershell -ExecutionPolicy Bypass -File scripts\fetch_data.ps1     # Windows
```

SHA-256 of the zip files: `scripts/fetch_data.sha256`. SHA-256 of each of the 181 tables:
`scripts/MIE_SHA256SUMS_data-v1.txt`. Rayleigh-only runs and the IPSS gates do not need the
Mie tables. All other small input data (AFGL atmospheres, gas cross sections, IOP tables, test
fixtures) are included in the repository.

Not tracked in git: build products (`build/`, `*.exe`), campaign outputs, delivery archives,
and the unpublished manuscript.

## Build

Linux (gcc ≥ 11 with OpenMP):

```
cd code/OCRT_C
make                      # Makefile default: -march=cascadelake (AVX-512)
# On a CPU without AVX-512:
gcc -std=c11 -O3 -march=x86-64-v3 -ffp-contract=fast -fassociative-math -fno-signed-zeros \
    -fno-trapping-math -DOCRT_FAST_KERNELS -fopenmp -Isrc $(find src -name '*.c') -o build/ocrt -lm
```

Windows: `code\OCRT_C\build_win.bat` (see `REBUILD_WIN_v1.11.1_KO.txt`; uses `-march=x86-64-v3`).

Gates: `cd code/OCRT_C && bash scripts/run_ipss_gates.sh`. They run without Mie tables and must
print `IPSS GATES: ALL PASS` and `IPSS REGRESSION: ALL PASS`.

## Run example

```
cd code/OCRT_C
OCRT_ADVANCED=1 ./build/ocrt --wavelength 555 --sza 40 --vza 55 --raa 90 \
    --surface black_fresnel_ocean --wind-speed 5 --mie inputs/M80C.mie --aod-555 0.1 \
    --n-layers 400 --pssa
```

`--pssa` switches the IPSS spherical-shell correction on (there is no `--pssa-mode` since v1.11).
Relative azimuth convention: 180° is the specular (sun-glint) direction and 0° is the
backscattering direction.

## Citation and license

Academic and non-commercial license (`LICENSE`). Commercial use requires a separate agreement
(brtnt@kiost.ac.kr). Copyright (c) 2026 Jae-Hyun Ahn, Korea Ocean Satellite Center (KOSC) / KIOST.

Spherical-shell correction: Zhai, P.-W., Hu, Y. (2022). An improved pseudo spherical shell
algorithm for vector radiative transfer. *JQSRT* 282, 108132.

---

# OCRT — 해색 벡터 복사전달 코드

저장소: https://github.com/brtnt/OCRT (공개, 학술·비상업 라이선스)

**문서:** [디렉토리 안내](docs/README.md) · [한글판 PDF](docs/manual/OCRT_Manual_KO.pdf) · [영문판 PDF](docs/manual/OCRT_Manual_EN.pdf) · [LaTeX 원본](docs/manual/latex/README.md)

OCRT 는 해색 위성(GOCI-II / GOCI-III) 대기보정을 위한 대기–해양 결합 벡터 복사전달 코드다.
연속차수산란(SOS) 방법으로 편광 복사전달 방정식을 푼다. 패키지는 두 부분으로 구성된다.
C 참조 구현과, 같은 물리를 배치 자료생산용으로 재현한 Python 패키지다. 파장 범위는 330–1100 nm 다.

현재 버전: **OCRT C v1.11.1 (2026-09-05)**. 구면 보정은 Zhai & Hu (2022) 의 IPSS 방식만 사용하며,
이전의 PSSA 방식은 삭제했다. 관측천정각은 지표 화소를 기준으로 하고, 혼합 대기에서는 보정계수 κ 를
성분별 위상함수로 가중한다. 상세는 `OCRT_MIGRATION_STATUS_2026-08-23.md` §14–§14.1 과
`code/OCRT_C/validation/ipss_2026-09-05/OCRT_IPSS_REPLACEMENT_RECORD_2026-09-05.md` 를 참조한다.

## 저장소 구성

| 경로 | 내용 |
|---|---|
| `code/OCRT_C/` | C 참조 구현(`src/`), 게이트·회귀 테스트(`tests/`, `scripts/`), 진단 도구(`tools/`), 기술 문서(`docs/`), 검증 기록(`validation/`), 변경 이력(`CHANGES_*.md`, `RELEASE_NOTES_*.md`), 파일 매니페스트(`SHA256SUMS.txt`) |
| `code/OCRT_Python/` | Python 배치 패키지(`ocrt_py/`, `produce_grid.py`, `ocrt_solve.py`, GPU 게이트)와 소형 입력 자료 |
| `code/campaign_runner/`, `code/run_tools/`, `code/validation_05/` | 캠페인 실행기, 실행 도구, 검증 하니스와 비트 기준선 |
| `docs/technical/`, `docs/validation/` | 기술 보고서와 검증 기록 |
| `docs/reports/` | 작업 세션 산출 기록(검증 문서·그림) |
| `docs/manual/` | 사용자 매뉴얼: 설치, 전체 CLI 옵션, RAA, Rayleigh LUT·시뮬레이션 예시, 출력 형식 |
| `scripts/` | 자료 설치 스크립트(`fetch_data.sh` / `.ps1`)와 SHA-256 매니페스트 |
| `OCRT_MIGRATION_STATUS_2026-08-23.md` | 개발 상태 문서(누적) |
| `MIGRATION_README_2026-09-05_KO.md` | **새 작업 세션 시작 절차와 자료 위치 목록** |

## 자료 설치 (Mie 표 — 에어로졸·하이드로졸 계산에 필요)

Mie 표 181개(FR631 형식, 4.3 GB)는 git 에 저장하지 않는다. GitHub Release
[`data-v1`](https://github.com/brtnt/OCRT/releases/tag/data-v1) 의 zip 4개(1.3 GB)로 배포한다.
클론 후 아래 명령을 한 번 실행하면 `code/OCRT_C/inputs/` 와 `code/OCRT_Python/data/` 에 설치된다.

```
bash scripts/fetch_data.sh                                          # Linux / macOS
powershell -ExecutionPolicy Bypass -File scripts\fetch_data.ps1     # Windows
```

zip 파일의 SHA-256 은 `scripts/fetch_data.sha256`, 표 181개 각각의 SHA-256 은
`scripts/MIE_SHA256SUMS_data-v1.txt` 에 있다. Rayleigh 단독 계산과 IPSS 게이트는 Mie 표 없이
실행된다. 그 밖의 소형 입력 자료(AFGL 대기, 기체 흡수 단면적, IOP 표, 테스트 픽스처)는 저장소에
들어 있다.

git 에 넣지 않은 것: 빌드 산출물(`build/`, `*.exe`), 캠페인 산출물, 배포 묶음, 투고 전 원고.

## 빌드

Linux (gcc ≥ 11, OpenMP):

```
cd code/OCRT_C
make                      # Makefile 기본값: -march=cascadelake (AVX-512)
# AVX-512 가 없는 CPU 에서는:
gcc -std=c11 -O3 -march=x86-64-v3 -ffp-contract=fast -fassociative-math -fno-signed-zeros \
    -fno-trapping-math -DOCRT_FAST_KERNELS -fopenmp -Isrc $(find src -name '*.c') -o build/ocrt -lm
```

Windows: `code\OCRT_C\build_win.bat` (`REBUILD_WIN_v1.11.1_KO.txt` 참조, `-march=x86-64-v3` 사용).

게이트: `cd code/OCRT_C && bash scripts/run_ipss_gates.sh`. Mie 표 없이 실행되며
`IPSS GATES: ALL PASS` 와 `IPSS REGRESSION: ALL PASS` 가 출력되어야 한다.

## 실행 예

```
cd code/OCRT_C
OCRT_ADVANCED=1 ./build/ocrt --wavelength 555 --sza 40 --vza 55 --raa 90 \
    --surface black_fresnel_ocean --wind-speed 5 --mie inputs/M80C.mie --aod-555 0.1 \
    --n-layers 400 --pssa
```

`--pssa` 는 IPSS 구면 보정을 켜는 스위치다(v1.11 부터 `--pssa-mode` 는 없다).
상대방위각 규약: 180° 가 경면(선글린트) 방향, 0° 가 후방산란 방향이다.

## 인용과 라이선스

학술·비상업 라이선스(`LICENSE`). 상업적 이용은 별도 협의가 필요하다(brtnt@kiost.ac.kr).
Copyright (c) 2026 Jae-Hyun Ahn, Korea Ocean Satellite Center (KOSC) / KIOST.

구면 보정: Zhai, P.-W., Hu, Y. (2022). An improved pseudo spherical shell algorithm for vector
radiative transfer. *JQSRT* 282, 108132.
