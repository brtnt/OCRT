# OCRT — Ocean Colour vector Radiative Transfer

저장소: https://github.com/brtnt/OCRT (공개, 학술·비상업 라이선스)

해색 위성(GOCI-II / GOCI-III) 대기보정을 위한 **대기–해양 결합 벡터 복사전달 코드**다.
연속차수산란(SOS) 기반의 C 참조 구현과, 같은 물리를 배치 자료생산용으로 재현한 Python 패키지로
구성된다. 정확도 기준(reference)은 OSOAA 이며, 파장 범위는 330–1100 nm(≤ 2000 nm 정책)다.

Coupled atmosphere–ocean vector radiative transfer (successive orders of scattering) for
ocean-colour atmospheric correction. C reference implementation plus a batch-oriented Python
port; validated against OSOAA.

현재 버전: **OCRT C v1.11.1 (2026-09-05)** — 구면 보정을 Zhai & Hu (2022) IPSS 로 단일화(구 PSSA 삭제),
화소 기준 VZA 앵커, 혼합 대기 위상함수 가중 κ. 상세는 `OCRT_MIGRATION_STATUS_2026-08-23.md` §14–§14.1 과
`code/OCRT_C/validation/ipss_2026-09-05/OCRT_IPSS_REPLACEMENT_RECORD_2026-09-05.md`.

## 저장소 구성

| 경로 | 내용 |
|---|---|
| `code/OCRT_C/` | C 참조 구현(`src/`), 게이트·회귀 테스트(`tests/`, `scripts/`), 진단 도구(`tools/`), 기술 문서(`docs/`), 검증 기록(`validation/`), 변경 이력(`CHANGES_*.md`, `RELEASE_NOTES_*.md`), 파일 매니페스트(`SHA256SUMS.txt`) |
| `code/OCRT_Python/` | Python 배치 패키지(`ocrt_py/`, `produce_grid.py`, `ocrt_solve.py`, GPU 게이트) |
| `code/OSOAA/` | OSOAA 참조 구성 기록·패치(diff)·매니페스트 (OSOAA 본체는 포함하지 않음) |
| `code/campaign_runner/`, `code/run_tools/`, `code/validation_05/` | 캠페인 실행기, 실행 도구, 검증 하니스·비트 기준선 |
| `docs/technical/`, `docs/validation/` | 기술 보고서, 검증 기록 |
| `docs/reports/` | 세션 산출 기록(검증 문서·그림) |
| `OCRT_MIGRATION_STATUS_2026-08-23.md` | 마이그레이션 상태 문서(누적) |
| `MIGRATION_README_2026-09-05_KO.md` | **새 세션 부트스트랩 절차와 인벤토리** |

## 저장소에 없는 것 (로컬 디스크에만 있음)

GitHub 용량 정책상 다음은 `.gitignore` 로 제외했다. 위치와 복구 절차는 `MIGRATION_README_2026-09-05_KO.md` §3.

- Mie 산란 자료 `*.mie` (182 파일, 4.3 GB; `code/OCRT_C/inputs/`) — 배포본은 `package/PART01–PART10*.zip`
- `code/OCRT_Python/data/` (4.3 GB), `code/validation_05/runs/` (OSOAA 참조 런 원문)
- `package/`(배포 묶음), `task/`(캠페인 산출물), `docs/paper/`(투고 전 원고)
- 빌드 산출물(`build/`, `*.exe`)

## 빌드

Linux (gcc ≥ 11, OpenMP):

```
cd code/OCRT_C
make                      # Makefile 기본 -march=cascadelake (AVX-512)
# AVX-512 가 없는 CPU 에서는:
gcc -std=c11 -O3 -march=x86-64-v3 -ffp-contract=fast -fassociative-math -fno-signed-zeros \
    -fno-trapping-math -DOCRT_FAST_KERNELS -fopenmp -Isrc $(find src -name '*.c') -o build/ocrt -lm
```

Windows: `code\OCRT_C\build_win.bat` (`REBUILD_WIN_v1.11.1_KO.txt` 참조, `-march=x86-64-v3`).

게이트: `cd code/OCRT_C && bash scripts/run_ipss_gates.sh` — Mie 자료 없이 실행되며
`IPSS GATES: ALL PASS` / `IPSS REGRESSION: ALL PASS` 가 나와야 한다.

## 실행 예

```
cd code/OCRT_C
OCRT_ADVANCED=1 ./build/ocrt --wavelength 555 --sza 40 --vza 55 --raa 90 \
    --surface black_fresnel_ocean --wind-speed 5 --mie inputs/M80C.mie --aod-555 0.1 \
    --n-layers 400 --pssa
```

`--pssa` 는 IPSS 구면 보정 on/off 스위치다(v1.11 부터 `--pssa-mode` 없음). RAA 규약: 180° 가 경면(글린트) 방위, 0° 가 후방산란.

## 인용·라이선스

학술·비상업 라이선스(`LICENSE`). 상업적 이용은 별도 협의(brtnt@kiost.ac.kr).
Copyright (c) 2026 Jae-Hyun Ahn, Korea Ocean Satellite Center (KOSC) / KIOST.

구면 보정 알고리즘: Zhai, P.-W., Hu, Y. (2022). An improved pseudo spherical shell algorithm for
vector radiative transfer. *JQSRT* 282, 108132.
