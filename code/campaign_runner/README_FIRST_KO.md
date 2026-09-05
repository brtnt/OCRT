# OCRT 편광 민감도 자료생성 실행 패키지

**패키지 버전:** 2026-08-22  
**실행 정책:** PSSA ON · 16 independent workers · worker당 OpenMP 1 thread · latest raw FR631 · broad truncation OFF

## 1. 이번에 고정한 실행 조건

```text
PSSA                         ON: 모든 실제 OCRT 명령에 --pssa 명시
Process workers              16
OMP_NUM_THREADS per worker   1
Water particle phase         latest exact FR631, direct linear kernel
Broad water truncation       OFF
OCRT_DEBUG                    사용하지 않음
OCRT_ADVANCED                 1
```

16개의 독립 case를 동시에 실행한다. 한 case 내부에 OpenMP 16개를 주는 방식이 아니라, 서로 독립인 OCRT 프로세스 16개를 띄우고 각 프로세스는 `OMP_NUM_THREADS=1`로 고정한다. 따라서 24-core PC에서 최대 16개의 계산 case만 동시에 실행된다.

`--n-mu-water`는 advanced option이다. 이 패키지는 water run에 값을 임의로 넣지 않는다. 모든 water row의 matrix에 값을 명시하지 않으면 runner가 중단한다.

## 2. 중요한 제한

과거 campaign의 전체 수량과 공통 격자는 확인되었으나, 정확한 8,045-run command matrix는 현재 실행환경에 없다. 다음 값을 추정해 넣지 않았다.

- 24 aerosol model 정확한 목록과 중앙 model
- AOD865 7개 값
- water state 64개 값
- 170개 wind run 조건
- coupled run별 advanced numerical option

따라서 다음은 즉시 실행 가능하다.

- 5 bands × 3 SZA의 pure-Rayleigh PSSA reference 15 runs
- 정확한 archived command file을 matrix로 변환한 뒤의 임의 campaign

빈 `HISTORICAL_RUN_MATRIX_REQUIRED.csv`를 채우지 않은 상태에서 full campaign은 fail-loud로 중단한다.

## 3. 패키지 구성과 용량

이 배포는 다운로드 실패를 피하기 위해 core와 data parts를 분리한다.

```text
Core ZIP                       약 200 MiB
PART00–PART10 compressed       약 1.41 GiB
총 전송량                      약 1.60 GiB
최종 runtime                   약 5.5–6 GiB
권장 여유공간                  runtime용 8 GiB + 결과 저장공간
```

Canonical 198개 Mie의 비압축 크기는 약 4.94 GB이다. 설치기는 PART를 병렬로 임시 추출한 뒤 C/Python 목적지에 hardlink하고 임시 추출명을 삭제한다. 따라서 C/Python data 중복 저장을 피하지만, 설치 중에는 임시 추출을 위한 추가 여유공간이 필요하다. 권장 여유공간 8 GiB는 이 peak를 포함한다.

전체 campaign 출력용량은 열 수와 case에 따라 크게 달라지므로 고정 추정하지 않는다. 정확한 matrix가 준비되면 처음 16개를 `--pilot 16`으로 실행해 실제 CSV.GZ 평균 크기와 시간을 산정한다.

## 4. 파일 배치

Core ZIP을 풀고, **그 상위 폴더 또는 지정한 parts 폴더**에 다음 11개를 둔다.

```text
PART00 ... PART10
```

PART11은 generator 재현자료이며 runtime 계산에는 필요하지 않다.

## 5. Windows 10/11 권장 실행: WSL2

PowerShell에서 core package 폴더로 이동한 뒤:

```powershell
.\runner\SETUP_AND_BUILD_WSL.ps1
.\runner\RUN_PURE_RAYLEIGH_16CORES_WSL.ps1
```

첫 명령은 archive SHA-256 확인 → migration 추출 → policy overlay → FR631 consumer patch → 198개 data streaming install → CPU별 release build → PSSA full-grid smoke를 수행한다.

## 6. Linux/WSL 직접 실행

```bash
./runner/setup_and_build_linux.sh <PARTS_DIR> ./runtime
./runner/run_pure_rayleigh_16cores_linux.sh ./runtime ./results/pure_rayleigh_5band_pssa
```

Full matrix가 확보된 후:

```bash
python3 runner/run_campaign.py \
  --runtime-root runtime/MIGRATION_PKG_2026-08-19 \
  --matrix matrices/HISTORICAL_RUN_MATRIX.csv \
  --output-dir results/polarization_full \
  --workers 16 --pilot 16
```

Pilot 확인 후 `--pilot 16`을 제거하고 resume 실행한다.

## 7. 결과 구조

```text
results/<campaign>/
  run_config.json
  run_commands.txt
  run_manifest.csv
  RUN_SUMMARY.json
  PILOT_OR_RUN_ESTIMATE.json
  logs/*.stdout.log
  logs/*.stderr.log
  <stage>/*.csv.gz
```

각 `.csv.gz`가 authoritative source이다. runner는 구조·행수·기하격자·finite 검사를 통과한 뒤에만 atomic rename한다. 재실행 시 유효한 결과는 재사용한다.

## 8. 과거 package 복원

원본 `OCRT_paper_migration_2026-08-04.tar.gz`가 확보되면:

```bash
python3 runner/inspect_historical_package.py --archive OCRT_paper_migration_2026-08-04.tar.gz
python3 runner/build_matrix_from_run_commands.py \
  --commands historical_import/.../run_commands.txt \
  --output matrices/HISTORICAL_RUN_MATRIX.csv
```

PSSA 누락, broad truncation 사용, water row의 `--n-mu-water` 누락이 있으면 변환을 중단한다.
