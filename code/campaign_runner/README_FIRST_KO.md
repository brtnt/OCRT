# OCRT polarization-sensitivity data-generation run package

**Package version:** 2026-08-22
**Run policy:** PSSA ON · 16 independent workers · 1 OpenMP thread per worker · latest raw FR631 ·
broad truncation OFF

## 1. Fixed run conditions

```text
PSSA                         ON: --pssa is given explicitly in every real OCRT command
Process workers              16
OMP_NUM_THREADS per worker   1
Water particle phase         latest exact FR631, direct linear kernel
Broad water truncation       OFF
OCRT_DEBUG                    not used
OCRT_ADVANCED                 1
```

Sixteen independent cases run at the same time. Instead of giving 16 OpenMP threads to one case,
16 independent OCRT processes are started and each is fixed to `OMP_NUM_THREADS=1`. On a 24-core
PC therefore at most 16 cases run at once.

`--n-mu-water` is an advanced option. This package does not insert a value for water runs on its
own; if the matrix does not give a value for every water row, the runner stops.

## 2. Important limitation

The total count and the common grid of the past campaign were confirmed, but the exact
8,045-run command matrix is not present in the current environment. The following values were
not guessed:

- the exact list of the 24 aerosol models and the central model
- the 7 AOD865 values
- the 64 water states
- the conditions of the 170 wind runs
- the advanced numerical options of each coupled run

Hence the following can run immediately:

- the pure-Rayleigh PSSA reference: 5 bands × 3 SZA = 15 runs
- any campaign after an archived command file has been converted into a matrix

With the empty `HISTORICAL_RUN_MATRIX_REQUIRED.csv` unfilled, a full campaign stops (fail-loud).

## 3. Package layout and size

The distribution separates the core and the data parts to avoid download failures.

```text
Core ZIP                       about 200 MiB
PART00–PART10 compressed       about 1.41 GiB
total transfer                 about 1.60 GiB
final runtime                  about 5.5–6 GiB
recommended free space         8 GiB for the runtime + space for results
```

The uncompressed size of the 198 canonical Mie tables is about 4.94 GB. The installer extracts
the PARTs in parallel to a temporary location, hard-links the files into the C / Python
destinations and deletes the temporary names, so the C / Python data are not stored twice; during
installation extra space is needed for the temporary extraction. The recommended 8 GiB includes
this peak.

The output size of a full campaign depends strongly on the columns and cases, so it is not fixed
here. Once the exact matrix is available, run the first 16 cases with `--pilot 16` to measure the
real mean CSV.GZ size and time.

## 4. File placement

Unpack the core ZIP and place the following 11 files **in its parent folder or in the parts folder
you specify**:

```text
PART00 ... PART10
```

PART11 is generator reproduction material and is not needed for runtime computation.

## 5. Windows 10/11 recommended run: WSL2

In PowerShell, from the core package folder:

```powershell
.\runner\SETUP_AND_BUILD_WSL.ps1
.\runner\RUN_PURE_RAYLEIGH_16CORES_WSL.ps1
```

The first command verifies the archive SHA-256 → extracts the migration package → applies the
policy overlay → applies the FR631 consumer patch → streams the 198 data files into place → builds
a release per CPU → runs the PSSA full-grid smoke test.

## 6. Direct run on Linux / WSL

```bash
./runner/setup_and_build_linux.sh <PARTS_DIR> ./runtime
./runner/run_pure_rayleigh_16cores_linux.sh ./runtime ./results/pure_rayleigh_5band_pssa
```

After the full matrix is available:

```bash
python3 runner/run_campaign.py \
  --runtime-root runtime/MIGRATION_PKG_2026-08-19 \
  --matrix matrices/HISTORICAL_RUN_MATRIX.csv \
  --output-dir results/polarization_full \
  --workers 16 --pilot 16
```

After checking the pilot, remove `--pilot 16` and resume.

## 7. Result layout

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

Each `.csv.gz` is the authoritative source. The runner renames a file atomically only after the
structure, row count, geometry grid and finite-value checks pass. On a rerun, valid results are
reused.

## 8. Restoring a past package

When the original `OCRT_paper_migration_2026-08-04.tar.gz` is available:

```bash
python3 runner/inspect_historical_package.py --archive OCRT_paper_migration_2026-08-04.tar.gz
python3 runner/build_matrix_from_run_commands.py \
  --commands historical_import/.../run_commands.txt \
  --output matrices/HISTORICAL_RUN_MATRIX.csv
```

The conversion stops if PSSA is missing, broad truncation was used, or a water row lacks
`--n-mu-water`.

---

# OCRT 편광 민감도 자료생성 실행 패키지

**패키지 버전:** 2026-08-22
**실행 정책:** PSSA ON · 독립 worker 16개 · worker 당 OpenMP 1 thread · 최신 raw FR631 · broad truncation OFF

## 1. 이번에 고정한 실행 조건

```text
PSSA                         ON: 모든 실제 OCRT 명령에 --pssa 명시
Process workers              16
OMP_NUM_THREADS per worker   1
Water particle phase         최신 exact FR631, 직접 선형 커널
Broad water truncation       OFF
OCRT_DEBUG                    사용하지 않음
OCRT_ADVANCED                 1
```

독립 케이스 16개를 동시에 실행한다. 케이스 하나에 OpenMP 16개를 주는 방식이 아니라, 서로 독립인 OCRT
프로세스 16개를 띄우고 각 프로세스를 `OMP_NUM_THREADS=1` 로 고정한다. 따라서 24코어 PC 에서는 최대 16개의
케이스만 동시에 실행된다.

`--n-mu-water` 는 고급 옵션이다. 이 패키지는 물 실행에 값을 임의로 넣지 않는다. 행렬의 모든 물 행에 값이
명시되어 있지 않으면 실행기가 중단한다.

## 2. 중요한 제한

과거 캠페인의 전체 수량과 공통 격자는 확인되었으나, 정확한 8,045회 실행 명령 행렬은 현재 실행 환경에
없다. 다음 값은 추정해 넣지 않았다.

- 에어로졸 모델 24개의 정확한 목록과 중심 모델
- AOD865 값 7개
- 해수 상태 64개
- 풍속 실행 170건의 조건
- 결합 실행별 고급 수치 옵션

따라서 다음은 즉시 실행할 수 있다.

- 순수 Rayleigh PSSA 기준: 5 밴드 × 3 SZA = 15회
- 보관된 명령 파일을 행렬로 변환한 뒤의 임의 캠페인

빈 `HISTORICAL_RUN_MATRIX_REQUIRED.csv` 를 채우지 않은 상태에서 전체 캠페인은 fail-loud 로 중단한다.

## 3. 패키지 구성과 용량

이 배포는 다운로드 실패를 피하기 위해 코어와 자료 파트를 분리한다.

```text
Core ZIP                       약 200 MiB
PART00–PART10 압축             약 1.41 GiB
총 전송량                      약 1.60 GiB
최종 실행 환경                 약 5.5–6 GiB
권장 여유 공간                 실행 환경용 8 GiB + 결과 저장 공간
```

표준 Mie 198개의 비압축 크기는 약 4.94 GB 다. 설치기는 PART 를 병렬로 임시 추출한 뒤 C/Python 목적지에
hardlink 하고 임시 추출본을 삭제한다. 따라서 C/Python 자료를 중복 저장하지 않지만, 설치 중에는 임시
추출을 위한 추가 공간이 필요하다. 권장 여유 공간 8 GiB 는 이 최대치를 포함한다.

전체 캠페인 출력 용량은 열 수와 케이스에 따라 크게 달라지므로 고정하지 않는다. 정확한 행렬이 준비되면
처음 16개를 `--pilot 16` 으로 실행해 실제 CSV.GZ 평균 크기와 시간을 잰다.

## 4. 파일 배치

코어 ZIP 을 풀고, **그 상위 폴더 또는 지정한 parts 폴더**에 다음 11개를 둔다.

```text
PART00 ... PART10
```

PART11 은 생성기 재현 자료이며 실행 계산에는 필요하지 않다.

## 5. Windows 10/11 권장 실행: WSL2

PowerShell 에서 코어 패키지 폴더로 이동한 뒤:

```powershell
.\runner\SETUP_AND_BUILD_WSL.ps1
.\runner\RUN_PURE_RAYLEIGH_16CORES_WSL.ps1
```

첫 명령은 archive SHA-256 확인 → 마이그레이션 패키지 추출 → 정책 overlay → FR631 소비 패치 → 자료 198개
스트리밍 설치 → CPU 별 release 빌드 → PSSA 전체격자 스모크 테스트를 수행한다.

## 6. Linux/WSL 직접 실행

```bash
./runner/setup_and_build_linux.sh <PARTS_DIR> ./runtime
./runner/run_pure_rayleigh_16cores_linux.sh ./runtime ./results/pure_rayleigh_5band_pssa
```

전체 행렬이 확보된 후:

```bash
python3 runner/run_campaign.py \
  --runtime-root runtime/MIGRATION_PKG_2026-08-19 \
  --matrix matrices/HISTORICAL_RUN_MATRIX.csv \
  --output-dir results/polarization_full \
  --workers 16 --pilot 16
```

pilot 확인 후 `--pilot 16` 을 빼고 재개 실행한다.

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

각 `.csv.gz` 가 정본이다. 실행기는 구조·행 수·기하 격자·유한값 검사를 통과한 뒤에만 원자적으로 이름을
바꾼다. 재실행 시 유효한 결과는 재사용한다.

## 8. 과거 패키지 복원

원본 `OCRT_paper_migration_2026-08-04.tar.gz` 가 확보되면:

```bash
python3 runner/inspect_historical_package.py --archive OCRT_paper_migration_2026-08-04.tar.gz
python3 runner/build_matrix_from_run_commands.py \
  --commands historical_import/.../run_commands.txt \
  --output matrices/HISTORICAL_RUN_MATRIX.csv
```

PSSA 누락, broad truncation 사용, 물 행의 `--n-mu-water` 누락이 있으면 변환을 중단한다.
