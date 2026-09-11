# C LUT 예시 도구

[make_c_luts.py](make_c_luts.py)는 **Python 표준 라이브러리로 C 실행 파일을 호출**한다. OCRT Python 물리 구현을 실행하는 도구가 아니다. 설치된 Python 3과 빌드된 OCRT C가 필요하다.

전체 옵션은 `python make_c_luts.py --help`로 본다. 기본 격자는 테스트용으로 작다. 생산 정확도나 최적 격자를 뜻하지 않는다.

현재 해양 출력기의 배열 크기에 맞춰 simulation의 RAA 표본은 512개 이하로 제한한다. 이 제약은 방위각 간격과 관련되며 내부 Fourier 차수와는 다르다.

저장소 최상위에서 먼저 계획을 만든다.

```bash
python3 docs/manual/examples/make_c_luts.py rayleigh \
  --ocrt-root code/OCRT_C --exe code/OCRT_C/build/ocrt \
  --outdir task/rayleigh_demo
```

생성된 `jobs.csv`와 `manifest.json`을 확인한 뒤 같은 명령에 `--run`을 붙인다. `--run` 없이도 버전 확인을 위해 `--version`만 호출하지만 복사전달 계산은 수행하지 않는다. Windows는 `python3` 대신 `python`, 실행 파일은 `code/OCRT_C/build/ocrt_x86-64-v3.exe`를 사용한다.

| 옵션 | 의미·기본값 |
|---|---|
| `rayleigh` / `simulation` | 필수 위치 인자. 해면 포함 Rayleigh / 결합 해양 CDOM 모의자료 |
| `--ocrt-root`, `--exe`, `--outdir` | 필수. C 작업 폴더, 실행 파일, 결과 폴더. 상대 경로는 호출 당시 폴더 기준 |
| `--wavelengths` | 쉼표 구분 nm, 기본 `443,555` |
| `--szas` | 쉼표 구분 도, 기본 `30` |
| `--pressures` | 쉼표 구분 hPa, 기본 `1013.25` |
| `--winds` | 쉼표 구분 m/s, 기본 `5` |
| `--aods` | simulation의 AOD555, 기본 `0`; 양수에는 `--mie` 필수 |
| `--adom440` | simulation의 aDOM440(m⁻¹), 기본 `0.1`; Chl/TSM은 이 도구에서 0 고정 |
| `--mie` | 에어로졸 파일. 상대 경로는 C 작업 폴더 기준 |
| `--vza-step`, `--vza-max`, `--raa-step` | 도 단위, 기본 `30`, `60`, `90`; VZA 간격은 최대각을, RAA 간격은 360을 정확하게 나누도록 선택 |
| `--pssa` | rayleigh에 IPSS 적용. 결합 해양에는 사용 불가 |
| `--include-direct-glint` | rayleigh의 직달 선글린트 포함. 기본은 decouple. simulation의 해양 격자 TOA는 직달 선글린트가 분리된 값이며 해당 항은 CSV에 별도 기록되지 않음 |
| `--gas-off` | simulation의 6기체 총량을 0으로 설정. rayleigh는 항상 0 |
| `--threads` | 각 프로세스 OMP 스레드 수, 기본 `1`; 사례는 순차 실행 |
| `--timeout` | 사례당 제한 초, 기본 `900` |
| `--run` | 실제 계산. 기존 사례 CSV가 있으면 덮어쓰지 않고 종료 |

각 사례는 CSV, 표준출력, 표준에러 로그를 저장한다. 행 수·주요 관측량의 유한수·해수 수렴 플래그를 검사하지만 수렴 정확도나 I/Q/U 대칭까지 자동 보장하지 않는다. 상향 전달비의 유효성은 사용자가 `T_up_rt_valid`로 확인한다. 실패한 폴더를 그대로 재개하지 않으며 새 폴더를 사용한다. 외부에서 설정한 `OCRT_*` 환경변수는 상속되고 manifest에 기록되므로, 진단·실험 변수가 켜진 셸인지 확인한다.

`jobs.csv`는 **이 도구의 계획표**이며 C의 `--batch` 또는 `--batch-full-grid` 입력 형식이 아니다. 전체 C 인자 배열은 `manifest.json`의 `argv`에 있다. 실행 파일 해시 외에도 사용 소스 커밋과 입력 자료 해시는 따로 보관한다.

## 핵심 사용 예시 10개 / 10 core worked examples

[run_recipes.py](run_recipes.py)는 LaTeX 매뉴얼 부록 C의 10가지 예시(13개 C 호출)를 선택·계획·실행한다. [recipe_catalog.json](recipe_catalog.json)은 **현재 매뉴얼에 남긴 C 인자의 실행 목록**이다. ID는 R01, R02, R07, R21, R26, R28, R36, R38, R39, R40이며 기본 `--select all`도 이 10개만 포함한다. 설명을 고칠 때는 `latex/ko/09_cookbook.tex`와 `latex/en/09_cookbook.tex`를 직접 수정하고, 명령을 바꾸면 카탈로그도 같은 내용으로 수정한다. Python 표준 라이브러리 외의 패키지는 필요 없다.

The runner executes the 10 retained examples (13 C invocations) documented in Appendix C; it does not invoke the OCRT Python solver. Its default `--select all` includes only these 10 examples. The Korean and English LaTeX sources are the authoritative narrative. Keep their command blocks and this catalogue synchronized when updating an example. The 30 deferred candidates are listed in the [review record](../review/README.md) and are not selectable in the current catalogue.

저장소 최상위에서 / From the repository root:

```bash
python3 docs/manual/examples/run_recipes.py --list
python3 docs/manual/examples/run_recipes.py \
  --select R01,R02,R07,R38 \
  --ocrt-root code/OCRT_C --exe code/OCRT_C/build/ocrt \
  --outdir task/cookbook_selected --threads 1 --timeout 900
```

계획 단계에서는 실행 파일도 호출하지 않는다. `manifest.json`을 확인한 뒤 같은 명령에 `--run`을 추가한다. Windows에서는 `python`과 Windows용 실행 파일 경로를 사용한다.

Planning never calls the executable. Review `manifest.json`, then repeat with `--run`. Native Windows uses `python` and a Windows executable; no Bash array or redirection is needed.

| 옵션 / Option | 동작 / Behavior |
|---|---|
| `--list` | ID·그룹·영문 제목을 표시. 다른 경로 인자는 불필요 / List IDs, groups and titles without running. |
| `--select` | 쉼표 구분 ID 또는 `all`; 기본 `all` / Comma-separated IDs or `all` (default). |
| `--group` | geometry, rayleigh, aerosol, ocean, workflow 중 쉼표 구분 선택; `--select`와 교집합 / Filter the selection by one or more groups. |
| `--ocrt-root` | C 작업 폴더 / C working directory, normally `code/OCRT_C`. |
| `--exe` | 빌드된 C 실행 파일 / Built C executable. |
| `--outdir` | 새 결과 폴더 또는 동일한 미실행 계획 폴더 / New output folder or the identical unexecuted plan folder. |
| `--timeout` | C 호출당 제한 초; 기본 900 / Seconds per C invocation, default 900. |
| `--threads` | C 프로세스 내부 OMP 스레드 수; 기본 1. 사례는 순차 실행 / Threads within each C process, default 1; cases run serially. |
| `--run` | 자료·기존 파일 확인 후 실제 실행 / Execute after preflight. |

CLI의 상대 경로는 도구를 호출한 폴더 기준이다. 카탈로그 내부 자료 경로는 C 작업 폴더 기준이다. `{out}`은 절대 결과 폴더로 바꾼다. R39는 BOM 없는 UTF-8 네이티브 배치 입력을 자동 생성하며 경로의 쉼표·따옴표·줄바꿈을 거부한다. 자료가 없으면 실제 실행 전에 필요한 파일 목록을 보여준다.

CLI paths are relative to the invocation directory; catalogue input paths are relative to the C working directory. The `{out}` token becomes the absolute output directory. R39 creates its native UTF-8 batch fixture without a BOM; batch output paths cannot contain commas, quotes or newlines. Missing required data are reported before execution.

이 도구는 상속된 `OCRT_*` 환경변수를 지우고 사례에 명시된 설정만 활성화한다. 소스 기준 커밋, 실행 파일·필수 자료·카탈로그 해시, 전체 인자, 환경 설정, 시간, 종료 코드, 검사 결과를 manifest에 기록한다. 기존 stdout·stderr·CSV를 덮어쓰지 않으며 실행한 폴더의 재개는 지원하지 않는다. 새 폴더를 사용한다.

The runner removes inherited `OCRT_*` settings and enables only documented per-case values. Its manifest stores the documented source commit, executable/input/catalogue hashes, complete arguments, environment, timing, exit code and basic checks. Existing stdout, stderr or CSV files are never overwritten. Use a new folder after any execution; resuming a partially executed folder is unsupported.

주요 관측량·격자 행/열/각도·수렴 플래그를 검사하지만 생산 정확도를 보증하지 않는다. 대기 단일 출력의 `conv`는 Fourier SOS 수렴, 결합 해양 단일 출력의 `conv`는 반환된 수중 계산 수렴을 보고한다. 해양 CSV는 `water_converged`를 검사한다. `T_up_rt_valid=0`은 별도로 기록한다. [현재 검증 요약](verification_core_examples.json)은 2026-09-11에 실행한 명령 중 유지한 13건을 선별한 기록이다. 원래 40개·49회 및 강한 기체흡수의 층 수 민감도 등의 과거 진단은 [검토 보관본](../review/recipes_2026-09-12/verification_40_2026-09-11.json)에 남긴다.

Basic checks cover primary finite observables, grid dimensions/angles and convergence flags, not production accuracy. Atmospheric single-text `conv` aggregates Fourier SOS convergence; coupled-ocean single-text `conv` reports the returned water solve. Ocean CSVs use `water_converged`. Invalid upward transmittance is recorded separately. The [current verification summary](verification_core_examples.json) selects 13 unchanged invocations executed on 2026-09-11. The original 40-example/49-invocation record, including historical strong-absorption diagnostics, is preserved in the [review archive](../review/recipes_2026-09-12/verification_40_2026-09-11.json).
