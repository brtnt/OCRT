# OCRT user manual

**LaTeX is the authoritative source.**

| Document | File |
|---|---|
| Korean PDF | [OCRT_Manual_KO.pdf](OCRT_Manual_KO.pdf) |
| English PDF | [OCRT_Manual_EN.pdf](OCRT_Manual_EN.pdf) |
| Sources and build instructions | [latex/README.md](latex/README.md) |
| Korean entry point | [latex/ocrt_manual_ko.tex](latex/ocrt_manual_ko.tex) |
| English entry point | [latex/ocrt_manual_en.tex](latex/ocrt_manual_en.tex) |
| Executable examples | [examples/](examples/) |

## Updating

Edit the corresponding chapters in `latex/ko/` and `latex/en/`, then run this command from the repository root. The build reads the LaTeX sources and stored figure data; it does not convert Markdown.

```bash
python docs/manual/latex/build.py
```

On Linux, use `python3` if `python` is unavailable. See the [build instructions](latex/README.md) for XeLaTeX, fonts and Windows/WSL support. Both PDFs are written to this directory. Keep option names, units and defaults consistent between editions, and review the rendered pages after updates.

## Contents

The opening **Quick start** covers minimal setup, one calculation, a 12-row Rayleigh LUT, and expected success checks, with Linux/WSL and Windows PowerShell commands.

1. Concepts and viewing geometry, including the RAA convention
2. Installation and first run
3. Complete C option reference
4. Input files and batch CSV contracts
5. Output formats and interpretation
6. Verification, troubleshooting and provenance

Detailed applications and the introduction and usage of code calling OCRT are collected in the appendices. Short option syntax and file-format samples remain with their definitions in the main chapters.

| Appendix | Contents |
|---|---|
| A | Worked examples and LUT generation |
| B | Python and batch data generation |
| C | 10 core worked examples |
| D | C LUT example utility |
| E | Environment-variable source index |

The manual includes 12 vector figures for observation geometry, key options and output normalization, including an azimuth plot calculated with the C solver. The editable figures and plotted data are under `latex/figures/`. The 10 core examples and their runner are documented in Appendix C and [examples/README.md](examples/README.md).

The current Appendix C and default runner catalogue contain **10 examples and 13 C invocations**. The [verification summary](examples/verification_core_examples.json) selects unchanged commands from execution records dated 2026-09-11; it does not claim new solver runs on 2026-09-12. Execution and output checks are distinct from establishing physical accuracy for all conditions.

The [review list](review/README.md) describes all original 40 examples. The other 30 are deferred and excluded from the current manual and default runner catalogue; add them after review. Original LaTeX, the full catalogue and the historical 49-invocation record are preserved under `review/recipes_2026-09-12/`. The `review/` directory is not a PDF build input.

**OCRT RAA=180° is the direct-sunglint branch.** Read Chapter 1 to confirm the angle conversions and the sign of Stokes U.

The source baseline is `54861c68e35ec5aaa3938fad53e93dc2dfa7eb1d`, C documentation release v1.11.1. The earlier numbered Markdown files and `03_environment_index.md` are retained as the initial draft record. Future updates belong in the LaTeX sources and rebuilt PDFs.

---

# OCRT 사용자 매뉴얼

**수정 원본은 LaTeX이다.**

| 문서 | 파일 |
|---|---|
| 한글판 PDF | [OCRT_Manual_KO.pdf](OCRT_Manual_KO.pdf) |
| 영문판 PDF | [OCRT_Manual_EN.pdf](OCRT_Manual_EN.pdf) |
| LaTeX 원본과 빌드 안내 | [latex/README.md](latex/README.md) |
| 한글판 진입 파일 | [latex/ocrt_manual_ko.tex](latex/ocrt_manual_ko.tex) |
| 영문판 진입 파일 | [latex/ocrt_manual_en.tex](latex/ocrt_manual_en.tex) |
| 실행 가능한 예시 | [examples/](examples/) |

## 업데이트

`latex/ko/`와 `latex/en/`에서 해당 장의 원본을 수정한다. 저장소 최상위에서 다음 명령을 실행하면 두 PDF를 각각 생성한다. 빌드는 LaTeX 원본과 저장된 그림 자료를 읽으며 Markdown 변환을 수행하지 않는다.

```bash
python docs/manual/latex/build.py
```

Linux에서 `python`이 없으면 `python3`를 사용한다. XeLaTeX와 글꼴 설치, Windows/WSL 실행 방법은 [빌드 안내](latex/README.md)를 따른다. 출력은 이 디렉토리의 한글·영문 PDF다. 옵션·단위·기본값을 바꾸면 두 언어판을 함께 수정하고 페이지 배치도 확인한다.

## 구성

매뉴얼 앞부분의 **빠른 시작(Quick start)**은 최소 준비, 단일 계산, 12행 Rayleigh LUT 생성과 정상 실행 확인을 안내한다. Linux/WSL과 Windows PowerShell 명령을 함께 제공한다.

1. 기본 개념과 관측 기하
2. 설치와 첫 실행
3. C 실행 프로그램 전체 옵션
4. 입력 파일과 배치 CSV 계약
5. 출력 형식과 결과 해석
6. 검증·문제 해결·문서 근거

상세 활용 예시와 OCRT를 호출하는 코드의 소개·사용법은 부록에 모았다. 본문의 짧은 옵션 문법과 파일 형식 샘플은 해당 기능 설명에 포함된다.

| 부록 | 내용 |
|---|---|
| A | 상세 실행 예시와 LUT 제작 |
| B | Python 사용과 배치 자료 생성 |
| C | 핵심 사용 예시 10개 |
| D | C LUT 예시 도구 |
| E | 환경변수 소스 색인 |

관측기하와 주요 옵션, 출력 높이·정규화를 설명하는 벡터 그림 12개를 포함한다. SZA/VZA와 RAA, 해면 경계, 직달 글린트, 압력과 기체, AOD 기준 파장, 수중 모델, 출력 격자와 내부 계산을 그림으로 구분한다. 실제 C 계산에서 얻은 RAA 곡선과 재생성 도구도 함께 제공한다.

현재 부록 C와 기본 실행 목록은 R01, R02, R07, R21, R26, R28, R36, R38, R39, R40의 **10개 예시·13개 C 호출**이다. [검증 요약](examples/verification_core_examples.json)은 2026-09-11의 실행 기록에서 변경 없는 13개 명령을 선별한 것으로, 2026-09-12에 새로 실행했다는 뜻은 아니다. 실행 및 출력 검사는 모든 조건의 물리적 정확도 검증과 구분한다.

원래 40개 각각의 설명과 유지·보류 상태는 [검토 목록](review/README.md)에 있다. 나머지 30개는 현재 매뉴얼과 기본 실행 목록에서 제외했으며 검토 후 추가한다. 원본 LaTeX, 전체 카탈로그와 49건의 과거 검증 기록은 `review/recipes_2026-09-12/`에 보존한다. `review/`는 PDF 빌드 입력이 아니다.

**OCRT RAA=180°는 선글린트 가지다.** 각도 변환과 Stokes U의 부호는 1장을 읽고 확인한다.

기준 소스는 `54861c68e35ec5aaa3938fad53e93dc2dfa7eb1d`, C 문서 릴리스는 v1.11.1이다. 기존 `01_*.md`부터 `08_*.md`, `03_environment_index.md`는 2026-09-11의 최초 Markdown 초안 기록이다. 앞으로는 LaTeX 원본과 PDF를 갱신한다.
