# OCRT documentation guide

Installation and execution of OCRT start with the [user manual](manual/README.md). A [Korean PDF](manual/OCRT_Manual_KO.pdf) and an [English PDF](manual/OCRT_Manual_EN.pdf) are provided; edit the [LaTeX sources](manual/latex/README.md) to regenerate the PDFs. The documentation directory of this repository is `docs/`.

| Directory | Contents | When to read |
|---|---|---|
| [manual/](manual/README.md) | Korean and English PDFs and the LaTeX sources. Basic concepts, angle conventions, installation, the complete option reference, examples, input/output formats, Python batch processing | When installing for the first time, or when producing calculations and LUTs |
| [technical/](technical/) | Technical records on the physical model, spectral data, algorithm changes and performance improvements | When the implementation rationale or the background of a model change is needed |
| [validation/](validation/) | Validation records by feature: comparison experiments, error analyses, IPSS and others | When checking the validation conditions and limits of a specific result |
| [reports/](reports/) | Result reports and figures for individual work items | When looking for the results of an individual experiment or development task |

`technical/`, `validation/` and `reports/` record the state at the time they were written. For the current usage, consult the manual together with the source of the corresponding version. Do not assume that the options and defaults in older documents match the current code.

`paper/`, which may exist in a local working folder, is a private area holding the manuscript and is not part of the public repository distribution.

## Related locations

- [Repository guide](../README.md)
- [Technical documents inside the C code](../code/OCRT_C/docs/)
- [Python code and documents](../code/OCRT_Python/)
- [Development status record](../OCRT_MIGRATION_STATUS_2026-08-23.md)
- [Work hand-over and data inventory](../MIGRATION_README_2026-09-05_KO.md)

When adding a new document, put execution procedures for users in `manual/`, implementation descriptions in `technical/`, comparison conditions and numerical results in `validation/`, and individual work reports in `reports/`.

---

# OCRT 문서 안내

OCRT의 설치와 실행은 [사용자 매뉴얼](manual/README.md)에서 시작한다. [한글판 PDF](manual/OCRT_Manual_KO.pdf)와 [영문판 PDF](manual/OCRT_Manual_EN.pdf)를 제공하며, [LaTeX 원본](manual/latex/README.md)을 수정해 PDF를 재생성한다. 이 저장소의 문서 디렉토리 이름은 `docs/`이다.

| 디렉토리 | 내용 | 읽는 시점 |
|---|---|---|
| [manual/](manual/README.md) | 한글·영문 PDF와 LaTeX 원본. 기본 개념, 각도 규약, 설치, 전체 옵션, 예시, 입출력 형식, Python 배치 | 처음 설치하거나 계산·LUT를 만들 때 |
| [technical/](technical/) | 물리 모델, 스펙트럼 자료, 알고리즘 변경, 성능 개선에 관한 기술 기록 | 구현 근거나 모델 변경의 배경이 필요할 때 |
| [validation/](validation/) | 비교 실험, 오류 분석, IPSS 등 기능별 검증 기록 | 특정 결과의 검증 조건과 한계를 확인할 때 |
| [reports/](reports/) | 작업별 결과 보고서와 그림 | 개별 실험·개발 작업의 결과를 찾을 때 |

`technical/`, `validation/`, `reports/`는 작성 당시 상태를 기록한 자료다. 현재 실행법은 매뉴얼과 해당 버전 소스를 함께 확인한다. 과거 문서의 옵션·기본값이 현재 코드와 같다고 가정하지 않는다.

로컬 작업 폴더에 있을 수 있는 `paper/`는 논문 원고를 보관하는 비공개 영역이며, 공개 저장소의 배포 대상에 포함되지 않는다.

## 관련 위치

- [저장소 전체 안내](../README.md)
- [C 코드 내부 기술 문서](../code/OCRT_C/docs/)
- [Python 코드와 문서](../code/OCRT_Python/)
- [개발 상태 기록](../OCRT_MIGRATION_STATUS_2026-08-23.md)
- [작업 인계 및 자료 목록](../MIGRATION_README_2026-09-05_KO.md)

새 문서를 추가할 때는 사용자에게 필요한 실행 절차는 `manual/`, 구현 설명은 `technical/`, 비교 조건·수치 결과는 `validation/`, 개별 작업 보고는 `reports/`에 둔다.
