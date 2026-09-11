# OCRT manual: authoritative LaTeX source

The Korean and English PDFs are built directly from LaTeX. **Edit the `.tex` chapter files, then rebuild.** The build does not read or convert the earlier Markdown drafts.

| File or directory | Purpose |
|---|---|
| `ocrt_manual_ko.tex` | Korean edition entry point and front matter |
| `ocrt_manual_en.tex` | English edition entry point and front matter |
| `ko/`, `en/` | Matching, editable chapter sources and appendices |
| `ko/00_quick_start.tex`, `en/00_quick_start.tex` | Opening Quick start: minimal setup, one calculation and a 12-row LUT |
| `common/preamble.tex` | Shared typography, tables, code layout and PDF navigation |
| `figures/*.tex` | Editable TikZ/PGFPlots figures with shared Korean/English labels |
| `figures/data/` | Actual C output and metadata used by the azimuth plot |
| `ko/09_cookbook.tex`, `en/09_cookbook.tex` | Ten core worked examples in Appendix C |
| `build.py` | Three-pass XeLaTeX build, both editions by default |
| `build/` | Local intermediate files and logs; ignored by the repository |

## Document structure

The opening Quick start is followed by six main chapters: concepts and geometry,
installation, C options, input formats, output interpretation, and troubleshooting
and verification. Detailed applications and supporting code are in Appendices
A (worked examples and LUTs), B (Python and batch workflows), C (10 core examples),
D (C LUT utility), and E (environment-variable index).

The numbered source filenames are stable identifiers from the initial edition;
they are not the current printed chapter numbers. The entry-point `.tex` files
set the reading order. Use LaTeX labels and references for internal links.

## Build

Python 3 and XeLaTeX are required. On Windows, `build.py` uses native XeLaTeX when present, otherwise the Ubuntu WSL distribution. Both editions use Noto Serif/Sans CJK KR, TeX Gyre Pagella/Heros and DejaVu Sans Mono fonts. The CJK fonts are also needed for occasional original Korean identifiers in the English edition.

On Ubuntu/WSL, install the dependencies if they are missing:

```bash
sudo apt-get update
sudo apt-get install texlive-xetex texlive-latex-extra texlive-pictures texlive-lang-korean fonts-noto-cjk fonts-texgyre fonts-dejavu-core
```

From this directory:

```bash
python3 build.py
python3 build.py --language ko
python3 build.py --language en
```

On Windows, from the repository root:

```powershell
python docs\manual\latex\build.py
```

PDFs are written to `docs/manual/OCRT_Manual_KO.pdf` and `docs/manual/OCRT_Manual_EN.pdf`. Use `--out-dir PATH` to choose another destination. Build logs remain under `latex/build/ko/` and `latex/build/en/`. Check them for missing glyphs, undefined references and overfull boxes, and inspect rendered pages before publishing.

## Updating the manual

1. Edit the appropriate chapter in `ko/` and its matching English chapter in `en/`.
2. Keep command names, output column names, units and numerical examples consistent across editions.
3. If the source revision changes, update the baseline statements and verification scope in both entry points and chapters.
4. Run `build.py` and inspect both PDFs. Commit the sources and PDFs together using the repository's normal publication workflow.

Executable example resources are kept in `../examples/`. The source package includes them alongside this LaTeX tree. Those scripts run OCRT; they are not required to typeset the PDFs.

### Updating figures and recipes

Figures use native TikZ/PGFPlots. Edit the `.tex` file in `figures/`; shared labels use `\FigLang{Korean}{English}`. Diagrams distinguish conceptual sketches from computed data. For the numerical azimuth figure, rerun `../examples/generate_figure_data.py` in the full OCRT repository and replace the inspected CSV/metadata files in `figures/data/`. A normal LaTeX build only reads those stored data.

Appendix C is the authoritative explanation of the 10 retained examples (13 C invocations). Keep both language sources and `../examples/recipe_catalog.json` aligned when changing commands. `../examples/run_recipes.py` reads that catalogue to plan or run the finite examples; it never generates the manual text. Its validation checks do not replace scientific convergence studies. The original 40-example sources and the 30 deferred candidates are documented under `../review/`; those files are not included in either PDF build. Retain the original R identifiers when reviewing or adding an example.

## 한글 안내

수정 원본은 `ko/`와 `en/`의 LaTeX 파일이다. 해당 장의 두 언어판을 수정한 뒤 `python build.py`를 실행하면 한글·영문 PDF가 각각 생성된다. 이전 Markdown 초안은 빌드 입력으로 사용하지 않는다. 공통 글꼴·여백·표·코드 스타일은 `common/preamble.tex`에서 관리한다.
