# OCRT 레퍼런스 논문 — 세션 마이그레이션 문서

**작성 시점**: 2026-07-20
**목적**: 이미지 100개로 세션이 꽉 차서, 새 세션으로 작업을 이어받기 위한 인수인계 문서.
**사용자(Jae)**: 한국해양과학기술원(KIOST) 해양위성센터(KOSC) 소속 해색/복사전달 연구자. 응답은 반말 논문체(~한다/~이다), 팩트 우선, 틀리면 공감 없이 지적, 모르면 모른다고 명확히.

---

## 1. 프로젝트 개요

**OCRT 레퍼런스 논문**을 작성 중이다. OCRT(Ocean Color Radiative Transfer)는 KOSC/KIOST에서 개발한 **결합 해양-대기 벡터 복사전달 모델**(C11 + Python 포팅)이다. 이 논문은 OCRT를 인용 가능한 기준 논문으로 서술하는 model-description 논문이다(Chami 2015의 OSOAA 논문, He 2010의 PCOART 논문과 같은 성격).

**논문 목적**: 물리 모델·수치 방법 명세, 독립 코드 대비 검증, 계산 성능 특성화. 새 물리를 제안하는 게 아니라 "검증되고 재현 가능한, GOCI-III 대기보정 알고리즘 개발을 위한 순방향 모델 도구"로 위치시킨다.

**동기(motivation)**: GOCI-III(정지궤도 해색 위성, 편광 관측) 대기보정 알고리즘 개발. 후속 논문에서 단일 관측각 편광이 흡수성 에어로솔 대기보정을 개선하는지 다룰 예정(이건 별도 논문, 이 논문 아님).

---

## 2. 산출물 파일 (사용자가 직접 옮길 것)

`/mnt/user-data/outputs/` 에 있다:
- **OCRT_reference_paper_draft.docx** — 최종 Word 원고 (약 949 KB, 22페이지, 그림 5개 임베드)
- **OCRT_reference_paper_draft.md** — 같은 내용, LaTeX 수식 소스 보존 (저널 변환용)

**중요**: docx/md는 사용자가 옮긴다고 했다. 새 세션에서 이 두 파일을 다시 만들려면 아래 빌드 스크립트를 실행해야 한다.

---

## 3. 빌드 시스템 (새 세션에서 논문을 재생성하는 방법)

### 3.1 논문 빌드 스크립트
- **위치**: `/home/user/build_ocrt_paper.js` (약 94 KB, Node docx-js v9.6.1)
- **실행**: `NODE_PATH=/home/user/.npm-global/lib/node_modules node build_ocrt_paper.js`
- 출력: `/mnt/user-data/outputs/OCRT_reference_paper_draft.docx`
- 헬퍼 함수: `para/seg/eq/h/bullet/simpleTable/fillNote/figure/figcaption`
- 폰트: FONT='Calibri', SERIF='Cambria'. 수식은 eq()로 중앙정렬 이탤릭 유니코드(OMML 아님).
- `[TO COMPLETE]` = 주황색 fillNote 박스.
- 그림: `figure(파일명, 종횡비, {widthIn})` + `figcaption(라벨, 텍스트)`. FIG_DIR='/home/user/paper_figs'.

### 3.2 빌드 검증 (PDF 렌더링)
```
rm -f /tmp/render/*.jpg /tmp/render/*.pdf
python /mnt/skills/public/docx/scripts/office/soffice.py --headless --convert-to pdf /mnt/user-data/outputs/OCRT_reference_paper_draft.docx --outdir /tmp/render
cd /tmp/render && pdftoppm -jpeg -r 72 OCRT_reference_paper_draft.pdf page
```
그다음 `view /tmp/render/page-NN.jpg`로 확인.

### 3.3 반복되는 빌드 버그 (주의)
1. **단일 문자열 문단**: `P(para('...'))`는 끝이 `');`여야 한다. `')]));`로 끝나면 문법 오류(배열 닫기 대괄호는 `P(para([seg(...)]))`에만 붙는다). 이 버그가 여러 번 발생했다.
2. **그림 ENOENT**: 그림은 빌드 전 반드시 `/home/user/paper_figs/`에 cp되어 있어야 한다. 없으면 빌드가 조용히 실패하고 이전 PDF가 렌더링된다.
3. **이중 이스케이프**: str_replace로 `\u2013`(en-dash) 등을 넣을 때 `\\u2013`(이중 백슬래시)로 들어가면 리터럴로 출력된다. 파이썬 문자열 치환으로 파일을 편집하는 게 안전하다(파이썬이 이스케이프를 투명 처리). str_replace는 view에 보이는 실제 문자를 그대로 써야 한다.

---

## 4. 그림 생성 시스템

### 4.1 그림 생성 스크립트
- **위치**: `/home/user/compare_data/make_figs.py` (약 17 KB)
- **실행**: `cd /home/user/compare_data && python3 make_figs.py`
- 출력: `figs_final/` → 반드시 `/home/user/paper_figs/`로 cp해야 빌드가 인식.
- **통일 스타일**: DejaVu Serif, 큰 폰트(축라벨 16, 틱 14, 제목 17, 범례 14 — 사용자가 캡션 수준 크기 요청함), 마커 크기 큼.
- **통계 박스**: `stats_box()`(MAPE/RMSE/bias)와 `stats_box_ratio()`(Q/U용 RMS/I·max/I)가 각 산포도 패널 안에 통계를 넣는다. 모든 산포도에 통계 박스 필수(사용자 요청).

### 4.2 현재 논문의 그림 5개 (순서 중요)
논리 순서 = Rayleigh 대기 → shadowing → 에어로솔 → 순수해수 → 구성모델:
- **Figure 1**: `fig_ray_scatter_OSOAA.png` — 전격자 Rayleigh 산포도 (3밴드×IQU, 2940점/패널, 풍속 색분류, 1:1선+±1%밴드, 통계박스)
- **Figure 2**: `fig_shadowing_p6.png` — shadowing 발견 (P6 이식 전 17.18% → 후 0.07%)
- **Figure 3**: `fig_aerosol.png` — 에어로솔 (M80, AOT865=0.1, TOA I/Q/U + rrs, case 색분류, 통계박스)
- **Figure 4**: `fig_pure_ocean_rrs.png` — 순수해수 (Rrs(0+) + rrs(0−), SZA 색분류, 통계박스) — 신규
- **Figure 5**: `fig_inwater.png` — 구성모델 (Chl/TSM rrs + CDOM TOA I, SZA 색분류, 통계박스)

**삭제된 그림**: Rayleigh 각도곡선(I/Q/U vs VZA)은 사용자 요청으로 제거했다. make_figs.py에 `stokes_vza_figure()` 함수는 남아있지만 호출 안 함.

### 4.3 데이터 소스 (그림 재생성에 필요)
- `/home/user/compare_data/` — master1(대기 Rayleigh IQU), csv/master2(P6 shadowing), csv/master3(aCDOM), csv/master4(Chl/TSM rrs), csv/master5(timing), summaryD_timing.csv
- `/home/user/compare_data/item6_aerosol_162.csv` — 에어로솔 (최신: 새 composite, AOT865=0.1)
- `/home/user/compare_data/pure_ocean_rrs.csv` — 순수해수 Rrs/rrs (신규, gate3 150셀)
- `/home/user/newdata/OCRT_Rrs_rrs_Aerosol_IQU_Comparison_20260720/` — 최신 첨부 원본(README_KO.md에 공식 수치)

---

## 5. OCRT 코드 사실 (검증됨, 소스에서 확인)

- **코드 패키지**: `/home/user/ocrt/OCRT_v1.2_SNF_PLUS_ACCURT_AHMAD2010_80_AEROSOL_MODELS_2026-07-19/ocrt/`
- **언어**: **순수 C (C11)**, C++ 아님. .c 파일만, gcc -std=c11, malloc/free만(class/template/namespace/new/delete 0개). 사용자가 초기에 "C++"라 했으나 착각으로 확정됨.
- **버전 번호 논문에 절대 넣지 않음**: v1.1과 v1.2 결과 동일. 논문 전체에서 버전 번호 제거 완료.
- **배포**: 오픈소스. GitHub: **https://github.com/brtnt/OCRT** (Code and data availability 절에 명시됨)
- **병렬화**: C = OpenMP 멀티코어(case 루프), Python = CuPy GPU. 둘 다 오픈소스.

### 5.1 물리 모델 출처 (논문에 인용된 것, 소스 검증 완료)
- 대기 연직: **U.S. Standard Atmosphere (1962)** 기본 + AFGL 6종(Anderson 1986)
- 레일리: **Bodhaine 1999** (τ_R), 공기 굴절률 **Peck & Reeder 1972**, King factor, depol 0.0279(Young 1980)
- 가스 흡광 (6종: H2O/O3/NO2/O2/CO2/CH4): **HITRAN2020 (Gordon et al. 2022)** LBL + O3/NO2 연속흡수 **Serdyuchenko 2014**(O3)·**Bogumil 2003**(NO2). Beer-Lambert. AFGL 프로파일 컬럼. 사용자 컬럼 오버라이드(예: 오존 DU). **HITRAN 버전은 2020 최신 확정, Bogumil은 2003 확정(주석의 2000은 오타)**.
- 에어로솔: 16종 OPAC(Shettle&Fenn 1979, Hess 1998) + 80종 Ahmad 2010/AccuRT. Mie(Bohren&Huffman 1983).
- 표면: Cox-Munk 1954, Sancer 1969 shadowing, 편광 프레넬. 해수 굴절률 기본 n=1.34(Quan&Fry 1995 옵션).
- 순수수: 흡광 **Pope&Fry 1997 + Kou 1993**, 산란 **Zhang et al. 2009**(T=20°C S=38.4), bb=0.5bw, depol δ_w=0.039 **Farinato&Rowell 1976**(Zhang 2009 채택값). NASA z09 파일 헤더로 확정.
- 구성모델: 식물플랑크톤 EAP(Lain&Bernard 2014/2023) + Huot 2008 후방산란, 광물 AHN, detritus Stramski 2001, CDOM 지수. CCRR Case-1 어댑터(Morel 1988/Morel&Maritorena 2001).
- 위상함수: Fournier-Forand 1994, 후방산란 폐형 **Mobley/Sundman/Boss 2002**.
- 전방피크 절단: δ-M(Wiscombe 1977), δ-fit(Hu 2000), log-linear(**Potter 1970** 정신), vector δ-M(Hu-Stamnes).
- PSSA: PCOART-SA(He et al. 2018) + Dahlback&Stamnes 1991 average-secant.

### 5.2 인용하지 않기로 한 코드 주석 레퍼런스 (사용 안 함/하위 참조)
List 1968(Bodhaine 하위), Sullivan 2006·Pegau/Gray/Zaneveld 1997(순수수 대체 데이터, 실제는 P&F/Kou 사용), Walter 2007(Cox-Munk 세부), Bricaud 1998(주석이 "not used"라 명시), Travis&Lacis 2002(Mishchenko 2002 공저자). 이들을 넣으면 "안 쓴 걸 인용"이 되므로 제외.

---

## 6. 검증 결과 (논문 5.1절)

- **검증 기준 코드**: (1) **Ahmad & Fraser 1982** (AF1982) — Rayleigh 대기 독립 검증(완전 독립 레거시 Fortran). (2) **OSOAA (Chami 2015)** — 에어로솔/구성모델/결합. 단, 사용자가 OSOAA에 shadowing 패치 등 수정함. 6S는 결합 해양 없어 제외.
- **핵심 발견 — shadowing**: cross-comparison에서 발견한 가장 중요한 결과. shadowing 없을 때 OSOAA가 최대 +17.2% 밝음(고천정각). 3중 확증(방위각평균 지배, 커널비율 1.0699≈Sancer 1/S 1.0700, 진단 빌드 4.81%→0.02%). OCRT=additive Sancer, AF1982=multiplicative Smith, public OSOAA=없음.
- **수치 결과** (발표 시점에 달라질 수 있음, 학회 초록에는 수치 뺌):
  - Rayleigh I: mean 0.030%/max 0.205%, Q 0.015/0.177, U 0.005/0.157
  - 에어로솔: TOA I MAPE 0.58%, Q RMS/I 0.52%, U RMS/I 0.13%, rrs 2.52%
  - 순수해수: Rrs(0+) MAPE 2.08%, rrs(0−) MAPE 2.59%
  - 구성모델: Chl 3.3%, TSM 2.3%, CDOM 0.2%
- **U 부호 규약**: U(OSOAA) = −U(OCRT), Stokes 기준면 정의 차이(수치 불일치 아님). 그림·본문에 명시.

---

## 7. 남은 작업 ([TO COMPLETE] 5개, 논문 완성용)

논문 build_ocrt_paper.js에 fillNote로 표시된 것:
1. **Table 5 속도 미측정 3열** (line 475): 멀티코어 C, 단일 Python, GPU Python. 현재 단일코어 C와 OSOAA만 실측(nadir ~3s 동급, 다중기하 OCRT 37.9s vs OSOAA 6.0s). **중요: 이건 future work 아니라 논문 완성용 TODO. 속도 최적화 진행 중이라 현재 격차는 최종 아님(single-solve multi-angle 작업 중).**
2. **7절 GOCI-III 응용 그림** (line 484): 흡수성 vs 비흡수성 에어로솔 분리도(DoLP 축 추가 효과). 예시 그림 1개.
3. **Zenodo DOI** (line 496): GitHub 릴리스에서 발급, 인용 영속성용.
4. **저자·자금** (line 499): 저자 목록, Acknowledgements(KIOST/KOSC, GOCI-III).
5. **저자 정보**: 제목 아래 [Author list — to be completed] 상태.

---

## 8. 문체 규칙 (전면 교정 완료됨, 유지할 것)

사용자가 자기 그룹 논문(KIOST, Optics Express/Remote Sensing) 문체로 전면 교정 요청함. 확인한 문체 특징:
- **em-dash(—) 삽입구 금지**: 별도 문장/쉼표/괄호/관계절로. 본문 삽입구 em-dash 0개로 교정 완료(표 빈 셀 `—`는 유지).
- **수사적 강조 배제**: "the single most consequential" → "the main", "not cosmetic; it originated" → "originated" 등.
- **담담한 서술, 짧은 문장, 능동태**: "This paper describes...", "We developed..." 형식.
- **약어는 정의 후 사용 OK**: "Ahmad and Fraser (1982), hereafter AF1982" 후 AF1982 사용(사용자 논문이 GW1994, SR660 식으로 씀). 단 초록에서는 정식 인용.
- 겸양 표현("not claimed as novel")은 유지(사용자 문체에 부합).

---

## 9. 학회 발표 (별도, 채팅에만 있음 — 파일로 안 만듦)

**학회**: "Advancing Climate Data Records — Creation and Application of Long-Term Ocean Color Products" 세션. 하지만 사용자가 **CDR 프레이밍 걷어내고** ocean color 알고리즘 개발을 위한 RT로 초점 맞추라 함(세션 주제 살짝 벗어남 감수).

**확정된 제목**:
> A coupled ocean–atmosphere vector radiative transfer model (OCRT) for ocean color algorithm development

**확정된 초록** (수치 제거, 발표 시점 변동 대비):
> Radiative transfer simulation of the coupled ocean–atmosphere system is a fundamental tool for developing ocean color algorithms such as atmospheric correction and in-water inversion. The accuracy of the radiative transfer code sets the accuracy limit of the algorithms built with it. This paper describes OCRT (Ocean Color Radiative Transfer), a vector radiative transfer model of the coupled ocean–atmosphere system developed at the Korea Institute of Ocean Science and Technology for ocean color algorithm development, with the Geostationary Ocean Color Imager III (GOCI-III) as the target sensor.
>
> OCRT computes the Stokes parameters I, Q, and U using the successive orders of scattering method over a coupled domain. It includes a molecular atmosphere, Mie aerosol models, a wind-roughened air–sea interface with polarized Fresnel reflection, and a multi-constituent ocean whose water-leaving reflectance is computed from the inherent optical properties of chlorophyll, suspended sediment, and colored dissolved organic matter.
>
> The code was validated against two independent vector radiative transfer codes, Ahmad and Fraser (1982) and OSOAA (Chami et al., 2015), for the Rayleigh atmosphere, the aerosol atmosphere, and the coupled ocean. The code is released as open source to support ocean color algorithm development.

**학회 초록 미확정**: 저자·소속, 단어 수 제한(현재 ~165단어), 구두/포스터 여부. 사용자가 원하면 별도 docx로 뽑을 수 있음.

---

## 10. 새 세션 시작 시 권장 순서

1. 이 문서(`/home/user/OCRT_paper_MIGRATION.md`)를 먼저 읽는다.
2. 사용자에게 docx/md 파일을 업로드받거나, 없으면 `build_ocrt_paper.js`가 유실됐는지 확인.
3. **만약 build_ocrt_paper.js도 유실됐다면**: 이 문서만으로는 94KB 스크립트 전체를 재현 불가. 사용자에게 docx를 받아 내용을 파악하거나, 처음부터 재작성해야 함. → **사용자가 build_ocrt_paper.js와 make_figs.py도 백업하는 것을 강력 권장** (docx/md만으로는 편집 재개가 어렵다).
4. 데이터 소스(compare_data, newdata, item6, ocrt 패키지)가 유실됐으면 그림 재생성 불가 → 사용자에게 재업로드 요청.
5. 남은 [TO COMPLETE] 5개 중 사용자가 자료를 주는 것부터 채운다.

---

## 11. 중요 경고

**docx/md만으로는 작업 재개가 제한적이다.** 논문은 `build_ocrt_paper.js`(스크립트)에서 생성되므로, 편집을 이어가려면 이 스크립트가 필요하다. docx를 직접 편집하면 스크립트와 동기화가 깨진다. 따라서:
- **사용자가 백업해야 할 파일** (docx/md 외 추가):
  - `/home/user/build_ocrt_paper.js` (논문 빌드 스크립트, 필수)
  - `/home/user/compare_data/make_figs.py` (그림 생성 스크립트)
  - `/home/user/paper_figs/*.png` (임베드된 그림 8개)
  - 데이터: `/home/user/compare_data/`, `/home/user/newdata/` (그림 재생성용)
- 이것들이 없으면 새 세션에서 논문 편집을 처음부터 다시 해야 한다.
