# OCRT migration README — v1.11.1 (2026-09-05)

This document lets a new work session continue from this repository alone. Order: bootstrap →
current state → inventory (what is in git and what stays local) → build and gates → GitHub
update procedure → rules → open items. The earlier package notes (2026-08-19) are accumulated in
`OCRT_MIGRATION_STATUS_2026-08-23.md` and are not repeated here.

`<LOCAL_ROOT>` denotes the local working folder into which this repository is checked out (the
repository root).

---

## 1. Bootstrap of a new session

1. Clone the repository (code, documents and validation records, about 110 MB; no Mie tables).

       git clone https://github.com/brtnt/OCRT.git ocrt
       cd ocrt/code/OCRT_C && make        # a sandbox with AVX-512 can use the Makefile default
       bash scripts/run_ipss_gates.sh     # runs without Mie tables; expect ALL PASS

2. Connect the local folder `<LOCAL_ROOT>`. It is the authoritative copy and the repository root
   (`.git` lives there). The local VM has git 2.34 and reaches github.com. Delete permission must
   be granted again in every session (git needs it for lock files).

3. Mie tables come from the GitHub Release `data-v1` (four zip files, 1.3 GB → 181 tables, 4.3 GB):
   `bash scripts/fetch_data.sh` (Windows: `scripts\fetch_data.ps1`). If release downloads are
   blocked in the sandbox, stage only the tables that are needed from the local folder — the IPSS
   bench and the SZA sweep need only `code/OCRT_C/inputs/M80C.mie` (24 MB); the coupled-ocean
   harness needs `C50.mie`, `tsm_ahn/Red_clay_AHN.mie`, `tsm_ahn/Brown_earth_AHN.mie` and
   `water_iop/Detritus_Stramski2001.mie`. The zip originals are kept locally in
   `package/release_data_v1/`. After installation check:
   `cd code/OCRT_C/inputs && sha256sum -c ../../../scripts/MIE_SHA256SUMS_data-v1.txt`.

4. The Windows executable is built locally by the user: `code\OCRT_C\build_win.bat`
   (`REBUILD_WIN_v1.11.1_KO.txt`). The local CPU (Core Ultra 9 275HX) has no AVX-512, so
   `-march=x86-64-v3` is required; a `cascadelake` build dies with "Illegal instruction".

## 2. Current state (v1.11.1, 2026-09-05)

- Spherical correction = IPSS only. Zhai & Hu (2022) Eq. (7): κ = I₁,ss / I₁,pp multiplies the
  final (I, Q, U). The old average-secant Chapman PSSA is deleted from the source tree; `--pssa` is
  a plain on/off switch (`--pssa-mode` exits with code 2).
- Viewing-angle anchor = surface pixel (satellite L1B senz). A TOA anchor leaves a false
  offset of −2.9 % (VZA 55) to −11 % (VZA 70) even with the sun at zenith. Gate G-I3 checks the
  surface-collapse limit κ → 1.
- κ is phase-function weighted (v1.11.1). In a mixed atmosphere (aerosol 2 km / Rayleigh ~8 km)
  the phase functions do not cancel. Layer source = β_j [w_R,j P_R(Θ) + w_A,j P_A(Θ)]; P_A is taken
  from the FR631 value table. Rayleigh-only results are bit-identical to v1.11. Gates G-I4a/b
  (independent scipy quadrature, |Δκ| ≤ 7.5e-6).
- The non-PSSA path is unchanged: smoke SHA `4ced4f18d50e4ac9`, ten-case bit battery.
- Independent cross-check: with the same Rayleigh profile, κ agrees with the IPSS authors'
  reference implementation within 0.1 %; the layer-count extrapolation of that code gives
  −0.239 % versus OCRT −0.240 % at VZA 55.
- The residual of −0.2 % at VZA 55 and low SZA is the view-path curvature term
  −(h̄/R_e) tan²θ_v; it does not depend on SZA and vanishes at nadir. The legacy scheme omitted it.
  The correction matters for SZA ≥ 75°.
- Details: status document §14–§14.1, the IPSS record (§1–§11.6), and figures figI1–figI5 in
  `docs/reports/`.
- Note: the `V2_VERSION` string in `src/main.c` follows the v1.2 lineage and differs from the
  document version (v1.11.1).

## 3. Inventory — in git / local only

In git (public): `code/OCRT_C` (src, tests including one fixture Mie table, scripts, tools, docs,
validation, patches, manifests), `code/OCRT_Python` (code, documents, tests, small data: AFGL,
cross sections, IOP tables), `scripts/` (data installation and manifests), `code/campaign_runner`,
`code/run_tools`, `code/validation_05` (harness, bit baselines, tools), `docs/technical`,
`docs/validation`, `docs/reports`, the status document, this README, `LICENSE`, `.gitignore`,
`.gitattributes` (`* -text`, byte-exact storage).

Local only (`.gitignore`):

| Item | Location | Size | Note |
|---|---|---|---|
| Mie tables `*.mie` (181, same files in both trees) | `code/OCRT_C/inputs/`, `code/OCRT_Python/data/` | 4.3 + 4.3 GB | GitHub Release `data-v1` (four zips, 1.3 GB); local originals in `package/release_data_v1/`; received originals `package/PART01–PART10*.zip` |
| Reference-run raw text | `code/validation_05/runs/` | 158 MB | consumed by re-analysis only |
| Delivery archives | `package/` | 2.3 GB | per-session tar.gz / zip |
| Campaign outputs | `task/` | 1.9 GB | sensitivity analyses etc. |
| Unpublished manuscript | `docs/paper/` | 69 MB | LaTeX, PDF, figures |
| Large data of technical reports | `docs/technical/**/*.gz` | 64 MB | polarization reference CSV |
| Build products | `build/`, `*.exe` | — | rebuild |

Third-party codes used for validation are not part of this package and are not stored here.

## 4. Build and gates (verbatim)

    # Linux / sandbox
    cd code/OCRT_C && make
    # host without AVX-512 (including the local VM)
    gcc -std=c11 -O3 -march=x86-64-v3 -ffp-contract=fast -fassociative-math -fno-signed-zeros \
        -fno-trapping-math -DOCRT_FAST_KERNELS -fopenmp -Isrc $(find src -name '*.c') -o build/ocrt -lm
    # gates (G-I0 … G-I4 and regression §2–§5)
    bash scripts/run_ipss_gates.sh
    # manifest (validation/ipss_2026-09-05, build/ and _RELEASE/ are excluded by convention)
    find . -type f ! -path "./SHA256SUMS.txt" ! -path "./build/*" ! -path "./_RELEASE/*" \
         ! -path "./validation/ipss_2026-09-05/*" -print0 | sort -z | xargs -0 sha256sum \
         | sed 's#  \./#  #' > SHA256SUMS.txt

Environment variables: `OCRT_ADVANCED=1` (advanced options), `OCRT_DEBUG=1`, `OCRT_IPSS_DUMP=1`
(κ diagnostics), `OCRT_IPSS_PROFILE_DUMP=<file>` (dump of the κ input profile and P11 table),
`OCRT_IPSS_NQUAD` (quadrature convergence gate); timing runs use `OMP_NUM_THREADS=1`.

## 5. GitHub update procedure (only when the user asks for a sync — not automatic)

1. Confirm ALL PASS on the gates, regenerate `SHA256SUMS.txt`, add a section to the status document.
2. Commit in the local VM (repository root = connected folder):

       cd ~/mnt/OCRT && git add -A && git status --short | head
       git commit -m "<summary>"

3. Push. The token stays outside the repository: one line in `.github_token` at the folder root
   (ignored by git). Use a fine-grained token limited to this repository (Contents: read/write),
   90-day expiry; revoke it on GitHub if it may have leaked.

       git -c credential.helper='!f() { echo "username=brtnt"; echo "password=$(cat .github_token)"; }; f' push origin main

   No token is stored in the remote URL (`origin` is `https://github.com/brtnt/OCRT.git`).
   Repository: https://github.com/brtnt/OCRT (public). First push 2026-09-05.
4. After the push, check the commit SHA and file count through the GitHub API and update §2.

## 6. Rules (fixed)

- Code changes only after explicit user approval. If a gate does not pass, revert.
- Phase representation policy: no higher moments and no increase of L; the standard is linear
  interpolation with a dense forward grid (FR631, truncation OFF is the scientific reference).
- The non-PSSA path stays bit-identical (proved by smoke SHA and bit battery).
- Consistency comparisons use 1:1 scatter plots with error bands; speed comparisons are valid only
  as "% of the reference" measured again on the same host.
- Underwater comparisons exclude wind speeds below 3 m/s (3 included).
- RAA convention: 180° = specular (glint), 0° = backscattering. RT wavelengths ≤ 2000 nm.
- Public repository: the manuscript (`docs/paper`), large data and tokens are never uploaded.
  License: academic and non-commercial (`LICENSE`).

## 7. Open items

1. Rebuild the Windows `.exe` from the v1.11.1 source (user).
2. Accuracy of the Eq. (7) approximation itself in aerosol-laden atmospheres (transfer of the
   single-scattering ratio to multiple scattering): published validation is Rayleigh-only.
3. Phase II — policy for applying Eq. (7) to the coupled atmosphere–ocean mode (`--water-model`);
   currently fail-loud.
4. Decision on re-running past `--pssa` campaigns (legacy + TOA anchor); slots with SZA ≥ 75° are
   affected.
5. Option re-classification (production / advanced / debug), and the open science gates of the
   status document §1–§13 (raw FR631 convergence freeze, coupled-LUT state recovery, Rrs-Q
   reference physics, Chl/detritus end-to-end validation).
6. Align the `V2_VERSION` string in `src/main.c` with the document version scheme.

---

# OCRT 마이그레이션 README — v1.11.1 (2026-09-05)

새 작업 세션이 이 저장소 하나로 작업을 이어가기 위한 문서다. 순서: 부트스트랩 → 현재 상태 →
자료 위치(무엇이 git 에 있고 무엇이 로컬에만 있는가) → 빌드·게이트 → GitHub 반영 절차 → 규칙 →
남은 과제. 이전 패키지 노트(2026-08-19)는 `OCRT_MIGRATION_STATUS_2026-08-23.md` 에 누적되어 있으므로
여기서 반복하지 않는다.

`<LOCAL_ROOT>` 는 이 저장소를 체크아웃한 로컬 작업 폴더(저장소 루트)를 뜻한다.

---

## 1. 새 세션 부트스트랩

1. 저장소를 받는다(코드·문서·검증 기록 약 110 MB, Mie 표 없음).

       git clone https://github.com/brtnt/OCRT.git ocrt
       cd ocrt/code/OCRT_C && make        # AVX-512 가 있는 샌드박스는 Makefile 기본값 사용 가능
       bash scripts/run_ipss_gates.sh     # Mie 표 없이 실행됨. ALL PASS 확인

2. 로컬 폴더 `<LOCAL_ROOT>` 를 연결한다. 이 폴더가 정본이며 저장소 루트다(`.git` 이 여기 있다).
   로컬 VM 에는 git 2.34 가 있고 github.com 에 닿는다. 삭제 권한은 세션마다 다시 받아야 한다
   (git 의 lock 파일 처리에 필요).

3. Mie 표는 GitHub Release `data-v1`(zip 4개, 1.3 GB → 181 표 4.3 GB)에서 받는다:
   `bash scripts/fetch_data.sh` (Windows: `scripts\fetch_data.ps1`). 샌드박스에서 릴리스 다운로드가
   막혀 있으면 로컬 폴더에서 필요한 표만 스테이징한다 — IPSS 벤치·SZA 스윕은
   `code/OCRT_C/inputs/M80C.mie` 하나(24 MB), 결합 해양 하니스는 `C50.mie`,
   `tsm_ahn/Red_clay_AHN.mie`, `tsm_ahn/Brown_earth_AHN.mie`, `water_iop/Detritus_Stramski2001.mie`.
   로컬에는 zip 원본이 `package/release_data_v1/` 에 있다. 설치 후 검증:
   `cd code/OCRT_C/inputs && sha256sum -c ../../../scripts/MIE_SHA256SUMS_data-v1.txt`.

4. Windows 실행파일은 사용자가 로컬에서 빌드한다: `code\OCRT_C\build_win.bat`
   (`REBUILD_WIN_v1.11.1_KO.txt`). 로컬 CPU(Core Ultra 9 275HX)는 AVX-512 가 없으므로
   `-march=x86-64-v3` 가 필수다. `cascadelake` 빌드는 Illegal instruction 으로 죽는다.

## 2. 현재 상태 (v1.11.1, 2026-09-05)

- 구면 보정 = IPSS 단일. Zhai & Hu (2022) Eq.(7) κ = I₁,ss / I₁,pp 를 최종 (I, Q, U) 에 곱한다.
  구 평균할선 Chapman PSSA 는 소스 트리에서 삭제됐고 `--pssa` 는 on/off 스위치뿐이다
  (`--pssa-mode` 는 exit 2).
- 관측각 기준 = 지표 화소(위성 L1B senz). TOA 기준은 태양이 천정에 있어도 −2.9 %(VZA 55)~
  −11 %(VZA 70)의 가짜 편차를 남긴다. 게이트 G-I3 가 지표 붕괴 극한 κ→1 을 검사한다.
- κ 는 위상함수 가중(v1.11.1). 혼합 대기(에어로졸 2 km / Rayleigh ~8 km)에서는 위상함수가
  상쇄되지 않는다. 층 원천 = β_j[w_R,j P_R(Θ) + w_A,j P_A(Θ)], P_A 는 FR631 값표를 쓴다.
  Rayleigh 단독 결과는 v1.11 과 비트 동일. 게이트 G-I4a/b(독립 scipy 적분과 |Δκ| ≤ 7.5e-6).
- 비-PSSA 경로 불변: 스모크 SHA `4ced4f18d50e4ac9`, 10케이스 비트 배터리.
- 독립 교차검증: 같은 Rayleigh 프로파일에서 κ 가 IPSS 저자의 참조 구현과 0.1 % 이내로 일치하고,
  그 코드의 층수 외삽값은 VZA 55 에서 −0.239 %(OCRT −0.240 %)다.
- VZA 55·저 SZA 에서 남는 −0.2 % 는 시선경로 곡률 항 −(h̄/R_e)tan²θ_v 이며, SZA 에 무관하고
  천저에서 0 이다. 구 방식이 빠뜨렸던 항이다. 구면 보정이 실제로 필요한 구간은 SZA 75° 이상이다.
- 상세: 상태 문서 §14·§14.1, IPSS 기록(§1–§11.6), `docs/reports/` 의 그림 figI1–figI5.
- 참고: `src/main.c` 의 `V2_VERSION` 문자열은 v1.2 계보 표기이고 문서 버전(v1.11.1)과 체계가 다르다.

## 3. 자료 위치 — git 에 있는 것 / 로컬에만 있는 것

git 에 있음(공개): `code/OCRT_C`(src·tests(픽스처 Mie 1개 포함)·scripts·tools·docs·validation·
patches·매니페스트), `code/OCRT_Python`(코드·문서·테스트·소형 data: AFGL·단면적·IOP 표),
`scripts/`(자료 설치·매니페스트), `code/campaign_runner`, `code/run_tools`,
`code/validation_05`(harness·bit_baselines·tools), `docs/technical`, `docs/validation`, `docs/reports`,
상태 문서, 이 README, `LICENSE`, `.gitignore`, `.gitattributes`(`* -text`, 바이트 보존).

로컬에만 있음(`.gitignore`):

| 항목 | 위치 | 크기 | 비고 |
|---|---|---|---|
| Mie 표 `*.mie` (181, 두 트리에 같은 파일) | `code/OCRT_C/inputs/`, `code/OCRT_Python/data/` | 4.3 + 4.3 GB | GitHub Release `data-v1`(zip 4개, 1.3 GB); 로컬 원본 `package/release_data_v1/`; 수령 원본 `package/PART01–PART10*.zip` |
| 참조 런 원문 | `code/validation_05/runs/` | 158 MB | 재분석용 소비 파일 |
| 배포 묶음 | `package/` | 2.3 GB | 세션별 tar.gz·zip |
| 캠페인 산출물 | `task/` | 1.9 GB | 민감도 분석 등 |
| 투고 전 원고 | `docs/paper/` | 69 MB | LaTeX·PDF·그림 |
| 기술 보고서 대형 자료 | `docs/technical/**/*.gz` | 64 MB | 편광 참조 CSV |
| 빌드 산출물 | `build/`, `*.exe` | — | 재빌드 |

검증에 쓴 제3자 코드는 이 패키지의 일부가 아니며 여기에 저장하지 않는다.

## 4. 빌드·게이트 (verbatim)

    # Linux / 샌드박스
    cd code/OCRT_C && make
    # AVX-512 없는 호스트(로컬 VM 포함)
    gcc -std=c11 -O3 -march=x86-64-v3 -ffp-contract=fast -fassociative-math -fno-signed-zeros \
        -fno-trapping-math -DOCRT_FAST_KERNELS -fopenmp -Isrc $(find src -name '*.c') -o build/ocrt -lm
    # 게이트 (G-I0~G-I4 + 회귀 §2–§5)
    bash scripts/run_ipss_gates.sh
    # 매니페스트 재생성 (validation/ipss_2026-09-05, build/, _RELEASE/ 는 제외 — 기존 규약)
    find . -type f ! -path "./SHA256SUMS.txt" ! -path "./build/*" ! -path "./_RELEASE/*" \
         ! -path "./validation/ipss_2026-09-05/*" -print0 | sort -z | xargs -0 sha256sum \
         | sed 's#  \./#  #' > SHA256SUMS.txt

환경변수: `OCRT_ADVANCED=1`(고급 옵션), `OCRT_DEBUG=1`, `OCRT_IPSS_DUMP=1`(κ 진단),
`OCRT_IPSS_PROFILE_DUMP=<file>`(κ 입력 프로파일·P11 표 덤프), `OCRT_IPSS_NQUAD`(구적 수렴 게이트),
속도 측정은 `OMP_NUM_THREADS=1`.

## 5. GitHub 반영 절차 (사용자가 동기화를 요청할 때만 — 자동·상시 동기화 아님)

1. 게이트 ALL PASS 확인, `SHA256SUMS.txt` 재생성, 상태 문서에 절 추가.
2. 로컬 VM 에서 커밋(저장소 루트 = 연결 폴더):

       cd ~/mnt/OCRT && git add -A && git status --short | head
       git commit -m "<요약>"

3. 푸시. 토큰은 저장소 밖에 둔다 — 폴더 루트의 `.github_token`(gitignore 됨)에 한 줄. 토큰은 이 repo
   한정 fine-grained(Contents: Read/Write), 만료 90일 권장, 유출 의심 시 GitHub 에서 즉시 폐기.

       git -c credential.helper='!f() { echo "username=brtnt"; echo "password=$(cat .github_token)"; }; f' push origin main

   원격 URL 에 토큰을 저장하지 않는다(`origin` 은 `https://github.com/brtnt/OCRT.git`).
   저장소: https://github.com/brtnt/OCRT (공개). 첫 푸시 2026-09-05.
4. 푸시 뒤 GitHub API 로 커밋 SHA 와 파일 수를 확인하고 §2 를 갱신한다.

## 6. 규칙 (불변)

- 코드 수정은 사용자 명시 승인 후에만. 게이트 미완이면 원복.
- 위상 표현 정책: 고차 모먼트·L 상향 금지, 표준 = 선형 보간 + 전방 조밀(FR631, 절단 OFF 가 과학 기준).
- 비-PSSA 경로는 비트 불변이 원칙(스모크 SHA·비트 배터리로 증명).
- 정합 비교는 1:1 산포도(오차 밴드) 표준, 속도 비교는 동일 호스트 재실측 "기준 대비 %" 만 유효.
- 수중 비교에서 풍속 3 m/s 미만 제외(3 포함).
- RAA 규약: 180° = 경면(글린트), 0° = 후방산란. RT 파장 ≤ 2000 nm.
- 공개 저장소: 원고(`docs/paper`)·대용량 자료·토큰은 올리지 않는다. 라이선스는 학술·비상업(`LICENSE`).

## 7. 남은 과제

1. Windows `.exe` 재빌드(v1.11.1 소스 반영) — 사용자 실행.
2. 에어로졸 대기에서 Eq.(7) 근사 자체(단일산란 비의 다중산란 전이)의 정확도 — 발표된 검증은 Rayleigh 뿐.
3. Phase II — 결합 대기–해양(`--water-model`)에 Eq.(7) 적용 정책(현재 fail-loud).
4. 과거 `--pssa` 캠페인 산출물(legacy + TOA 앵커) 재실행 판정 — SZA 75° 이상 슬롯이 영향권.
5. 옵션 재분류(production/advanced/debug), 상태 문서 §1–§13 의 미결 과학 게이트(raw FR631 convergence
   freeze, coupled-LUT state 복구, Rrs-Q 참조물리 정합, Chl/detritus end-to-end 검증).
6. `src/main.c` `V2_VERSION` 문자열 정리(문서 버전 체계와 일치).
