# OCRT 마이그레이션 README — v1.11.1 (2026-09-05)

새 작업 세션이 이 저장소 하나로 작업을 재개하기 위한 문서다. 순서: 부트스트랩 → 현재 상태 →
인벤토리(무엇이 git 에 있고 무엇이 로컬에만 있는가) → 재빌드·게이트 → 세션 종료 시 GitHub 반영 절차 →
규칙 → 남은 과제. 이전 패키지(2026-08-19, `00_README_MIGRATION_KO.md`)의 §1–§5 는 상태 문서
`OCRT_MIGRATION_STATUS_2026-08-23.md` 에 누적돼 있으므로 여기서는 반복하지 않는다.

---

## 1. 새 세션 부트스트랩

1. 샌드박스에서 저장소를 받는다(코드·문서·검증 기록 약 130 MB, Mie 자료 없음).

       git clone https://github.com/brtnt/OCRT.git ocrt
       cd ocrt/code/OCRT_C && make        # 샌드박스는 AVX-512 지원 → Makefile 기본값 사용 가능
       bash scripts/run_ipss_gates.sh     # Mie 자료 없이 실행됨. ALL PASS 확인

2. 로컬 폴더 `<LOCAL_ROOT>` 를 연결한다(`<LOCAL_ROOT>` = 이 저장소를 체크아웃한 로컬 작업 폴더, 이하 동일).
   이 폴더가 **정본**이며 저장소 루트와 같다
   (`.git` 이 여기 있다). 로컬 VM 에는 git 2.34 가 있고 github.com 에 닿는다. 삭제 권한은 세션마다
   다시 받아야 한다(git 의 lock 파일 처리에 필요).

3. Mie 표는 GitHub Release `data-v1`(zip 4개, 1.3 GB → 181 표 4.3 GB)에서 받는다:
   `bash scripts/fetch_data.sh` (Windows: `scripts\fetch_data.ps1`). 샌드박스에서 릴리스 다운로드가 막혀 있으면
   로컬 폴더에서 필요한 표만 스테이징한다 — IPSS 벤치·SZA 스윕은 `code/OCRT_C/inputs/M80C.mie` 하나(24 MB),
   OSOAA 정합 하니스는 `C50.mie`, `tsm_ahn/Red_clay_AHN.mie`, `tsm_ahn/Brown_earth_AHN.mie`,
   `water_iop/Detritus_Stramski2001.mie`. 로컬에는 zip 원본이 `package/release_data_v1/` 에 있다.
   설치 후 검증: `cd code/OCRT_C/inputs && sha256sum -c ../../../scripts/MIE_SHA256SUMS_data-v1.txt`.

4. RTSOS 교차검증이 필요하면 `<LOCAL_ROOT>\..\RTSOS` 를 따로 연결한다(이 저장소에 없음, 저자 코드).
   빌드 주의는 상태 문서 §14 와 IPSS 기록 §11.2(링크 순서, `fSnowBRDF_INPUT` 초기화).

5. Windows 실행파일은 사용자가 로컬에서 빌드한다: `code\OCRT_C\build_win.bat`
   (`REBUILD_WIN_v1.11.1_KO.txt`). 이 PC(Core Ultra 9 275HX)는 AVX-512 가 없으므로
   `-march=x86-64-v3` 필수 — `cascadelake` 산출물은 Illegal instruction 으로 죽는다.

## 2. 현재 상태 (v1.11.1, 2026-09-05)

- **구면 보정 = IPSS 단일화.** Zhai & Hu (2022) Eq.(7) κ = I₁,ss/I₁,pp 를 최종 (I,Q,U) 에 곱한다.
  구 평균할선 Chapman PSSA 는 소스 트리에서 삭제됐고 `--pssa` 는 on/off 뿐이다(`--pssa-mode` 는 exit 2).
- **VZA 앵커 = 지상 화소(L1B senz).** TOA 앵커는 태양 천정에서도 −2.9 %(VZA 55)~−11 %(VZA 70) 의
  가짜 편차를 남긴다. 게이트 G-I3 가 지표 붕괴 극한 κ→1 을 상시 검사한다.
- **κ 는 위상함수 가중(v1.11.1).** 혼합 대기(에어로졸 2 km / Rayleigh ~8 km)에서 위상함수는
  상쇄되지 않는다. 층 원천 = β_j[w_R,j P_R(Θ) + w_A,j P_A(Θ)], P_A 는 FR631 값표 우선.
  Rayleigh 단독 산출은 v1.11 과 비트 동일. 게이트 G-I4a/b(독립 scipy 적분과 |Δκ| ≤ 7.5e-6).
- **비-PSSA 경로 불변.** 스모크 SHA `4ced4f18d50e4ac9`, 10케이스 비트 배터리.
- **RTSOS 교차검증.** 같은 Rayleigh 프로파일에서 κ 가 0.1 % 이내 일치, RTSOS 층수 외삽값 −0.239 % vs
  OCRT −0.240 %(VZA 55). RTSOS 결함 2건(미초기화 플래그, 링크 순서)과 이산화 인공항 exp(−Δτ)−1 기록.
- **저 SZA 잔차의 정체.** VZA 55 에서 남는 −0.2 % 는 시선경로 곡률 −(h̄/R_e)tan²θ_v 이며 legacy 가
  빠뜨렸던 항이다(SZA 무관, 천저에서 0). 구면 보정이 실제로 필요한 구간은 SZA 75° 이상.
- 상세: 상태 문서 §14·§14.1, IPSS 기록 문서(§1–§11.6), `docs/reports/` 의 그림 figI1–figI5.
- 참고: `src/main.c` 의 `V2_VERSION` 문자열은 v1.2 계보 표기이고 문서 버전(v1.11.1)과 체계가 다르다.

## 3. 인벤토리 — git 에 있는 것 / 로컬에만 있는 것

git 에 있음(공개): `code/OCRT_C`(src·tests(픽스처 Mie 1개 포함)·scripts·tools·docs·validation·patches·매니페스트),
`code/OCRT_Python`(코드·문서·테스트·소형 data: AFGL·xsec·IOP 표), `scripts/`(자료 설치·매니페스트), `code/campaign_runner`,
`code/run_tools`, `code/validation_05`(harness·bit_baselines·tools), `docs/technical`, `docs/validation`,
`docs/reports`, 상태 문서, 이 README, `LICENSE`, `.gitignore`, `.gitattributes`(`* -text`, 바이트 보존).

로컬에만 있음(`.gitignore`):

| 항목 | 위치 | 크기 | 비고 |
|---|---|---|---|
| Mie 표 `*.mie` (181, 두 트리에 같은 파일) | `code/OCRT_C/inputs/`, `code/OCRT_Python/data/` | 4.3 + 4.3 GB | **GitHub Release `data-v1`** zip 4개(1.3 GB), 로컬 원본 `package/release_data_v1/`; 수령 원본 `package/PART01–PART10*.zip` |
| OSOAA 참조 런 원문 | `code/validation_05/runs/` | 158 MB | 재분석용 소비 파일 |
| 배포 묶음 | `package/` | 2.3 GB | 세션별 tar.gz·zip |
| 캠페인 산출물 | `task/` | 1.9 GB | sens/analysis 등 |
| 투고 전 원고 | `docs/paper/` | 69 MB | RSE 리뷰 논문·OCRT 논문 LaTeX/PDF/그림 |
| 기술 문서 대형 자료 | `docs/technical/**/*.gz` | 64 MB | 편광 legacy 참조 CSV |
| 빌드 산출물 | `build/`, `*.exe` | — | 재빌드 |
| RTSOS | `<LOCAL_ROOT>\..\RTSOS` | — | 저자 코드(CC BY-NC 4.0), 별도 폴더. 교차검증용 수정 드라이버 `main_*.f90`·`ray.pmtx` 도 저장소에서 제외(변경 내용은 `xcheck_rtsos/NOTICE_RTSOS.md`) |
| OSOAA 작업 사본 | `<LOCAL_ROOT>\code\OSOAA` | 4 MB | CNES 코드 — 제3자 자산이라 저장소에서 제외(이력 포함 삭제, 2026-09-06) |

## 4. 재빌드·게이트 (verbatim)

    # Linux / 샌드박스
    cd code/OCRT_C && make
    # AVX-512 없는 호스트(로컬 VM 포함)
    gcc -std=c11 -O3 -march=x86-64-v3 -ffp-contract=fast -fassociative-math -fno-signed-zeros \
        -fno-trapping-math -DOCRT_FAST_KERNELS -fopenmp -Isrc $(find src -name '*.c') -o build/ocrt -lm
    # 게이트 (G-I0~G-I4 + 회귀 §2–§5)
    bash scripts/run_ipss_gates.sh
    # 매니페스트 재생성 (validation/ipss_2026-09-05 와 build/_RELEASE 는 제외 — 기존 규약)
    find . -type f ! -path "./SHA256SUMS.txt" ! -path "./build/*" ! -path "./_RELEASE/*" \
         ! -path "./validation/ipss_2026-09-05/*" -print0 | sort -z | xargs -0 sha256sum \
         | sed 's#  \./#  #' > SHA256SUMS.txt

환경변수: `OCRT_ADVANCED=1`(고급 옵션), `OCRT_DEBUG=1`, `OCRT_IPSS_DUMP=1`(κ 진단),
`OCRT_IPSS_PROFILE_DUMP=<file>`(κ 입력 프로파일·P11 표 덤프), `OCRT_IPSS_NQUAD`(구적 수렴 게이트),
속도 측정은 `OMP_NUM_THREADS=1`.

## 5. 세션 종료 시 GitHub 반영 절차

1. 게이트 ALL PASS 확인, `SHA256SUMS.txt` 재생성, 상태 문서에 절 추가.
2. 로컬 VM 에서 커밋(저장소 루트 = 연결 폴더):

       cd ~/mnt/OCRT && git add -A && git status --short | head
       git commit -m "<요약>"

3. 푸시. 토큰은 저장소 밖에 둔다 — 폴더 루트의 `.github_token`(gitignore 됨) 에 한 줄로 두면
   다음 세션이 대화에 붙여넣지 않고 읽어 쓴다. 토큰은 해당 repo 한정 fine-grained(Contents: Read/Write),
   만료 90일 권장, 유출 의심 시 GitHub 에서 즉시 폐기.

       git -c credential.helper='!f() { echo "username=brtnt"; echo "password=$(cat .github_token)"; }; f' push origin main

   원격 URL 에 토큰을 저장하지 않는다(`origin` 은 `https://github.com/brtnt/OCRT.git`, 토큰 없음).
   저장소: https://github.com/brtnt/OCRT (공개). 첫 푸시 2026-09-05, main = a5655e2 → 9d7b63e.
4. 푸시 뒤 GitHub API 로 커밋 SHA 와 파일 수를 확인하고 이 README §2 의 상태를 갱신한다.

## 6. 규칙 (불변)

- 코드 수정은 사용자 명시 승인 후에만. 게이트 미완이면 원복.
- 위상 표현 정책: 고차 모먼트·L 상향 금지, 표준 = 선형 보간 + 전방 조밀(FR631, 절단 OFF 가 과학 기준).
- 비-PSSA 경로는 비트 불변이 원칙(스모크 SHA·비트 배터리로 증명).
- 정합 비교는 1:1 산포도(오차 밴드) 표준, 속도 비교는 동일 호스트 재실측 "OSOAA 대비 %" 만 유효.
- OSOAA 수중 비교에서 풍속 3 m/s 미만 제외(3 포함).
- RAA 규약: 180° = 경면(글린트), 0° = 후방산란. RT 파장 ≤ 2000 nm.
- 공개 저장소: 원고(`docs/paper`)·자료·토큰은 올리지 않는다. 라이선스는 학술·비상업(`LICENSE`).

## 7. 남은 과제

1. Windows `.exe` 재빌드(v1.11.1 소스 반영) — 사용자 실행.
2. 에어로졸 대기에서 Eq.(7) 근사 자체(단일산란 비의 다중산란 전이)의 전체 구면 RT 대비 정확도 — 논문·RTSOS
   검증은 Rayleigh 뿐. 선택지: RTSOS 에 M80C 위상행렬(.pmtx)을 넣어 2성분 교차검증.
3. Phase II — 결합 대기–해양(`--water-model`)에 Eq.(7) 적용 정책(현재 fail-loud).
4. 과거 `--pssa` 캠페인 산출물(legacy + TOA 앵커) 재실행 판정 — SZA 75° 이상 슬롯이 실질 영향권.
5. 옵션 재분류(production/advanced/debug) P0, 상태 문서 §1–§13 의 미결 과학 게이트(raw FR631 convergence
   freeze, coupled-LUT state 복구, Rrs-Q 참조물리 정합, Chl/detritus OSOAA end-to-end).
6. `src/main.c` `V2_VERSION` 문자열 정리(문서 버전 체계와 일치).
