# OCRT v1.19 CCRR Chl 자료 통합 완료 보고서

## 릴리스

```text
OCRT-v1.19-2026-07-18-KST-ccrr-chl-morel-mm01
```

기준선은 OCRT v1.18 해수 입력 분기 계약이다. RT solver, TSM, OCRT organic Chl,
aDOM 및 직접 IOP 수치 경로는 변경하지 않았다.

## 입력자료 판정

첨부 파일은 440 nm에서 1로 정규화된 `A_chl(lambda)` 스펙트럼이다. 파일 내부
설명은 Morel (1988, Fig. 10c; Prieur and Sathyendranath 1981)을 출처로 명시하며,
300–350 nm와 700–1000 nm를 외삽 구간으로 경고한다.

파일명에 `bricaud_2011`이 들어가지만, 이 자료를 Bricaud et al. (1998)의
파장별 `A(lambda), E(lambda)` 표로 취급하지 않았다.

## OCRT 인터페이스 변환

현재 CCRR Chl 로더가 소비하는 식은 다음이다.

```text
a_p(lambda) = Aphi(lambda) * Chl^Ephi(lambda)
```

classic Case-1 폐쇄를 보존하도록 다음으로 변환했다.

```text
Aphi(lambda) = 0.06 * A_chl_normalized(lambda)
Ephi(lambda) = 0.65
```

결과적으로:

```text
a_p(lambda) = 0.06 * A_chl_normalized(lambda) * Chl^0.65
```

5열 파일의 `Ap`, `Ep`는 현재 코드가 읽지 않으므로 명시적 0 placeholder다.

## 파일

```text
canonical : ocrt/inputs/water_iop/aph_ccrr_morel1988_mm01.txt
legacy    : ocrt/inputs/water_iop/aph_bricaud_1998.txt
source    : ocrt/inputs/water_iop/source/apstarchl_morel1988_normalized_user_supplied.txt (LF-normalized copy)
generator : ocrt/scripts/build_ccrr_chl_table.py
```

canonical과 legacy 사본은 byte-identical이다.

## 대표 수치

| 조건 | 기대 `a_pig` | OCRT 출력 |
|---|---:|---:|
| 440 nm, Chl=1.0 | 0.060000000 | 0.060000000 |
| 443 nm, Chl=0.3 | 0.0257791915 | 0.0257791915 |
| 555 nm, Chl=0.3 | 0.00568146277 | 0.00568146277 |
| 670 nm, Chl=0.3 | 0.0162680223 | 0.0162680223 |

443 nm 값은 440–445 nm `A_chl` 선형 보간을 사용한다.

## 검증

| 검사 | 결과 |
|---|---:|
| 141개 원자료 행 변환 | PASS |
| 440 nm 정규화와 계수식 | PASS |
| generator 재실행 byte identity | PASS |
| canonical/legacy 파일 byte identity | PASS |
| canonical-only 로드 | PASS |
| legacy fallback 로드 | PASS |
| CCRR Chl runtime anchor | 4/4 PASS |
| CCRR Chl 자동시험 | 11/11 PASS |
| production CCRR Chl 실행 | PASS, 26 orders, converged |
| batch-full-grid CCRR Chl | PASS |
| v1.18 대 v1.19 수치 비교 | 10/10 stdout·stderr byte-exact |
| water branch 계약 | 14/14 PASS |
| Ahn TSM | 6/6 PASS |
| OCRT organic Chl | 10/10 PASS |
| TSM phase-cache | PASS, 기존 trace SHA 유지 |
| 수정 `rt_water_iop.c` strict compile | PASS |
| ASan/UBSan 대표 실행 | finding 0 |
| production build | 오류·경고 출력 0줄 |

TSM phase-cache trace SHA-256:

```text
89146542445cd110649942f6b09ce559f1a95225ed0fbc9ca9b5e13d1994b833
```

## 결과 영향

v1.18 바이너리와 v1.19 바이너리를 같은 새 자료파일로 실행하여 다음 경로를
비교했다.

```text
OCRT pure water
CCRR pure water
OCRT Chl
OCRT TSM
OCRT Chl+TSM+aDOM
CCRR TSM+aDOM
CCRR Chl
CCRR Chl+TSM+aDOM
direct IOP
CCRR Chl at 670 nm
```

모든 조건에서 종료코드, stdout, stderr가 byte-identical이었다. 즉 소스 로더의
canonical-name 정리와 유효성 검사는 수치 결과를 바꾸지 않았다. 새로 가능해진
물리 경로는 종전에 자료 누락으로 실패하던 양수 CCRR Chl뿐이다.

## 남은 과학적 경계

1. 이 통합은 CCRR Chl **흡광자료 누락**을 해결한다. CCRR 입자 산란은 기존 MM01
   `b`, `bb` 식을 사용한다.
2. CCRR 입자 위상은 기존 scalar positive Fournier–Forand P11 폐쇄를 유지한다.
   OCRT organic 분기의 EAP/유기 detritus P11/P12/P33 vector Mie와는 별도다.
3. 300–350 nm와 700–1000 nm 값은 첨부 원자료의 외삽값이다. 해당 구간을 독립
   검증된 흡광자료로 간주하면 안 된다.
4. 향후 실제 Bricaud 1998 `A(lambda), E(lambda)` 표를 도입한다면 별도 모델명과
   별도 자료파일로 추가해야 하며, 현재 Morel/MM01 표를 덮어쓰면 안 된다.
