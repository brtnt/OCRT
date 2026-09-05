# OCRT v1.19 변경 기록 — 2026-07-18

실행 식별자:

```text
OCRT-v1.19-2026-07-18-KST-ccrr-chl-morel-mm01
```

## 변경

- 사용자 제공 Morel (1988) 440-nm 정규화 Chl 흡광 형상을 OCRT 5열 CCRR
  인터페이스로 변환했다.
- 변환식은 `Aphi=0.06*A_chl`, `Ephi=0.65`이며, 따라서
  `a_p=0.06*A_chl*Chl^0.65`이다.
- canonical 파일 `aph_ccrr_morel1988_mm01.txt`를 추가했다.
- v1.18 이하 호환을 위해 byte-identical `aph_bricaud_1998.txt` 사본을 유지했다.
- 로더는 canonical 파일을 먼저 읽고 역사적 이름을 fallback으로 사용한다.
- 파장 단조성 및 계수 유효성 검사를 추가했다.
- CCRR 양수 Chl 경로의 “자료 미패키징” 제약을 제거했다.

## 비변경

- OCRT organic Chl 모델은 변경하지 않았다.
- Ahn TSM 4종과 vector P11/P12/P33 경로는 변경하지 않았다.
- aDOM 모델과 OCRT vector SOS solver는 변경하지 않았다.
- 기존 OCRT/TSM/IOP 및 순수해수 수치 경로는 변경하지 않았다.

## 명칭

원자료 파일명에 `bricaud_2011` 문자열이 있으나, 파일 내부 설명은 Morel (1988)
정규화 스펙트럼이다. 따라서 새 자료를 Bricaud 1998의 파장별 A/E 데이터로
표기하지 않는다.
