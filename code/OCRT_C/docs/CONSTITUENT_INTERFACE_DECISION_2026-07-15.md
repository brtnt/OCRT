# OCRT 구성성분 입력 인터페이스 결정 자료

- 작성일: 2026-07-15 KST
- 상태: **설계 승인 대기**
- 적용 기준선: OCRT v1.14 migration G5-2

## 1. 확인된 요구와 현재 구현의 차이

첨부된 `OCRT_methods_algorithm_v02`는 해양 입력 파라미터로 Chl, TSM,
`a_cdom`, `S_cdom`, PSD slope와 complex refractive index를 제시하고,
대표 명령 예시에 `--chl`, `--tsm`, `--acdom`, `--sdom`, `--psd`, `--m`을
사용한다. 반면 v1.13 이관 체크포인트는 통합 `--chl/--tsm` 경로의 구현
방식을 별도 의사결정 항목으로 남겼다.

현재 검증된 실행 경로는 다음 두 가지다.

1. SIMPLE constituent 경로: `--simple-chl`, `--simple-min`,
   `--simple-adom440`, `--simple-adom-s`가 수중 IOP/위상으로 변환된다.
   v1.14에서는 명시적 `--surface ocean`을 요구하고, 대기·풍속·선글린트
   설정을 덮어쓰지 않는다.
2. 전처리 + fixed-bulk 주입 경로: `scripts/osoaa_chl_compare.py`와
   `scripts/osoaa_tsm_compare.py`가 IOP/위상 LUT를 생성하고,
   `--fixed-bulk-iop` 및 `--fixed-bulk-phase-lut`로 주입한다. 이 경로는
   골든 항목3·5에서 각각 5/5 셀을 통과했다.

또한 현재 패키지에서 `--simple-chl > 0` 실행은
`aph_bricaud_1998.txt` 부재로 종료된다. 따라서 CLI 이름만 `--chl`로
노출하는 것은 방법론 문서의 통합 경로가 동작한다고 오해하게 만들 수 있다.

## 2. 선택지

### A. 전처리 + fixed-bulk를 운영 경로로 유지

- 장점: 현재 골든과 OSOAA 비교가 재현되며, IOP와 위상함수 산출물을
  케이스별로 감사할 수 있다.
- 단점: 사용자 관점에서 `--chl/--tsm` 단일 명령이 아니며, 중간 산출물
  관리가 필요하다.

### B. OCRT 내부에 `--chl/--tsm` 통합

- 장점: 방법론 문서와 CLI가 직접 일치하고 사용성이 단순해진다.
- 필수 선행조건:
  - `aph_bricaud_1998.txt`의 출처·라이선스·버전과 패키징 위치 확정
  - Chl/TSM별 IOP 및 입자위상 구성식 단일화
  - `--simple-*`와 새 canonical 옵션의 관계(별칭/폐기/독립 모델) 확정
  - 412/490/555/660/865 nm, SZA 0/40/80° 골든 추가
  - OSOAA 전처리 경로와의 수치 동등성 허용오차 확정

## 3. v1.14 권고

**v1.14에서는 A를 운영 경로로 유지하고, SIMPLE은 검증·개발 경로로
명시한다.** `--chl/--tsm` public 별칭은 필요한 광학자료와 골든 게이트가
준비되기 전에는 추가하지 않는다. 이 권고는 이관 완료를 위한 보수적
기준선이며, 통합 인터페이스의 최종 선택은 별도 승인 사항이다.

## 4. 통합 경로 승인 시 완료 기준

- 패키지 단독 실행에서 외부 미동봉 파일 오류가 없어야 한다.
- `--chl/--tsm` 도움말, 사용 가이드, 방법론 문서의 단위와 기본값이
  일치해야 한다.
- 기존 전처리 fixed-bulk 결과와 새 내장 결과의 IOP, phase moments,
  Rrs, TOA Stokes 비교표를 생성해야 한다.
- Tier-0, 골든 항목2·3·5·6, 해양 격자 13행, OSOAA 11/11을 모두
  재통과해야 한다.
