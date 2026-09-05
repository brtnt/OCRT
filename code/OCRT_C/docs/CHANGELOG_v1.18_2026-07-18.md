# OCRT v1.18 변경사항 — 해수 입력 분기 계약

- `--surface ocean`에서 `--water-model ocrt|ccrr|iop`를 정확히 한 번 요구한다.
- OCRT 구성성분 옵션을 `--ocrt-*`로 통일했다.
- CCRR 구성성분 옵션을 `--ccrr-*`로 통일했다.
- 직접 IOP 옵션을 `--iop-*`로 통일했다.
- OCRT/CCRR의 Chl, TSM, aDOM440은 모두 명시해야 하며 0은 유효하다.
- OCRT 또는 CCRR에서 세 값이 모두 0이면 동일한 native pure-water path를 실행한다.
- legacy SIMPLE, unprefixed constituent, native CDOM, fixed-bulk 이름은 오류로 종료한다.
- full-grid CSV에도 동일한 분기·접두사 규칙을 적용했다.
- public mode provenance를 `rt_water_input_mode_t`로 solver 입력까지 전달한다.
- 기존 OCRT/CCRR/IOP 수치 kernel은 변경하지 않고 parser/validation lowering만 정리했다.
- CCRR 양수 Chl 경로의 기존 `aph_bricaud_1998.txt` 미패키징 제약을 문서화했다.
