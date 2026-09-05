# OCRT Mie Stage 3C 릴리스 노트

버전: `OCRT-v1.2-2026-08-16-KST-mie-direct-cache-allview`

## 변경 사항

- 현재 worker에서 사용하는 완전한 `.mie` 모델을 모든 bulk/phase 파장과 함께 RAM에 보존한다.
- 수중 입자의 direct θ-linear vector Fourier kernel을 다중-entry worker-private cache에 보존한다.
- 고정 방향격자의 Fourier 삼각함수와 FR631 geometry mapping을 worker-local cache로 재사용한다.
- 한 case는 native coupled all-view 경로에서 모든 VZA×RAA를 산출한다.
- production에서 VZA×RAA cell별 solver replay를 금지한다. 진단용 override만 `OCRT_ALLOW_LEGACY_CELL_REPLAY=1`로 허용한다.
- 독립 batch row가 이전 row의 대기 S7/S7b all-view cache를 상속하지 않도록 case 시작 시 해당 cache를 초기화한다.
- 병렬 batch의 각 row는 private solver/work/cache 상태를 사용하며 serial/parallel 결과가 byte-identical인지 회귀시험한다.

## 물리 결과

이번 Stage 3C는 cache 및 실행구조 변경이다. 현재 active Mie 자료와 동일 입력에서 Stage 3B 대비 443/555/660 nm full-grid CSV가 byte-identical이다.

## 검증

- C cache ON/OFF 결과 byte-identical
- native all-view `replay=0`
- 독립 4-row batch serial/parallel 결과 byte-identical
- active post-Z09 pure-water reference byte-identical
- 2026-07-26 OSOAA handoff summary 재계산 결과 max absolute difference 0
- Python 전체 시험 50 passed, 알려진 VZA=0 편광 추출 경고 1건

## 범위 제한

- Stage 2에서 재생성한 FR631 aerosol/EAP 자료의 active-tree 설치와 forward truncation 최종정책은 후속 단계이다.
- 2026-07-26 TSM campaign은 현재 AHN TSM scalar IOP 및 phase contract와 matched input이 아니므로 현재 AHN 결과의 절대 기준으로 사용하지 않는다.
