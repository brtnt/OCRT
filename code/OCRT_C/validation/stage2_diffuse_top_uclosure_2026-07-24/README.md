# Stage-2 diffuse-top incoming-U closure validation

This directory records the validation of
`OCRT-v1.2-2026-07-24-KST-stage2-interface-fix123-dtp-uclosure`.

The production change negates only the molecular incoming-U coefficients `ruI`, `ruQ` and `ruU`
inside `ocrt_add_diffuse_top_primary`. The sign of the existing U-output wrapper is kept.

Acceptance evidence:

- reference-free source/operator closure for m = 0..4, every incident node and the I/Q/U basis;
- pure Rayleigh, maximum absolute closure error: 2.776e-17;
- mixed Rayleigh–aerosol, maximum absolute closure error: 5.551e-17;
- 72/72 physical runs and 600/600 matched cells;
- no runtime increase relative to FIX123;
- six independent pure-water Rrs/rrs I/Q/U scatter plots.

The attached harness from another session remains diagnostic, because it does not pass an
identical physical water depth to both codes and its 412 nm water SOS run can stop at order 20.
The production 72-run gate uses the earlier exact input contract and raises the water order cap
explicitly.

---

# Stage-2 확산 상단 입사-U 폐합 검증

이 디렉터리는 `OCRT-v1.2-2026-07-24-KST-stage2-interface-fix123-dtp-uclosure` 의 검증을 기록한다.

생산 변경은 `ocrt_add_diffuse_top_primary` 안의 분자 입사-U 계수 `ruI`, `ruQ`, `ruU` 의 부호만 뒤집는다.
기존 U 출력 래퍼의 부호는 유지한다.

승인 근거:

- m = 0..4, 모든 입사 절점, I/Q/U 기저에 대한 기준값 없는 원천/연산자 폐합;
- 순수 Rayleigh 최대 절대 폐합 오차: 2.776e-17;
- Rayleigh–에어로졸 혼합 최대 절대 폐합 오차: 5.551e-17;
- 물리 실행 72/72, 대응 셀 600/600;
- FIX123 대비 실행 시간 증가 없음;
- 독립 순수해수 Rrs/rrs I/Q/U 산포도 6종.

첨부된 다른 세션의 하니스는 두 코드에 같은 물리적 수심을 주지 않고 412 nm 수중 SOS 실행이 20차에서 멈출 수
있으므로 진단용으로 남긴다. 생산 72회 게이트는 이전의 정확한 입력 계약을 쓰고 수중 차수 상한을 명시적으로
올린다.
