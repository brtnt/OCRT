# OCRT 안전 커널·LUT 불변량 업데이트 — 2026-07-29

버전: `OCRT-v1.2-2026-07-29-KST-safe-kernel-invariant-update`

## 목적

`MIGRATION_WORKORDER_2026-07-29_KO.md`에서 제안된 변경 중, 현재 2026-07-28 canonical 코드의 물리 출력과 부동소수점 결과를 바꾸지 않는 항목만 반영한다.

## 적용 항목

### 1. SOS 영가중 입력열 건너뛰기

관측방향으로 삽입된 zero-weight ordinate는 source contraction에 정확히 0을 기여한다. 스칼라·벡터 SOS source operator가 해당 입력열을 즉시 건너뛰도록 했다.

진단용 비활성화:

```bash
OCRT_SOS_ZERO_COL_SKIP_OFF=1
```

새 회귀시험 `scripts/test_sos_zero_col_skip.sh`는 skip ON/OFF의 스칼라 및 I/Q/U source 배열을 double bit pattern으로 비교한다.

### 2. Native coupled angular LUT의 Beer 항 hoist

`rt_solve_case_ocean_lut()`에서 다음 항을 실제 의존 차원으로 이동했다.

- `exp(-tau/mu_air[VZA])`: RAA loop 밖, VZA loop 수준
- `exp(-tau/mu_sun)` 및 PSSA direct-glint scale: 전체 VZA×RAA loop 밖

연산식과 각 셀에 저장되는 값은 변경하지 않았다.

## 적용하지 않은 항목

### 수중 SOS geometric-tail acceleration

근사이며, 첨부 문서에는 이식에 필요한 전체 guard·state 정의와 원본 diff가 없다. 기본 결과 보존 요구 때문에 반영하지 않았다.

### `--precision standard` 기본값

첨부 문서가 명시하듯 옵션 없는 실행 결과를 변경한다. 현재 기본 동작은 기존 reference 설정으로 유지했다. 또한 각 mode의 정확한 노브 매핑이 첨부 문서에 완전하게 제공되지 않았다.

### Python zero-column skip

Python/NumPy prototype은 물리적으로 roundoff 수준이지만 4개 행 15개 필드에서 마지막 비트가 바뀌었다. 최대 절대차는 `2.168404344971009e-19`였다. byte-identical 조건을 충족하지 않아 원복했다. Python coupled LUT의 Beer 항은 이미 VZA 배열로 RAA loop 밖에서 계산된다.

## 검증 요약

- 최고 탁도 12-cell coupled LUT: 기준본·수정본·skip OFF CSV SHA-256 동일
- PSSA coupled LUT: 기준본·수정본 CSV SHA-256 동일
- native LUT 진단: `water_cold=1`, `water_views=3`, `cells=12`
- S17A cache ON/OFF: 물리 출력 byte-identical; `T_diff_dn_dir`는 cache path 전용 diagnostic이므로 비교에서 제외
- 3회 교차 성능 측정 중앙값: 3.298239 s → 3.258898 s, 약 1.19% 단축
- Release 및 non-fast build: PASS
- 주요 수중·계면·RAA·PCHIP·IOP/Kd 회귀시험: PASS

