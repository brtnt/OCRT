# OCRT v1.18 해수 입력 분기 구현·검증 기록

## 기준선과 범위

- 기준선: OCRT v1.17 organic-Chl candidate (`78f02fc`)
- 구현 버전: `OCRT-v1.18-2026-07-18-KST-water-branch-contract`
- 변경 범위: CLI/CSV 입력 계약, 입력 검증, 기존 실행 필드로의 lowering
- 변경하지 않은 범위: 순수해수 IOP, OCRT organic/TSM 광학모델, CCRR 변환식,
  직접 IOP closure, vector SOS 수치 kernel

## 구현 계약

`--surface ocean`에서는 다음 중 하나를 정확히 한 번 지정한다.

```text
--water-model ocrt
--water-model ccrr
--water-model iop
```

OCRT·CCRR은 각각 `--ocrt-*`, `--ccrr-*` 입력만 허용하고 직접 IOP는
`--iop-*`만 허용한다. 모델 미지정, 중복 지정, 다른 분기 접두사 혼용은
종료코드 2다.

OCRT 또는 CCRR에서 Chl, TSM, aDOM440을 모두 명시적으로 0으로 주면
구성성분 변환기를 실행하지 않고 native pure-water 경로로 하강한다.

## 빌드

생산용 빌드 명령:

```bash
gcc -std=c11 -O3 -march=native -ffp-contract=fast -fassociative-math \
  -fno-signed-zeros -fno-trapping-math -DOCRT_FAST_KERNELS -fopenmp \
  -Isrc $(find src -name '*.c' -print | LC_ALL=C sort) \
  -o build/v2_solver_vk_v1.18 -lm
```

결과:

```text
build log lines = 0
elapsed_seconds = 15.77
binary SHA-256  = 487ceb9e6fade32731b04d39f064642e639ea5d49101c3ddcfd098169ec1ee4d
```

전체 트리를 `-Wall -Wextra -Wpedantic -Werror`로 빌드하는 것은 기준선부터
존재하는 GNU statement-expression, misleading-indentation, unused-symbol 경고
때문에 통과하지 않는다. 이번에 추가된 organic/TSM 독립 모듈 strict 검사와
실제 production build는 통과했다.

## 분기 계약 회귀시험

`scripts/smoke_cli_water_branch_v118.sh`:

```text
no model                         PASS
model duplicated                 PASS
cross-prefix options             PASS
missing required constituent     PASS
legacy unprefixed option         PASS
water options on non-ocean       PASS
OCRT species without TSM         PASS
incomplete IOP                   PASS
multiple IOP phase sources       PASS
orphan IOP Mie control           PASS
OCRT 0/0/0 pure-water lowering   PASS
CCRR 0/0/0 pure-water lowering   PASS
OCRT/CCRR zero stdout identity   PASS
IOP branch                       PASS
TOTAL                            14/14 PASS
```

OCRT 0/0/0과 CCRR 0/0/0의 stdout SHA-256은 동일하다.

```text
59e8956d5e5532030e66606ebe2646cea054e94379d7c48ff31f90c51e8524bc
```

## v1.17 수치 결과 보존

기존 v1.17 옵션과 대응하는 v1.18 접두사 옵션을 동일한 물리 입력으로 실행하고
stdout을 byte 단위로 비교했다. 모두 동일했다.

| 물리 경로 | stdout SHA-256 | 결과 |
|---|---|---:|
| 기존 native pure water → OCRT 0/0/0 | `59e8956d5e5532030e66606ebe2646cea054e94379d7c48ff31f90c51e8524bc` | exact |
| OCRT Chl | `3a4e5b97d7742e0e9495db7adc1e27d2200a7734eb7e094ba5347e3e1b115920` | exact |
| OCRT Ahn TSM | `dc7e64a8991c63212ab0758d9c7d6d58fc5f881d88b6a47a38f63bcb951283ee` | exact |
| CCRR TSM+aDOM (`Chl=0`) | `ab4c7b5ae895b746b6a95eb91bb5a8f88b0b89d16e1d6029f133c66acedaed9d` | exact |
| 직접 IOP | `da4f58307a7a084db3d8381533b3cd17570ff8ae9f8a96f097bcd08e69307c3b` | exact |

주의: v1.17의 “구성성분 모드 0/0/0”을 보존한 것이 아니라, 사용자 요구대로
v1.17의 **native pure-water** 결과에 하강하도록 바꾼 것이다. 이것은 의도된
분기 의미 변경이며 위 표의 첫 행으로 검증했다.

## 기존 광학모델 회귀시험

| 시험 | 결과 |
|---|---:|
| OCRT organic Chl smoke | 10/10 PASS |
| Ahn TSM smoke | 6/6 PASS |
| TSM phase-cache direct/fallback equivalence | PASS |
| phase-cache build count의 SOS 차수 독립성 | PASS |
| phase-cache trace SHA-256 | `89146542445cd110649942f6b09ce559f1a95225ed0fbc9ca9b5e13d1994b833` |
| shell syntax | PASS |
| modified Python scripts `py_compile` | PASS |
| PowerShell syntax execution | 미실행 (`pwsh` 미설치) |

## batch full-grid

- OCRT base row, CCRR 0/0/0 row, IOP row: 3/3 성공
- OCRT 0/0/0 CSV와 CCRR 0/0/0 CSV: byte-identical
- 서로 다른 접두사를 섞은 행: 종료코드 2로 거부

## 확인된 잔여 제약

CCRR의 양수 Chl 변환 함수는 구현되어 있으나 다음 자료가 현재 패키지에 없다.

```text
inputs/water_iop/aph_bricaud_1998.txt
```

따라서 `--ccrr-chl > 0`은 현재 다음 오류로 종료한다.

```text
rt_iop_ccrr: cannot open aph_bricaud_1998.txt
```

이번 변경에서는 출처·버전·라이선스가 확정되지 않은 계수표를 임의 생성하거나
OCRT organic 자료로 대체하지 않았다. CCRR `Chl=0`의 TSM/aDOM 경로와 CCRR
0/0/0 순수해수 하강 경로는 정상 동작한다.
