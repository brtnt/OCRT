# S_CDOM / S_detritus 제어 인터페이스 — 현황 정리 및 Python 페어 옵션 추가 (2026-08-29)

질의: "해수 광특성 입력 시 S_detritus·S_CDOM을 제어하는 옵션이 있나? 없으면 추가."

결론: **C에는 이미 있었고(추가 불필요), Python 페어는 라이브러리에만 있고 CLI·배치 경로에 미노출이어서 이번에 추가했다.**

공통 물리 규약(양 코드 동일): `a(λ) = a(440)·exp(−S·(λ−440))`, 기준파장 440 nm 고정.

## 1. C — 기존 옵션 (변경 없음)

| 양 | 옵션 | 기본값 |
|---|---|---|
| S_CDOM | `--ocrt-adom-slope S` / `--ccrr-adom-slope S` | 0.014 nm⁻¹ (Bricaud 1981) |
| S_detritus | `--ocrt-detritus-slope S` | 0.0109 nm⁻¹ (Bricaud–Stramski) |
| 크기 | `--ocrt-adom440`, `--ocrt-detritus-a440` | 0 |

`--help` 문서화 완료, `OCRT_ADVANCED` 게이트 없음(일반 옵션). 배치 CSV 열 `ocrt_adom_slope`·`ocrt_detritus_slope` 지원. 제약: detritus 계열은 `--ocrt-chl > 0` 필요(Chl 연동 설계, 아니면 fail-loud), CCRR 패밀리에는 detritus 없음, `--water-model iop`에서는 무관. 구 옵션명(`--cdom-slope`, `--adom-slope`, `--detritus-slope`)은 제거되어 사용 시 오류.

실측(412 nm): CDOM a440=0.1 → S=0.014에서 a_dom 0.147994, S=0.020에서 0.175067 (이론값과 6자리 일치). detritus a440=0.05 → S=0.0109에서 기여 0.067844, S=0.015에서 0.076098 (이론 0.067845 / 0.076098).

## 2. Python — 진단 및 추가한 것

진단: `constituent.py`는 C와 동일 기본값·동일 법칙으로 `evaluate(..., cdom_slope=, det_slope=)`를 이미 지원했고 `ocrt_solve.py`의 solve 함수도 인자를 받았으나, **argparse에 플래그가 없어 CLI로는 조절 불가**였고 **`batch_driver.py`가 슬로프를 전달하지 않아 배치 경로에서는 지정해도 무시**됐다.

추가 (3파일, diff 211줄 — `docs\technical\py_constituent_slope_options_2026-08-29.diff`):

- `ocrt_solve.py`
  - CLI 플래그 `--adom-slope`(별칭 `--cdom-slope`), `--detritus-a440`, `--detritus-slope` — 기본값은 C와 동일(0.014 / 0 / 0.0109).
  - `cmd_single`·`cmd_grid`에서 solve 경로로 전달.
  - grid CSV **선택 열** 인식(행별 우선, 없으면 CLI값): `ocrt_adom_slope|adom_slope|cdom_slope|s_cdom`, `ocrt_detritus_a440|detritus_a440|det_a440`, `ocrt_detritus_slope|detritus_slope|det_slope|s_detritus` — C 배치 CSV 철자를 그대로 수용.
  - 출력 CSV 말미에 `adom_slope, detritus_a440, detritus_slope` 열 추가(사용값 기록, 기존 열 순서 불변).
- `ocrt_py/batch_driver.py`
  - `_case_slopes()` 헬퍼 신설: 케이스 dict의 `cdom_slope`/`det_a440`/`det_slope`를 읽고 없으면 C 기본값.
  - 배치 IOP 고속경로는 **파장 + 슬로프가 모두 동일할 때만** 사용하고, 혼합 슬로프는 케이스별 `evaluate()`로 폴백.
  - `evaluate_batch()`·`evaluate()` 호출 양쪽에 슬로프 전달.
- `produce_grid.py`
  - grid CSV의 선택 열을 읽어 케이스 dict에 실어 보냄.

## 3. 검증

| 항목 | 결과 |
|---|---|
| C↔Py IOP 정합 (412 nm) | CDOM S=0.014/0.020, detritus S=0.0109/0.015 전부 **차 ≤2.5e-7** (C 출력 6자리 반올림 한계) |
| `evaluate_batch` vs `evaluate` | a_det 0.0760980778 **완전 일치** |
| CLI 플러밍 (R1 결합해 실행, 412 nm, aDOM440=0.1) | S=0.014 → Rrs 1.1043e-03, S=0.020 → 9.3946e-04 (흡수 증가 → Rrs 감소, 물리 방향 정상) |
| **C↔Py 슬로프 응답비** (0.014→0.020) | Rrs C 0.850720 / Py 0.850724 (**상대차 0.0004%**), rrs(0−) C 0.851134 / Py 0.851135 (**0.0002%**) |
| **무회귀** | 패치 전 코드와 기본값 실행 결과 **자릿수 동일** (Rrs 1.1043042958e-03, rrs 2.0921902692e-03) |
| 기본값 명시/생략 동등성 | `evaluate` 전 키 완전 일치 |
| 문법 | 3파일 `ast.parse` OK |

응답비 비교는 절대값 오프셋(파이썬 실행이 저해상 설정: n_mu_water 24 / m_max 2 / n_layers 40 / gas off)을 상쇄한 지표로, 슬로프 옵션이 양 코드에서 동일하게 작용함을 보인다.

## 4. 사용 예

```powershell
# C
ocrt.exe --water-model ocrt --ocrt-chl 0.5 --ocrt-adom440 0.1 --ocrt-adom-slope 0.019 `
         --ocrt-detritus-a440 0.05 --ocrt-detritus-slope 0.012 ...

# Python
python ocrt_solve.py single --wl 412 --sza 40 --vza 20 --raa 90 --wind 3 `
       --chl 0.5 --acdom440 0.1 --adom-slope 0.019 `
       --detritus-a440 0.05 --detritus-slope 0.012 --phyto-group micro
```

grid CSV에 `ocrt_adom_slope`, `ocrt_detritus_a440`, `ocrt_detritus_slope` 열을 넣으면 행별로 다른 슬로프를 쓸 수 있다(열이 없으면 기존과 동일).

## 5. 남은 항목

- Python `batch_driver` 혼합-슬로프 폴백은 정확하지만 배치 고속경로를 못 쓴다. 슬로프를 벡터화하려면 `evaluate_batch`의 `cdom_fac`/`det_fac`를 (B,) 배열로 확장해야 한다 — 슬로프 스윕 캠페인을 대량으로 돌릴 때만 필요.
- C 쪽 detritus 제약(`--ocrt-chl > 0`)은 설계 의도이나, CDOM처럼 Chl 독립으로 쓸 수요가 있으면 별도 논의 필요.
