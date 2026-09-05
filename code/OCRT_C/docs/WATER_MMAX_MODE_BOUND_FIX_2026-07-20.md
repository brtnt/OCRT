# OCRT v1.2 water Fourier mode-bound hotfix

- Issue: `OCRT-RT-C3D-EXTBOTTOM-001`
- Date: 2026-07-20
- Release identifier: `OCRT-v1.2-2026-07-20-KST-water-mmax-mode-bound-fix`
- Base: `OCRT-v1.2-2026-07-19-KST-pssa-numeric-fixes`

## 1. 결론

첨부된 결함 보고는 타당하다. Coupled-ocean C3d pass에서 수중 water-leaving
Fourier source는 `min(atmosphere_m_max, water_m_max)+1`개 mode만 할당했으나,
두 번째 대기 SOS pass는 atmospheric `m_max`까지 그대로 순회했다. 따라서
`water_m_max < atmosphere_m_max`이면 external bottom-source 배열 끝을 넘는
heap read가 발생할 수 있었다.

Aerosol-on 기본 대기 상한은 `m_max=16`이므로 `--water-m-max 2`, `4`, `8`에서
TOA I/Q/U가 폭주하면서도 return code 0과 `conv=1`을 반환할 수 있었다.

## 2. 직접 재현

대표 조건:

```text
wavelength          = 443 nm
SZA/VZA/RAA          = 30/30/90 deg
atmosphere n_mu      = 48
water n_mu           = 48
M80C AOD555          = 0.15
OCRT Chl/TSM/aDOM440 = 0.3 / 0 / 0.02
atmospheric m_max    = 16 (aerosol auto)
```

수정 전 production convergence 결과:

| water_m_max | TOA I | TOA Q | TOA U | orders/conv |
|---:|---:|---:|---:|---:|
| 2 | `-3.173e238` | `-1.868e239` | `8.947e244` | `64/1` |
| 4 | `1.433e40` | `-2.234e41` | `-5.090e40` | `64/1` |
| 8 | `-3.050e19` | `9.040e19` | `1.279e20` | `64/1` |
| 16 | `1.2681773953e-01` | `3.4351814993e-03` | `-2.5667208100e-02` | `64/1` |

수정 전 ASan 실행은 `rt_surface_boundary.c:add_external_bottom_source()`에서
heap-buffer-overflow를 검출했다.

## 3. 수정 내용

### 3.1 Full atmospheric allocation and zero padding

C3d water-leaving source 배열을 항상 다음 크기로 할당한다.

```text
(atmospheric_m_max + 1) * atmospheric_n_mu
```

수중 solver가 실제로 생성한 `0..water_source_m_max`만 채우고 더 높은 대기
mode는 `calloc`의 정확한 0으로 유지한다. Fourier mode는 서로 독립이므로
수중 source가 없는 상위 mode를 0으로 처리하는 것이 물리적으로 맞다.

### 3.2 External source shape contract

`rt_options_t`에 다음 metadata를 추가했다.

```c
int ext_bottom_n_mu;
int ext_bottom_m_max;
```

I/Q/U 포인터는 하나의 shaped object로 취급한다. 대기 solver는 mode slice를
만들기 전에 다음을 검사한다.

- I/Q/U 세 포인터가 모두 존재하거나 모두 없어야 함
- source stride와 실제 atmospheric direction count가 같아야 함
- source mode bound가 유효해야 함

불일치는 `RT_SOLVER_ERR_BOTTOM_SOURCE_SHAPE (-7)`로 종료하며 배열을 읽지 않는다.

### 3.3 Per-mode bounds

각 Fourier mode에서 external source pointer를 먼저 `NULL`로 초기화한다.
`m <= ext_bottom_m_max`인 경우에만 slice를 연결한다. 더 높은 mode는 exact-zero
source로 계산된다.

### 3.4 S7b cache shape

S7b atmosphere LUT solve는 `view_as_node=0`인 순수 Gauss ring을 사용한다. 일반
row path가 zero-weight view node를 포함하는 경우에는 별도의 Gauss-only
water-to-air projection을 한 번 구성한다. 이 source와 shape metadata 전체를
fingerprint에 포함한다. 따라서 현재 VZA가 cache source stride를 바꾸지 않는다.

Cache fill 실패 시 water-leaving source를 누락하지 않고 exact row path로
fallback한다.

## 4. 수정 후 production convergence 결과

| water_m_max | TOA I | TOA Q | TOA U | orders/conv |
|---:|---:|---:|---:|---:|
| 2 | `1.2650810514e-01` | `3.4298673631e-03` | `-2.5676164203e-02` | `64/1` |
| 4 | `1.2681917028e-01` | `3.4353972982e-03` | `-2.5667430761e-02` | `64/1` |
| 8 | `1.2671644563e-01` | `3.4332396094e-03` | `-2.5667450738e-02` | `64/1` |
| 16 | `1.2681773953e-01` | `3.4351814993e-03` | `-2.5667208100e-02` | `64/1` |
| 30 | `1.2681773953e-01` | `3.4351814993e-03` | `-2.5667208100e-02` | `64/1` |

모든 값은 finite이고 broad corruption guard `abs(TOA I/Q/U)<100`을 통과한다.
`m=16,30,48,64`는 수정 전 정상 baseline과 stdout/stderr byte-identical이다.

## 5. 검증

- Production build: PASS, compiler stderr 0 byte
- Dedicated mode-bound unit test: PASS
- Bad stride and partial I/Q/U rejection: PASS, error `-7`
- Baseline ASan low-mode reproduction: heap-buffer-overflow 확인
- Fixed ASan/UBSan low-mode and cached full-grid: findings 0
- `water_m_max=2,4,8,16,30,48,64`: finite
- Water branch: 14/14 PASS
- CCRR Chl: 11/11 PASS
- OCRT organic Chl: 10/10 PASS
- Ahn TSM: 6/6 PASS
- TSM phase-cache: 168 arrays PASS; build trace unchanged
- US62 atmosphere: PASS
- PSSA numeric fixes: 4/4 PASS
- Radiometry flux split: PASS
- PSSA layer warning: PASS

정상 `water_m_max=30` 경로의 5회 paired benchmark에서 median은
`1.73 s -> 1.70 s`로 측정되어 성능 회귀가 없었다.

## 6. Cache ON/OFF 관련 별도 상태

Shape mismatch와 OOB는 제거됐다. 다만 기존 S7/S7b LUT cache는 비-Gauss VZA를
PCHIP로 재구성하고 cache-off 경로는 zero-weight exact node를 사용한다. 따라서
비-Gauss 극단 VZA에서는 두 경로가 원래부터 bit-identical하지 않다. 대표
VZA=85 deg에서 cache ON/OFF `TOA_rho_I` 차이는 약 0.7%였고, 수정 전 기준선도
동일한 유형의 차이를 보였다.

이는 이번 mode-bound 결함과 별개의 기존 multi-view interpolation 문제이며,
진행 중인 integration-node/output-node 분리 구조에서 해소해야 한다. 이
hotfix는 cache source의 shape와 fingerprint를 정확히 맞추고, cache 사용 시
잘못된 stride 또는 OOB가 발생하지 않도록 한다.

## 7. 성능 불변조건

이번 수정은 다음을 추가하지 않는다.

```text
추가 SOS order
추가 per-order phase 계산
추가 per-order 메모리 할당
사용자 water_m_max의 강제 상향
```

S7b에서 row ring과 cache ring이 다를 때 water-to-air projection을 한 번 더
구성하지만 SOS 반복 밖에서 수행하며, 일반 단일기하 production 경로에는
적용되지 않는다.
