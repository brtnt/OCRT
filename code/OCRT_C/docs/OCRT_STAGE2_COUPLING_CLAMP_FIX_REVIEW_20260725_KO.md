# OCRT 2단계 water-to-air coupling 64절점 silent-clamp 수정 및 검증 보고서

문서 기준일: 2026-07-25  
수정 전 기준선: `OCRT-v1.2-2026-07-24-KST-stage2-interface-fix123-dtp-uclosure-pole-limit`  
수정 후 기준선: `OCRT-v1.2-2026-07-25-KST-stage2-coupling-clamp-fix`

## 1. 결론

`src/rt_air_water_coupling.c`의 water-to-air 역결합 보간에서 입력 방향 수를 64개로 조용히 절단하던 결함을 제거하였다. 이 결함은 수중 solver가 양의 Gauss 절점 뒤에 추가하는 zero-weight 관측·태양·천저 방향을 버려, 실제로 필요한 보간을 외삽으로 바꾸었다. LUT 생산 기본값인 `n_mu_water=96`은 결함 발동범위에 직접 포함된다.

수정은 다음 범위로 제한하였다.

- 64개 고정 배열 및 `if (n > 64) n = 64` 제거
- 검증된 solver 상한을 포괄하는 320개 stack workspace 사용
- 미래의 더 큰 격자에는 heap fallback 적용
- 할당 실패 시 조용한 절단 대신 즉시 실패
- 분위기·대기 하네스·FIX1–FIX4·exact-pole·in-water SOS는 변경하지 않음

본 수정 후 `n_mu_water=64/96`에서 천저 Q의 금지된 Fourier `m=0` 성분이 0으로 복원되었다. `rrs(0−)`는 수정 전후 완전 불변이며, 결함이 in-water solver 이후의 water-to-air coupling 단계에 국소화됨을 확인하였다.

## 2. 결함의 코드 구조

기존 함수는 다음과 같이 local sorting workspace를 64개로 고정하였다.

```c
int ord[64];
double ms[64], tmp[64];
if (n > 64) n = 64;
```

그러나 함수에 전달되는 `n`은 사용자가 요청한 기본 Gauss node 수가 아니라 관측·태양·천저 등 zero-weight 특수방향을 포함한 filled grid 크기이다. 배열 뒤쪽의 특수방향부터 잘리므로, 천저 및 지정 view 방향의 정확한 node가 보간 표본에서 사라진다.

수정 후 구현은 다음 계약을 사용한다.

```c
enum { AW_INTERP_STACK_N = 320 };

/* n <= 320: stack workspace
 * n > 320 : heap workspace
 * allocation failure: abort, never silent truncation */
```

이미 정렬된 64개 이하 경로의 산술순서는 유지된다.

## 3. 모듈 단위 검증

### 3.1 source audit

다음을 자동 검사하였다.

- `int ord[64]` 부재
- `double ms[64]` 부재
- `if (n > 64) n = 64` 부재
- `AW_INTERP_STACK_N = 320` 존재
- heap fallback 및 fail-loud `abort()` 존재

### 3.2 ASAN/UBSAN workspace 시험

48, 64, 96, 320, 321, 513개 방향에 대해 비정렬 입력, 보간 및 양측 외삽을 검사하였다. 320 이하 stack 경로와 321 이상 heap 경로가 모두 해석적 선형함수를 재현하였다.

```text
PASS: air-water interpolation workspace stack/heap paths
ASAN/UBSAN error: none
```

### 3.3 기존 동결 모듈 비회귀

다음 자동회귀를 재실행하고 모두 통과하였다.

- public water RAA source audit
- FIX1+FIX2+FIX3 interface analytic regression
- FIX4 diffuse-top source/operator closure
  - pure Rayleigh: `2.776e-17`
  - Rayleigh–aerosol: `5.551e-17`
- exact-pole continuity and spin-2 covariance
- direct-glint RAA convention
- external-bottom mode bounds

따라서 이미 정합성이 확보된 물리모듈에는 재진입하지 않았다.

## 4. 전용 full-grid 게이트

조건은 443 nm, SZA 40°, VZA 0/30/60°, RAA 45° 간격, 24개 full-grid 셀이다. 대기 설정은 기존 동결조건을 그대로 유지하였다.

### 4.1 `n_mu_water=48`

수정 전·후 full-grid CSV가 byte-identical이다. 64절점 clamp가 발동하지 않는 기존 72실행·600셀 검증영역은 완전 비회귀이다.

### 4.2 `n_mu_water=64`

- `rrs(0−)` I/Q/U: 완전 불변
- `Rrs(0+)` 및 TOA: 잘못된 역결합 외삽만 교정
- 천저 `Rrs_Q`의 `|m0|/|m2|`: `11.4118% → 0%`
- 천저 `TOA_Q`의 `|m0|/|m2|`: `1.3891% → 0%`

### 4.3 `n_mu_water=96`

LUT 생산 기본 해상도에서 다음 변화가 확인되었다.

- 천저 `Rrs_Q`의 `|m0|/|m2|`: `107.8935% → 0%`
- 천저 `TOA_Q`의 `|m0|/|m2|`: `7.4066% → 0%`
- `rrs(0−)` I/Q/U: 완전 불변
- `Rrs I` 최대 변화 / 수정 전 성분최대: `9.568%`
- `Rrs Q`: `30.066%`
- `Rrs U`: `35.718%`

이 변화량은 OCRT–OSOAA 정확도 수치가 아니라 기존 silent clamp가 생산 LUT 출력을 왜곡한 크기이다.

## 5. 단일기하와 full-grid 일관성

`n_mu_water=64`, VZA 30°, RAA 45°에서 single-target와 full-grid 결과를 비교하였다. view-list에 따라 추가되는 zero-weight node 집합이 달라 TOA에는 약 `1e-6` 수준의 기존 샘플링 차이가 존재한다. 반면 Rrs/rrs는 `2e-8` 이내로 일치하였다. 이 차이는 clamp 수정의 실패로 분류하지 않으며, 전용 회귀에서는 “bit parity”가 아니라 정해진 tolerance의 일관성으로 관리한다.

## 6. 72실행·600셀 비회귀

기존 동결 manifest를 변경하지 않고 72개 물리 실행을 모두 재실행하였다.

```text
success                  = 72/72
expected row count       = 72/72
water convergence        = 72/72
byte-identical CSV       = 72/72
maximum water SOS order  = 29
```

600셀 병합표의 57개 수치열은 pole-limit 기준표와 최대 절대차 `0`이다. `OCRT_source_version`만 새 버전으로 갱신하였다. 이 표는 최신 corrected-depth OSOAA 기준 이전의 자료이므로 물리 정확도 최종 승인용이 아니라 코드 비회귀용으로만 사용한다.

## 7. 계산시간

동일 compiler, fixed ISA, `OMP_NUM_THREADS=1`, 443 nm 24셀 full-grid에서 1회 예열 후 variant별 3회 교차 측정하였다.

```text
n_mu=48  mean -1.765%  median -1.974%
n_mu=64  mean -7.751%  median -1.122%
n_mu=96  mean +0.638%  median +0.254%
```

64절점 mean은 수정 전 1회의 10.55초 outlier 때문에 음의 변화가 크게 나타났으며 중앙값 변화는 -1.12%이다. 생산 기본 96절점의 mean 변화는 +0.64%, 중앙값은 +0.25%이다. 5% 병목조사 기준에 미달하므로 추가 프로파일링은 필요하지 않다.

## 8. C/Python 영향범위

현재 배포 트리에는 이 C water-to-air coupling 보간의 독립 Python 물리 구현이 없다. Python 파일은 실행 하네스·분석도구이며 동일한 `aw_interp_on_unsorted()` 경로를 별도로 구현하지 않는다. 따라서 C 물리코드만 변경하였고, Python 쪽에는 포팅할 대응 구현이 없다.

## 9. 최종 승인 상태

### 동결 유지

- 1단계 대기 모듈과 대기 하네스
- public water RAA/U reconstruction
- FIX1+FIX2+FIX3
- FIX4 diffuse-top incoming-U closure
- rough-interface exact-pole limit
- in-water SOS

### 이번 단계에서 추가 동결

- water-to-air coupling interpolation의 전체 filled-grid 보존
- 320 stack + fail-loud heap fallback 계약
- `n_mu_water=48/64/96` 회귀
- 천저 Q forbidden-`m0` 회귀

### 계속 미해결

- corrected-depth OSOAA SZA 80° `rrs(0−)`
- corrected-depth OSOAA `Rrs(0+)`의 두 extraction branch
- 동일 수심·층경계·`n_mu`·`m_max`·SOS tolerance 최종 게이트
- rough direct-solar boundary와 in-water solver의 분리 폐합
- `R_ww` 전 성분·전 mode 감사

## 10. 최종 버전

```text
OCRT-v1.2-2026-07-25-KST-stage2-coupling-clamp-fix
```

승인 바이너리 SHA-256:

```text
328f28b8f50b6d98207eabe57309dde3538e4658a1ceec8f53d970b9dda3587a
```
