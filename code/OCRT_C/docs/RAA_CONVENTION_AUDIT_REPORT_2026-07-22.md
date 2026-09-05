# OCRT `rt_raa_convention.h` 주석·동작 감사 보고서

작성일: 2026-07-22

## 결론

신고 내용인 “`rt_raa_convention.h`의 최상위 RAA 주석이 실제 동작과 반대”는 현재 OCRT v1.2 기준으로 사실이 아니다.

현재 OCRT 공개 RAA의 실제 동작은 다음과 같다.

- `RAA_OCRT = 180°`: `SZA = VZA`일 때 direct rough-Fresnel glint/specular branch
- `RAA_OCRT = 0°`: 반대 principal-plane branch
- `RAA_OCRT = 90°`: `180°-RAA` 변환의 고정점으로, OCRT–OSOAA 하네스에서 convention-invariant branch

따라서 기존 최상위 주석의 0°/180° 순서는 물리 동작과 일치했다. 다만 `same/opposite half-plane`이라는 표현은 태양을 향하는 look azimuth와 하향 광자 propagation azimuth를 혼동할 수 있고, “atmosphere Fourier field uses public RAA directly”도 내부 `+pi` 변환을 생략해 설명했기 때문에 오독 가능성이 있었다.

## 코드 데이터 흐름

### 공개 입력

CLI의 `--raa` 값은 `rt_case_t::raa_deg`에 그대로 저장된다.

### 대기 Fourier 재구성

`rt_raa_to_atm_fourier_phi()`는 공개 RAA를 radian으로 전달하지만, `rt_fourier_reconstruct_*()`가 내부에서 `+pi`를 적용한다. 따라서 helper 반환값을 물리 propagation azimuth로 해석하면 안 된다.

### 대기 산란각

```text
cosTheta = -mu_s*mu_v + sin(SZA)*sin(VZA)*cos(RAA_OCRT)
```

`SZA=VZA=30°`일 때:

- RAA=0°: `cosTheta=-0.5`
- RAA=180°: `cosTheta=-1.0`

### 수중 Fourier 좌표

수중 SOS 내부 좌표는 다음을 쓴다.

```text
phi_water = pi - RAA_OCRT
```

이는 OCRT 내부 좌표이며 OSOAA 수중 출력의 `Phi`와 동일한 이름으로 해석하면 안 된다.

### 거친 Fresnel 수면

surface-local azimuth는 다음과 같다.

```text
phi_surface = RAA_OCRT - pi
```

따라서 `RAA=180°`에서 `phi_surface=0`이고 direct-glint peak가 위치한다.

## 직접 실행 검증

조건:

```text
surface = black_fresnel_ocean
wind = 3 m/s
SZA = VZA = 30°
wavelength = 865 nm
Rayleigh only
n_mu = 24
m_max = 16
```

| OCRT RAA | TOA direct-glint I |
|---:|---:|
| 0° | 8.5719628224e-09 |
| 90° | 5.8064546592e-05 |
| 180° | 3.8886180428e-01 |

`RAA=180°`가 `RAA=0°`보다 약 4.54e7배 크므로 실제 direct-glint branch는 명백히 180°다.

## 하네스와의 정합

동결 OCRT–OSOAA 하네스는 TOA에서 다음을 사용한다.

```text
Phi_OSOAA = (180° - RAA_OCRT) mod 360°
U_OCRT = -U_OSOAA
```

수중 field 비교는 다음을 사용한다.

```text
Phi_OSOAA,water = RAA_OCRT
Q/U sign change 없음
```

TOA와 수중에서 mapping이 다른 것은 각 코드·출력 레벨의 내부 azimuth와 Stokes reference-plane 정의가 다르기 때문이다. 이를 OCRT 공개 RAA 정의 자체가 바뀐 것으로 해석하면 안 된다.

## 변경 내용

물리 계산은 수정하지 않았다.

- `src/rt_raa_convention.h`
  - 운영 정의를 glint branch 기준으로 명시
  - TOA 및 수중 하네스 매핑 명시
  - atmosphere/water/surface 내부 좌표를 구분
- `src/rt_fourier.h`
  - 모호한 “forward-scattering reference at zero” 문구 제거
  - 내부 `+pi` 변환 명시
- `src/main.c`
  - `--raa` 도움말에 `180=direct-glint`, `0=opposite branch` 표시
- `tests/test_raa_convention.c`
  - 공개 RAA helper의 0°/180° 단위시험 추가
- `scripts/regression_raa_convention.sh`
  - direct-glint peak가 180°인지 실행 회귀시험 추가

## 회귀시험

- RAA helper 단위시험: PASS
- direct-glint 순서 `g180 > g90 > g0`: PASS
- direct-glint 비율 `g180 > 10^6*g0`: PASS
- 수정 전/후 RAA 0°, 90°, 180° 계산 stdout/stderr: byte-identical
- water internal-reflection msign 회귀시험: PASS
- release build: PASS

## 판단

- **코드 물리 동작 오류:** 없음
- **최상위 0°/180° 주석 반전:** 없음
- **주석의 오독 가능성:** 있음
- **필요 조치:** 설명 명확화 및 회귀시험 고정
