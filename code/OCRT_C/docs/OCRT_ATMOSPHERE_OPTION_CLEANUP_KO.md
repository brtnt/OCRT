# 대기 옵션 정리

## 권장 이름

| 권장 옵션 | 기존 호환 별칭 | 의미 |
|---|---|---|
| `--surface-pressure P` | `--pressure P` | 표면기압 hPa |
| `--gas-profile NAME` | `--atm-profile NAME` | AFGL 흡광성 가스 연직 프로파일 |
| `--aod-band X` | `--aod X` | 현재 계산 파장의 AOD |
| `--rayleigh on|off` | `--no-rayleigh` | bulk molecular scattering 활성 여부 |

기존 스크립트 보존을 위해 별칭은 유지한다. 새 문서와 예제는 권장 이름을 사용한다.

## 물리적 역할 분리

- molecular vertical profile: US Standard Atmosphere 1962, 내부 고정.
- Rayleigh column: Bodhaine 1999, `--surface-pressure`로 scaling.
- absorbing gases: AFGL profile, `--gas-profile`과 gas-column 옵션.
- aerosol: `.mie` + `--aod-band` + 2 km exponential profile.
- PSSA: `--pssa`, 태양 직달 하향경로만 구면 보정.
