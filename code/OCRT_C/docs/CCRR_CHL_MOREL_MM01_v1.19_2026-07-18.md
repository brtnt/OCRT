# OCRT v1.19 CCRR Chl 흡광자료 통합

## 1. 목적

v1.18의 CCRR 구성성분 어댑터는 양수 `--ccrr-chl`에서 5열 계수표를 요구했지만
배포 파일이 없어 실행되지 않았다. v1.19는 사용자 제공 정규화 스펙트럼을 현재
OCRT 인터페이스에 맞게 변환하여 이 경로를 완결한다.

RT 적분기는 변경하지 않는다. CCRR 분기는 계속 구성성분 입력을 OCRT의 IOP와
입자 위상 표현으로 변환하고, 방사전달은 OCRT vector SOS solver가 수행한다.

## 2. 원자료와 모델 식

원자료는 440 nm에서 1로 정규화된 무차원 스펙트럼 `A_chl(lambda)`이다. 파일
헤더는 Morel (1988, Fig. 10c; Prieur and Sathyendranath 1981)의 스펙트럼임을
명시한다. 300–350 nm와 700–1000 nm는 원자료 설명대로 외삽 구간이다.

CCRR Chl 흡광은 HydroLight classic Case-1 / Morel-Maritorena 형태로 정의한다.

```text
a_p(lambda) = 0.06 * A_chl(lambda) * Chl^0.65    [m^-1]
```

`Chl` 단위는 `mg m^-3`이다.

## 3. OCRT 5열 인터페이스 변환

현재 로더는 다음 다섯 열을 읽고 마지막 두 열만 사용한다.

```text
wavelength_nm  Ap  Ep  Aphi  Ephi
```

따라서 변환은 다음과 같다.

```text
Ap    = 0                 # 현재 로더에서 미사용
Ep    = 0                 # 현재 로더에서 미사용
Aphi  = 0.06 * A_chl
Ephi  = 0.65
```

대표값:

| 파장 | 정규화 `A_chl` | `Aphi` | `Ephi` |
|---:|---:|---:|---:|
| 440 nm | 1.0000 | 0.060000 | 0.65 |
| 555 nm | 0.2071 | 0.012426 | 0.65 |
| 670 nm | 0.5930 | 0.035580 | 0.65 |

443 nm는 440–445 nm 사이에서 선형 보간되며 `Aphi(443)=0.056382`이다.

## 4. 파일

Canonical 파일:

```text
inputs/water_iop/aph_ccrr_morel1988_mm01.txt
```

v1.18 이하 로더 호환 사본:

```text
inputs/water_iop/aph_bricaud_1998.txt
```

두 파일은 byte-identical이다. 후자의 이름은 역사적 호환성만을 위한 것이며,
내용을 Bricaud et al. (1998)의 파장별 `A(lambda), E(lambda)` 표로 해석해서는 안 된다.

원자료 보존본:

```text
inputs/water_iop/source/apstarchl_morel1988_normalized_user_supplied.txt
```

재생성 도구:

```text
scripts/build_ccrr_chl_table.py
```

## 5. 로더 정책

v1.19 로더는 canonical 파일을 먼저 찾고, 없을 때만 역사적 파일명으로 fallback한다.
각 행의 파장 단조 증가, 양의 파장, 비음수 `Aphi`, 유한 `Ephi`를 검사한다.
파장 보간은 기존과 동일한 선형 보간이고, 표 범위 밖에서는 nearest-edge 값이다.

## 6. 사용자 인터페이스

```bash
./build/v2_solver_vk_v1.19 \
  --surface ocean --water-model ccrr --wind-speed 0 \
  --ccrr-chl 0.3 --ccrr-tsm 0 --ccrr-adom440 0 \
  --sza 30 --vza 20 --raa 90 --wavelength 443 --pressure 0
```

`--ccrr-chl`, `--ccrr-tsm`, `--ccrr-adom440`은 계속 모두 명시해야 한다. 세 값이
모두 0이면 native pure-water 경로로 하강한다.

## 7. 과학적 명칭 주의

사용자 제공 파일명에는 `bricaud_2011`이 포함되어 있으나 파일 내부 설명과 실제
스펙트럼은 Morel (1988) 정규화 형상이다. 따라서 v1.19 문서와 canonical 파일명은
`Morel1988/MM01 classic Case-1`로 기록한다. Bricaud 1998 비선형 `A(lambda),
E(lambda)` 자료를 새로 추가한 것으로 기술하지 않는다.
