# OCRT v1.19 해수 광특성 입력 분기 명세

## 1. 최상위 계약

`--surface ocean`에서는 다음 중 하나를 정확히 한 번 지정한다.

```text
--water-model ocrt
--water-model ccrr
--water-model iop
```

선택하지 않거나, 두 번 이상 지정하거나, 다른 분기의 접두사 옵션을 혼용하면
종료코드 2로 실패한다. 세 분기 모두 방사전달 적분은 OCRT vector SOS solver가
수행한다.

## 2. OCRT 분기

필수:

```text
--ocrt-chl X
--ocrt-tsm X
--ocrt-adom440 X
```

선택:

```text
--ocrt-adom-slope X
--ocrt-phyto-group pico|nano|micro
--ocrt-tsm-species red_clay|brown_earth|yellow_clay|calcareous_sand
--ocrt-detritus-a440 X
--ocrt-detritus-slope X
```

세 필수 값이 모두 0이면 native pure water로 하강한다.

## 3. CCRR 분기

필수:

```text
--ccrr-chl X
--ccrr-tsm X
--ccrr-adom440 X
```

선택:

```text
--ccrr-adom-slope X
--ccrr-phase-moments F
```

### 3.1 Chl 변환

```text
a_p(lambda) = 0.06 * A_chl(lambda) * Chl^0.65
```

`A_chl(lambda)`는 440 nm에서 1인 Morel (1988) 정규화 스펙트럼이다. 현 5열
로더에는 다음과 같이 표현한다.

```text
Aphi(lambda) = 0.06 * A_chl(lambda)
Ephi(lambda) = 0.65
```

자료 파일은 `inputs/water_iop/aph_ccrr_morel1988_mm01.txt`이다. 역사적
`aph_bricaud_1998.txt`는 fallback 호환 사본이며 과학적 출처명을 뜻하지 않는다.

### 3.2 TSM과 aDOM

CCRR TSM 경험식과 공통 aDOM 지수식은 기존 동작을 유지한다. aDOM slope를
생략하면 `0.014 nm^-1`이다.

### 3.3 순수해수

`--ccrr-chl 0 --ccrr-tsm 0 --ccrr-adom440 0`은 OCRT 0/0/0과 동일한 native
pure-water 경로를 사용한다.

## 4. IOP 분기

필수:

```text
--iop-a A       # A >= 0
--iop-b B       # B > 0
--iop-bb BB     # 0 < BB < 0.5 B
```

직접 phase source는 최대 하나만 선택한다.

## 5. 공통 순수해수 입력

```text
--water-temperature C
--water-salinity G
--wind-speed W
```

일반 사용자는 순수해수 모델 자체를 선택하지 않는다. 저수준 모델·수치 옵션은
기존 `OCRT_ADVANCED=1` 또는 `OCRT_DEBUG=1` 게이트를 유지한다.
