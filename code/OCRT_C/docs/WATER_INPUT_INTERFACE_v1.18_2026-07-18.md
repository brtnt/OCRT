# OCRT v1.18 해수 광특성 입력 분기 명세

## 1. 적용 범위

이 문서는 `--surface ocean`에서 사용하는 해수 광특성 입력 계약을 정의한다.
수중 방사전달 계산기는 세 분기 모두 **기존 OCRT vector SOS RT solver**를 사용한다.
분기 차이는 OCRT solver에 전달할 IOP와 입자 위상자료를 만드는 방법뿐이다.

## 2. 최상위 분기 선택

해양 수면 실행은 다음 중 하나를 **정확히 한 번** 선택해야 한다.

```text
--water-model ocrt
--water-model ccrr
--water-model iop
```

규칙:

1. `--surface ocean`인데 `--water-model`이 없으면 종료코드 2로 실패한다.
2. 같은 값이든 다른 값이든 `--water-model`을 두 번 이상 지정하면 실패한다.
3. 선택한 분기와 다른 접두사의 옵션을 함께 쓰면 실패한다.
4. `black`, `flat`, `coxmunk` surface에서는 water branch 옵션을 사용할 수 없다.
5. 옵션 순서는 의미가 없다. 파서가 입력 출처를 기록한 뒤 한 번에 검증한다.

## 3. OCRT 구성성분 분기

### 3.1 필수 입력

다음 세 값을 모두 명시해야 한다. 값 `0`은 유효하다.

```text
--ocrt-chl X       # chlorophyll-a, mg m^-3
--ocrt-tsm X       # inorganic TSM dry mass, g m^-3
--ocrt-adom440 X   # aDOM absorption at 440 nm, m^-1
```

선택 입력:

```text
--ocrt-adom-slope X       # default 0.014 nm^-1
--ocrt-phyto-group G      # pico | nano | micro; default micro
--ocrt-tsm-species S      # red_clay | brown_earth | yellow_clay | calcareous_sand
--ocrt-detritus-a440 X    # default 0 m^-1
--ocrt-detritus-slope X   # default 0.0109 nm^-1
```

`--ocrt-phyto-group`와 `--ocrt-detritus-*`는 `--ocrt-chl > 0`일 때만 유효하다.
`--ocrt-tsm-species`는 `--ocrt-tsm > 0`일 때만 유효하다.

### 3.2 순수해수 하강 규칙

다음 조건이면 구성성분 변환기를 실행하지 않고 native pure-water 경로로 하강한다.

```text
--water-model ocrt
--ocrt-chl 0
--ocrt-tsm 0
--ocrt-adom440 0
```

즉, `a_w`, `b_w`, 순수해수 vector Rayleigh-like phase만 사용한다.
온도와 염분은 공통 순수해수 입력인 `--water-temperature`,
`--water-salinity`로 지정한다.

### 3.3 OCRT 입자 모델

양의 Chl는 EAP phytoplankton과 Chl-covarying organic detritus를 생성한다.
양의 TSM은 Ahn 4종 광물자료를 사용한다. 각 입자의 P11/P12/P33은 산란계수로
가중 혼합된 뒤 OCRT vector solver로 전달된다.

## 4. CCRR 구성성분 어댑터 분기

### 4.1 필수 입력

```text
--ccrr-chl X       # mg m^-3
--ccrr-tsm X       # g m^-3
--ccrr-adom440 X   # m^-1 at 440 nm
```

선택 입력:

```text
--ccrr-adom-slope X       # default 0.014 nm^-1
--ccrr-phase-moments F    # optional CCRR particle-moment override
```

CCRR은 별도의 RT solver가 아니다. 위 입력을 OCRT solver가 소비하는 IOP/phase
representation으로 변환하는 전처리 어댑터다.

### 4.2 순수해수 하강 규칙

다음 조건도 OCRT 분기의 0/0/0과 동일한 native pure-water 경로를 사용한다.

```text
--water-model ccrr
--ccrr-chl 0
--ccrr-tsm 0
--ccrr-adom440 0
```

OCRT 0/0/0과 CCRR 0/0/0은 동일한 실행 플래그와 동일한 수치 solver path로
하강한다. 분기 선택은 provenance로만 유지된다.

### 4.3 현재 패키지의 CCRR Chl 자료 제약

CCRR의 양수 Chl 변환 함수는 존재하지만 다음 계수표를 런타임에 요구한다.

```text
inputs/water_iop/aph_bricaud_1998.txt
```

이 파일은 현재 패키지에 포함되어 있지 않다. 따라서 `--ccrr-chl > 0`은 자료 로드
오류로 종료한다. `--ccrr-chl 0`인 CCRR TSM/aDOM 경로와 CCRR 0/0/0 순수해수
하강 경로는 정상 동작한다. 이 제약은 이번 CLI 분기 정리에서 수치자료를 임의로
추가하지 않고 명시적으로 보존했다.

## 5. 직접 IOP 분기

총 bulk IOP를 직접 입력한다.

```text
--water-model iop
--iop-a A
--iop-b B
--iop-bb BB
```

검증식:

```text
A >= 0
B > 0
0 < BB < 0.5 B
```

선택 phase source는 최대 하나만 지정할 수 있다.

```text
--iop-phase-lut F
--iop-mie-phase F
```

phase source를 생략하면 입력 `BB/B`에 맞는 기존 analytic scalar closure를 쓴다.
`--iop-mie-*` 변환 옵션은 `--iop-mie-phase`가 있을 때만 허용한다.

## 6. 공통 순수해수 입력과 advanced 입력

다음은 OCRT/CCRR/IOP 분기 접두사를 사용하지 않는 공통 매질·solver 입력이다.

```text
--water-temperature C
--water-salinity G
--wind-speed W
```

순수해수 모델 자체 또는 수치 변환을 바꾸는 저수준 옵션은 기존
`OCRT_ADVANCED=1`/`OCRT_DEBUG=1` 게이트를 유지한다. 일반 사용자는 농도와 온도,
염분 이외의 순수해수 모델 선택을 지정할 필요가 없다.

## 7. 폐기된 모호한 옵션

다음 계열은 자동 alias로 변환하지 않고 오류로 종료한다.

```text
--simple-*
--simple-mode
--ccrr-mode
--chl --tsm --adom440 --adom-slope
--phyto-group --tsm-species --detritus-*
--cdom-a440 --cdom-slope --cdom-ref-lambda
--fixed-bulk-iop
--water-mie-*
```

자동 alias를 두지 않은 이유는 잘못된 분기에 입력이 조용히 들어가는 것을 막기
위함이다. 오류 메시지는 대응하는 `--ocrt-*`, `--ccrr-*`, `--iop-*` 이름을 안내한다.

## 8. 배치 full-grid CSV

CLI base에서도 `--water-model`과 해당 필수 입력을 지정해야 한다. 각 행은 base를
상속하며, `water_model`을 바꾸는 행은 새 분기의 필수 값을 모두 제공해야 한다.

지원되는 해수 열:

```text
water_model
ocrt_chl, ocrt_tsm, ocrt_adom440, ocrt_adom_slope
ocrt_phyto_group, ocrt_tsm_species
ocrt_detritus_a440, ocrt_detritus_slope
ccrr_chl, ccrr_tsm, ccrr_adom440, ccrr_adom_slope
iop_a, iop_b, iop_bb, iop_mie_phase
```

한 행에서 서로 다른 분기의 접두사를 섞으면 해당 batch 전체를 입력 오류로 종료한다.

## 9. 내부 lowering

public mode는 `rt_water_input_mode_t` 하나로 보존한다.

```text
RT_WATER_INPUT_OCRT
RT_WATER_INPUT_CCRR
RT_WATER_INPUT_IOP
```

검증 완료 후에만 기존 실행 필드로 lowering한다.

```text
OCRT/CCRR nonzero -> ccrr_mode=1 + constituent model enum
OCRT/CCRR zero    -> ccrr_mode=0 native pure-water path
IOP               -> fixed_bulk_iop_mode=1
```

이 구조는 사용자 입력 계약과 기존 검증된 OCRT 수치 kernel을 분리한다.
