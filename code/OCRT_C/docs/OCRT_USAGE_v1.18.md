# OCRT v1.18 사용 가이드

## 빌드

```bash
cd ocrt
gcc -std=c11 -O3 -march=native -ffp-contract=fast -fassociative-math \
  -fno-signed-zeros -fno-trapping-math -DOCRT_FAST_KERNELS -fopenmp \
  -Isrc $(find src -name '*.c' -print | LC_ALL=C sort) \
  -o build/v2_solver_vk_v1.18 -lm
```

실행은 `ocrt/` 디렉터리에서 한다. bitwise 회귀시험은
`OMP_NUM_THREADS=1`을 기준으로 한다.

## 해수 입력의 핵심 규칙

`--surface ocean`에서는 아래 세 분기 중 정확히 하나를 선택한다.

```text
--water-model ocrt | ccrr | iop
```

분기 미선택, 중복 선택, 접두사 혼용은 모두 종료코드 2다.

## OCRT 예시

```bash
./build/v2_solver_vk_v1.18 \
  --surface ocean --water-model ocrt --wind-speed 0 \
  --ocrt-chl 0.3 --ocrt-tsm 1.0 --ocrt-adom440 0.05 \
  --ocrt-adom-slope 0.014 \
  --ocrt-phyto-group micro --ocrt-tsm-species red_clay \
  --sza 30 --vza 20 --raa 90 --wavelength 443 --pressure 0
```

`--ocrt-adom-slope`를 생략하면 `0.014 nm^-1`이다. TSM species 기본값은
`red_clay`, phytoplankton group 기본값은 `micro`다.

## CCRR 예시

```bash
./build/v2_solver_vk_v1.18 \
  --surface ocean --water-model ccrr --wind-speed 0 \
  --ccrr-chl 0 --ccrr-tsm 1.0 --ccrr-adom440 0.05 \
  --sza 30 --vza 20 --raa 90 --wavelength 443 --pressure 0
```

CCRR은 구성성분을 OCRT RT 입력으로 변환하는 어댑터다. RT 적분은 OCRT가 한다.

현재 패키지에는 CCRR 양수 Chl 경로가 요구하는
`inputs/water_iop/aph_bricaud_1998.txt`가 없다. 따라서 `--ccrr-chl > 0`은 아직
실행되지 않으며, 위 예시는 `Chl=0`인 CCRR TSM/aDOM 경로를 사용한다.

## 순수해수 예시

OCRT 선택을 유지하면서 순수해수만 계산:

```bash
./build/v2_solver_vk_v1.18 \
  --surface ocean --water-model ocrt --wind-speed 0 \
  --ocrt-chl 0 --ocrt-tsm 0 --ocrt-adom440 0 \
  --water-temperature 20 --water-salinity 38.4 \
  --sza 30 --vza 20 --raa 90 --wavelength 443 --pressure 0
```

CCRR에서도 `--ccrr-chl 0 --ccrr-tsm 0 --ccrr-adom440 0`을 주면 동일한
native pure-water path를 실행한다.

## 직접 IOP 예시

```bash
./build/v2_solver_vk_v1.18 \
  --surface ocean --water-model iop --wind-speed 0 \
  --iop-a 0.1762 --iop-b 2.618309 --iop-bb 0.038578 \
  --sza 30 --vza 20 --raa 90 --wavelength 443 --pressure 0
```

vector `.mie` phase를 직접 지정할 때만 다음을 추가한다.

```text
--iop-mie-phase path/to/phase.mie
```

## 대표 입력 오류

```text
--surface ocean, no --water-model         -> error
--water-model ocrt --water-model ccrr    -> error
--water-model ocrt + --ccrr-*            -> error
--water-model ocrt without one of the three required values -> error
--water-model iop without a/b/bb         -> error
```

## 회귀시험

```bash
scripts/smoke_cli_water_branch_v118.sh build/v2_solver_vk_v1.18
scripts/smoke_cli_organic_chl.sh build/v2_solver_vk_v1.18
scripts/smoke_cli_tsm_ahn.sh build/v2_solver_vk_v1.18
scripts/regression_tsm_phase_cache_v116.sh build/v2_solver_vk_v1.18
```

상세 계약은 `docs/WATER_INPUT_INTERFACE_v1.18_2026-07-18.md`를 참조한다.
