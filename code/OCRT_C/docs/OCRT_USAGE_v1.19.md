# OCRT v1.19 사용 가이드

## 빌드

```bash
cd ocrt
gcc -std=c11 -O3 -march=native -ffp-contract=fast -fassociative-math \
  -fno-signed-zeros -fno-trapping-math -DOCRT_FAST_KERNELS -fopenmp \
  -Isrc $(find src -name '*.c' -print | LC_ALL=C sort) \
  -o build/v2_solver_vk_v1.19 -lm
```

실행은 `ocrt/` 디렉터리에서 한다. bitwise 회귀시험은 `OMP_NUM_THREADS=1`을
기준으로 한다.

## 해수 입력 분기

`--surface ocean`에서는 아래 셋 중 정확히 하나를 한 번만 선택한다.

```text
--water-model ocrt | ccrr | iop
```

분기 미선택, 중복 선택, 다른 분기 접두사 혼용은 종료코드 2다.

## OCRT 구성성분 모델

```bash
./build/v2_solver_vk_v1.19 \
  --surface ocean --water-model ocrt --wind-speed 0 \
  --ocrt-chl 0.3 --ocrt-tsm 1.0 --ocrt-adom440 0.05 \
  --ocrt-phyto-group micro --ocrt-tsm-species red_clay \
  --sza 30 --vza 20 --raa 90 --wavelength 443 --pressure 0
```

필수 입력은 `--ocrt-chl`, `--ocrt-tsm`, `--ocrt-adom440`이다. 값 0은 유효하다.

## CCRR 구성성분 어댑터

```bash
./build/v2_solver_vk_v1.19 \
  --surface ocean --water-model ccrr --wind-speed 0 \
  --ccrr-chl 0.3 --ccrr-tsm 1.0 --ccrr-adom440 0.05 \
  --ccrr-adom-slope 0.014 \
  --sza 30 --vza 20 --raa 90 --wavelength 443 --pressure 0
```

필수 입력은 `--ccrr-chl`, `--ccrr-tsm`, `--ccrr-adom440`이다. CCRR은 별도 RT
solver가 아니라 OCRT 입력 변환 어댑터다.

양수 Chl 흡광은 bundled Morel (1988) 정규화 스펙트럼과 classic Case-1 폐쇄를
사용한다.

```text
a_p(lambda) = 0.06 * A_chl(lambda) * Chl^0.65
```

canonical 자료 파일:

```text
inputs/water_iop/aph_ccrr_morel1988_mm01.txt
```

300–350 nm와 700–1000 nm는 원자료의 외삽 구간이므로 해석에 주의한다.

## 순수해수

OCRT 또는 CCRR에서 세 구성성분 입력을 모두 명시적으로 0으로 주면 동일한 native
pure-water 경로를 실행한다.

```bash
./build/v2_solver_vk_v1.19 \
  --surface ocean --water-model ccrr --wind-speed 0 \
  --ccrr-chl 0 --ccrr-tsm 0 --ccrr-adom440 0 \
  --water-temperature 20 --water-salinity 38.4 \
  --sza 30 --vza 20 --raa 90 --wavelength 443 --pressure 0
```

## 직접 IOP

```bash
./build/v2_solver_vk_v1.19 \
  --surface ocean --water-model iop --wind-speed 0 \
  --iop-a 0.1762 --iop-b 2.618309 --iop-bb 0.038578 \
  --sza 30 --vza 20 --raa 90 --wavelength 443 --pressure 0
```

## 회귀시험

```bash
scripts/smoke_cli_water_branch_v118.sh build/v2_solver_vk_v1.19
scripts/smoke_cli_ccrr_chl_v119.sh build/v2_solver_vk_v1.19
scripts/smoke_cli_organic_chl.sh build/v2_solver_vk_v1.19
scripts/smoke_cli_tsm_ahn.sh build/v2_solver_vk_v1.19
scripts/regression_tsm_phase_cache_v116.sh build/v2_solver_vk_v1.19
```

세부 자료 변환은 `docs/CCRR_CHL_MOREL_MM01_v1.19_2026-07-18.md`를 참조한다.
