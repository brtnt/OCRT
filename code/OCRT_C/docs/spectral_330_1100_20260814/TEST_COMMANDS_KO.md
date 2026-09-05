# 대표 빌드·검증 명령

## C build

```bash
cd OCRT_C

gcc -std=c11 -O3 -march=native -ffp-contract=fast \
  -fassociative-math -fno-signed-zeros -fno-trapping-math \
  -DOCRT_FAST_KERNELS -fopenmp -Isrc \
  $(find src -name '*.c' | sort) -o build/ocrt -lm
```

## C data tests

```bash
bash scripts/test_pure_water_z09_table.sh
bash scripts/test_spectral_iop_330_1100.sh
bash scripts/test_detritus_candidate_gate.sh
bash scripts/test_phase_wavelength_pchip.sh
bash scripts/test_stage2_eap_disabled.sh
```

## C endpoint coupled LUT

```bash
OCRT_ADVANCED=1 OMP_NUM_THREADS=1 ./build/ocrt \
  --surface ocean --wind-speed 3 --sza 30 --pressure 1013.25 \
  --water-model ocrt --ocrt-chl 0 --ocrt-tsm 1 --ocrt-adom440 0 \
  --ocrt-tsm-species red_clay \
  --wavelength 330 \
  --n-mu-water 24 --water-m-max 4 \
  --lut-vza-step 30 --lut-vza-max 60 --lut-raa-step 90 \
  --output-full-grid out_330.csv
```

1100 nm는 `--wavelength 1100`으로 반복한다.

Authoritative replay:

```bash
OCRT_NATIVE_COUPLED_LUT_OFF=1 ...same command...
```

## C Chl range gate

350–850 nm는 active detritus를 사용한다.

```bash
OCRT_ADVANCED=1 ./build/ocrt ... \
  --water-model ocrt --ocrt-chl 1 --ocrt-tsm 0 --ocrt-adom440 0 \
  --wavelength 443
```

330 또는 1100 nm에서 같은 Chl 조건은 명시적으로 실패해야 한다.

## Python tests

```bash
python -m compileall -q .
python -m pytest -q \
  tests/test_spectral_contract_330_1100.py \
  tests/test_spectral_data_hygiene.py \
  tests/test_native_coupled_lut.py \
  tests/test_fullgrid_aerosol_object_fix.py \
  tests/test_rrs_iqu_output.py \
  tests/test_water_intrefl_msign.py \
  tests/test_dtpsign_and_component_diagnostic.py \
  tests/test_pure_water_z09_table.py \
  tests/test_stage2_eap_disabled.py
```

## Data manifest

```bash
python build_installation_manifest.py
python validate_new_payload.py
```
