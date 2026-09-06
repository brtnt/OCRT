# SnF aerosol Mie generator

This directory holds the reproducible generator used for the 2026-07-19 integration of the SnF
atmospheric aerosol models. It implements the OCRT V1 / 6SV-lineage Bohren–Huffman Mie
calculation, log-normal radius integration, component mixing and the OPAC/OCRT P11/P12/P33
writer.

## Canonical output

- wavelengths: the 20 OPAC wavelengths, 0.350–3.750 µm;
- angles: 180° to 0° at 0.5° spacing, 361 points;
- radius step: `rlogpas = 0.011`;
- output blocks: P11, P12 and P33;
- model set: T50/T80/T90/T95, C50/C70/C80/C90/C95, M50C/M70C/M80C/M90C/M95C/M98C, O99.

`O99` is the canonical model name supplied by the source archive; no `O99C` input exists.

## Build and generate

Run from `ocrt/`:

```bash
make -C tools/snf_mie_generator
OMP_NUM_THREADS=32 ./scripts/generate_snf_aerosols.sh /tmp/snf_mie
```

The source `.inp` files are kept under `inputs/aerosol_snf_inp/`.

## Runtime smoke test

```bash
python3 tools/snf_mie_generator/run_snf_smoke.py \
  --exe ./build/ocrt \
  --mie-dir /tmp/snf_mie \
  --out /tmp/snf_smoke.csv \
  --workers 16
```

## Controlled regression result

The pre-existing high-resolution T50, C50, M80C and O99 tables were regenerated at file
precision. Spectral values, P11, P12 away from 0°/180° and P33 were exact; the only differences
were floating-point residuals below 2.3e-15 at the algebraic-zero P12 end points. The runtime
package therefore keeps these four pre-existing files byte for byte. The other twelve canonical
tables are direct generator outputs.

The full validation archive is under `validation/snf_aerosol_2026-07-19/`.

---

# SnF 에어로졸 Mie 생성기

이 디렉터리는 2026-07-19 SnF 대기 에어로졸 모델 통합에 쓴 재현 가능한 생성기를 담는다. OCRT V1 / 6SV
계보의 Bohren–Huffman Mie 계산, 로그정규 반지름 적분, 성분 혼합, OPAC/OCRT P11/P12/P33 출력기를
구현한다.

## 표준 출력

- 파장: OPAC 파장 20개, 0.350–3.750 µm;
- 각도: 180° 에서 0° 까지 0.5° 간격, 361점;
- 반지름 간격: `rlogpas = 0.011`;
- 출력 블록: P11, P12, P33;
- 모델 집합: T50/T80/T90/T95, C50/C70/C80/C90/C95, M50C/M70C/M80C/M90C/M95C/M98C, O99.

`O99` 는 원본 아카이브가 제공한 표준 모델 이름이다. `O99C` 입력은 없다.

## 빌드와 생성

`ocrt/` 에서 실행한다.

```bash
make -C tools/snf_mie_generator
OMP_NUM_THREADS=32 ./scripts/generate_snf_aerosols.sh /tmp/snf_mie
```

원본 `.inp` 파일은 `inputs/aerosol_snf_inp/` 에 보존한다.

## 실행 스모크 테스트

```bash
python3 tools/snf_mie_generator/run_snf_smoke.py \
  --exe ./build/ocrt \
  --mie-dir /tmp/snf_mie \
  --out /tmp/snf_smoke.csv \
  --workers 16
```

## 통제 회귀 결과

기존 고해상도 T50, C50, M80C, O99 표를 파일 정밀도로 재생성했다. 분광값, P11, 0°/180° 를 제외한 P12, P33 은
정확히 같았고, 유일한 차이는 P12 의 대수적 0 끝점에서 2.3e-15 미만의 부동소수 잔차였다. 따라서 실행
패키지는 이 네 파일을 바이트 그대로 유지한다. 나머지 표준 표 12개는 생성기 직접 출력이다.

전체 검증 아카이브는 `validation/snf_aerosol_2026-07-19/` 에 있다.
