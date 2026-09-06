# PSSA numeric fixes — validation assets

Date: 2026-07-19

## Scope

The matrices in this directory quantify four changes relative to the preceding
`pssa-aodref-radiometry-split-layerwarn` package:

1. sign of the reflected-beam lower-end-point attenuation,
2. gas-only ocean PSSA direct-transmittance helper,
3. coupling of the reflected downward first order into the higher SOS orders,
4. layer-local reflected-beam `beta` phase source.

## Source sequence

```text
b92a476 baseline package
f07a0e7 endpoint sign
4a831e5 gas-only ocean helper
8a6a380 downward first-order coupling
ff46823 layer-local beta source
f284cb4 direct spherical reflected amplitude
27134f1 pure-Rayleigh beta setup optimization
```

Intermediate binaries were used to produce the differential matrices and were removed from the
distribution. Their SHA-256 identifiers were:

```text
baseline        db3eba472d4fd63e94268d0a649c958181590ed734e36dbc1e48801ba6dc8a95
endpoint        7366ede52ce97e24930d3be62aa49b02f52f97e257aa2d383fb29dc084916ad2
gas-only        a15449178e15faf1e4ab2c5755384bc979860bc71c6971f8470304072d19a6eb
downward-order  2c3627a69d46bb4bcb44c5ec85df0d3291392ac09a1e289390855a281816695b
local-beta      6bf513bf28705302d44a65c32a1c4efe952788eecfc9808af597cae6a48b157a
stable-amplitude bc178a40488e819369922ae6721717a92ba2cb2e34a3700142c187b093564093
basis-opt       f8a332e3806bf27c5af9acbd0f11e47c05f71723ade51b57301580ec711f486e
rayleigh-fast   53d68444e1f62ccf36309be3669d2083ba0b698be5348de87510113c00bee933
```

## Matrices

```text
PSSA_NUMERIC_FIXES_RAYLEIGH_MATRIX_2026-07-19.csv       72 cases
PSSA_NUMERIC_FIXES_AEROSOL_MATRIX_2026-07-19.csv        12 cases
PSSA_NUMERIC_FIXES_PP_FLAT_MATRIX_2026-07-19.csv        27 cases
PSSA_NUMERIC_FIXES_GAS_ONLY_OCEAN_MATRIX_2026-07-19.csv 24 cases
PSSA_CORRECTION_BEFORE_AFTER_2026-07-19.csv               8 cases
```

The JSON summaries keep the maxima and the metadata of representative cases. For Q and U use
absolute differences or the normalized `|delta Q,U| / |I|` values in
`PSSA_NUMERIC_FIXES_ROBUST_SUMMARY_2026-07-19.json`; raw relative percentages are unstable where
Q or U crosses zero.

## Reproducible checks on the current code

```bash
./scripts/smoke_pssa_numeric_fixes.sh ./build/ocrt
./scripts/test_pssa_internal.sh
```

The historical differential matrices need the previous / intermediate binaries listed above and
are therefore audit records, not a standalone rerun harness. (Note: the legacy PSSA scheme these
records refer to was removed in v1.11; the records are kept as history.)

---

# PSSA 수치 수정 — 검증 자료

날짜: 2026-07-19

## 범위

이 디렉터리의 행렬은 직전 `pssa-aodref-radiometry-split-layerwarn` 패키지 대비 네 가지 변경을 정량화한다.

1. 반사빔 하단 끝점 감쇠의 부호,
2. 기체 전용 해양 PSSA 직달 투과 보조 함수,
3. 반사 하향 1차 산란을 고차 SOS 차수에 결합,
4. 층 국소 반사빔 `beta` 위상 원천.

## 소스 순서

```text
b92a476 baseline package
f07a0e7 endpoint sign
4a831e5 gas-only ocean helper
8a6a380 downward first-order coupling
ff46823 layer-local beta source
f284cb4 direct spherical reflected amplitude
27134f1 pure-Rayleigh beta setup optimization
```

중간 바이너리는 차분 행렬을 만드는 데 쓴 뒤 배포에서 제거했다. 그 SHA-256 은 다음과 같다.

```text
baseline        db3eba472d4fd63e94268d0a649c958181590ed734e36dbc1e48801ba6dc8a95
endpoint        7366ede52ce97e24930d3be62aa49b02f52f97e257aa2d383fb29dc084916ad2
gas-only        a15449178e15faf1e4ab2c5755384bc979860bc71c6971f8470304072d19a6eb
downward-order  2c3627a69d46bb4bcb44c5ec85df0d3291392ac09a1e289390855a281816695b
local-beta      6bf513bf28705302d44a65c32a1c4efe952788eecfc9808af597cae6a48b157a
stable-amplitude bc178a40488e819369922ae6721717a92ba2cb2e34a3700142c187b093564093
basis-opt       f8a332e3806bf27c5af9acbd0f11e47c05f71723ade51b57301580ec711f486e
rayleigh-fast   53d68444e1f62ccf36309be3669d2083ba0b698be5348de87510113c00bee933
```

## 행렬

```text
PSSA_NUMERIC_FIXES_RAYLEIGH_MATRIX_2026-07-19.csv       72 케이스
PSSA_NUMERIC_FIXES_AEROSOL_MATRIX_2026-07-19.csv        12 케이스
PSSA_NUMERIC_FIXES_PP_FLAT_MATRIX_2026-07-19.csv        27 케이스
PSSA_NUMERIC_FIXES_GAS_ONLY_OCEAN_MATRIX_2026-07-19.csv 24 케이스
PSSA_CORRECTION_BEFORE_AFTER_2026-07-19.csv               8 케이스
```

JSON 요약은 최댓값과 대표 케이스 메타데이터를 담는다. Q/U 는 절대차 또는
`PSSA_NUMERIC_FIXES_ROBUST_SUMMARY_2026-07-19.json` 의 정규화 `|delta Q,U| / |I|` 값을 쓴다. Q 나 U 가 0 을
지나는 곳에서는 원시 상대 백분율이 불안정하다.

## 현재 코드에서 재현 가능한 점검

```bash
./scripts/smoke_pssa_numeric_fixes.sh ./build/ocrt
./scripts/test_pssa_internal.sh
```

과거 차분 행렬은 위의 이전/중간 바이너리가 필요하므로 독립 재실행 하니스가 아니라 감사 기록이다.
(참고: 이 기록이 가리키는 옛 PSSA 방식은 v1.11 에서 제거됐다. 기록은 이력으로 보존한다.)
