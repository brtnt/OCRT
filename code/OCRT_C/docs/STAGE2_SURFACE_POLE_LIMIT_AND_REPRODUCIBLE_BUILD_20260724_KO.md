# 2단계 수면 연산자 극한 및 재현 빌드 규약

## 1. 적용 범위

본 규약은 다음 네 rough-interface Mueller kernel의 자오면 회전에 적용한다.

- 공기측 반사 `surface_R_coxmunk_trig`
- 물→공기 투과 `surface_T_coxmunk_trig`
- 수중 내부반사 `surface_R_ww_coxmunk_trig`
- 공기→물 투과 `surface_T_aw_coxmunk_trig`

기존 FIX1–FIX4, 공개 water RAA/U 재구성 및 대기 하네스는 변경하지 않는다.

## 2. 정확한 극 한계

자기 방향의 `sin(theta)`가 0으로 접근할 때 일반 회전식은 0/0이 된다. 분모에 작은 바닥값만 적용하면 극에서 잘못된 기준면 회전이 생성된다. 해석 극한은 다음과 같다.

```text
cos(sigma) -> -sign(mu_self) cos(phi)
sin(sigma) -> sin(phi)
```

`OCRT_SURF_POLE_EPS = 1e-8` 이하에서 위 식을 적용한다. 비극 구간의 기존 식은 변경하지 않는다.

## 3. 양쪽 극 축퇴

입사·출사 방향이 모두 연직이면 상대방위각 자체가 정의되지 않는다. 생산 격자에서 이 배치는 삽입 특수절점끼리의 쌍이며 적어도 한쪽 구적 가중치가 0이므로 경계 적분에 기여하지 않는다. 따라서 임의의 추가 물리 보정을 적용하지 않는다.

## 4. 빌드 재현성

릴리스 빌드의 기본 ISA를 `-march=cascadelake`로 고정한다. `-march=native`는 실행 호스트가 바뀌면 벡터화·FMA·연산 결합 순서를 바꾸므로 바이너리 해시와 최하위 비트를 재현할 수 없다.

```bash
./scripts/build_release_v1.2.sh build/ocrt
```

다른 CPU 기준이 필요한 경우에만 명시적으로 지정한다.

```bash
OCRT_MARCH=x86-64-v3 ./scripts/build_release_v1.2.sh build/ocrt
```

동일 해시 주장은 동일 소스, 컴파일러·링커 버전, 플래그, 환경변수, `OCRT_MARCH`에 한정한다.

## 5. 비트 동일 CSV 검사

`pandas`의 고속 CSV 파서는 매우 작은 수의 마지막 비트를 변경할 수 있다. 비트 동일 판정에는 다음 중 하나를 사용한다.

```bash
cmp file_a.csv file_b.csv
python scripts/compare_csv_bitexact.py file_a.csv file_b.csv
```

후자는 Python 내장 `float()`로 다시 파싱한 binary64 비트열을 비교한다.
