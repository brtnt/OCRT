# OCRT 누적 작업지시서 — 2026-07-25 ~ 07-26 세션

대상: 메인 코드 개발 세션
기준선: `MIGRATION_OCRT_2026-07-24.tar.gz` 를 푼 v3 상태에 v4·v5 패치를 적용한 트리
최종 빌드: `ocrt_v11`, md5 `1c0d88e7e31175b9cb94fd13b8216c74`
빌드 명령: `gcc -std=c11 -O3 -march=cascadelake -ffp-contract=fast -fassociative-math -fno-signed-zeros -fno-trapping-math -DOCRT_FAST_KERNELS -fopenmp -Isrc $(find src -name '*.c') -o build/<bin> -lm`

---

## 0. 요약

이 세션은 OCRT 와 OSOAA 의 수중 복사전달 정합성을 세 계열(CDOM, 총부유물, 엽록소)에서 확인하고 닫았다. 그 과정에서 코드 결함 다섯 건과 검증 규약 결함 두 건을 찾아 고쳤다. 가장 중요한 조치는 **종별 식물플랑크톤 산란을 비활성화한 것**이다(3절).

---

## 1. 검증 정책 변경 — 단일 기하에서 전각도로

**변경 전.** 기존 검증기(`harness_prev/scripts/stage_run.py`)는 관측천정각을 60도 한 점으로 고정하고 방위각만 4종을 돌았다.

**문제.** v4·v5 수정의 효과는 전부 천저(관측천정각 0도)에서 나타난다. 60도만 보면 수정의 효과가 전혀 보이지 않는다. 실제로 천저의 `Rrs(0⁺)` Q 오차는 v3 에서 8.366 %, v5 에서 0.062 % 로 135 배 차이가 나는데 이것이 보이지 않았다.

**변경 후.** 두 코드 모두 한 번 실행에 전 각도를 산출하므로 전각도 비교로 바꿨다.

- OCRT: `--output-full-grid` 에 `--lut-vza-step`, `--lut-vza-max`, `--lut-raa-step` 을 주면 한 번의 solve 로 전 격자를 낸다.
- OSOAA: 상세출력을 관측천정각에 대해 PCHIP 로 보간하면 한 실행에 전 관측천정각이 나온다. 방위각만 실행당 하나다.
- 실행 부담이 오히려 줄었다. 조건당 OCRT 1 회 + OSOAA (방위각 수 × 2) 회다.

**확인해 둔 사항.** 격자 모드와 단일 기하 모드의 값이 최대 상대차 4.37e−07 로 일치한다. 비트 동일은 아니다. 각도 링 구성이 달라 푸리에 재구성의 연산 순서가 달라지기 때문이다. 따라서 **두 모드의 값을 비트 단위로 섞어 쓰면 안 된다.**

---

## 2. 코드 수정 — 신규 EAP 위상 계산기 통합

첨부 `OCRT_EAP_PHASE_20260726` 묶음을 트리에 넣었다.

**배치한 파일**

| 파일 | 성격 |
|---|---|
| `src/rt_eap_mie_phase.{c,h}` | 신규 |
| `src/internal/rt_eap_coated_mie.{c,h}` | 신규 |
| `src/internal/rt_eap_species_catalog.{c,h}` | 신규 |
| `src/generated/rt_eap_species_data.inc` | 신규(생성물) |
| `tools/eap_mie_writer.c`, `tools/eap_generate_all.sh` 등 | 신규 |
| `inputs/eap_source/*` | 신규(원자료) |
| `inputs/water_iop/EAP_00..16_*.mie` | 생성물 17 종 |

**생성 방법과 소요 시간**

```
gcc -std=c11 -O3 -march=cascadelake -Isrc -o build/eap_mie_writer \
    tools/eap_mie_writer.c src/internal/rt_eap_species_catalog.c \
    src/internal/rt_eap_coated_mie.c -lm
sh tools/eap_generate_all.sh representative inputs/water_iop
```

유효직경 0.5~6 µm 는 종당 1~20 초이나, 24 µm 두 종(와편모조류·라피도조류)은 종당 276 초가 걸린다. 전체 생성에 약 15 분을 잡아야 한다. 동봉 표본 4 종과 생성물이 바이트 단위로 일치함을 확인했다.

**통합 문서의 오류 한 건 — 반드시 확인할 것**

통합 문서 3.2절 (c) 는 `rt_iop_organic_group_parse()` 가 이름표를 순회하므로 새 이름을 자동으로 받는다고 적고 있다. **이 기준선에서는 사실이 아니다.** 해당 함수는 `strcmp`/`strstr` 로 세 그룹을 직접 판정하도록 박혀 있었다. 그대로 두면 다음이 오류 없이 잘못 걸린다.

- `eap_synechococcus` → `strstr(norm,"synechococcus")` 에 걸려 **pico** 로 간다
- `eap_diatoms_pennate` → `strstr(norm,"diatom")` 에 걸려 **micro** 로 간다

이름표 정확 일치를 먼저 확인하고 기존 별칭을 뒤로 물리도록 다시 썼다(`src/rt_iop_organic.c`).

---

## 3. 코드 수정 — 종별 식물플랑크톤 산란 비활성화 **(가장 중요)**

### 3.1 왜 껐는가

수중 솔버는 위상함수를 **L=200 차 르장드르 전개**로 받는다(`rt_water_rt.c` 의 `FIXED_BULK_BETAL_LMAX`). 종별 EAP 위상은 전방 첨두가 가팔라 이 차수에 담기지 않는다. 그 결과 전개를 되살린 위상이 **중간 산란각에서 음수**가 된다. 물리적으로 불가능한 값이다.

측정 결과(412 nm, L=200 전개 되살림, 30~150도 241 점 중 음수 개수 / |β₂₀₀/β₂|):

| 그룹 | 유효직경 | 음수 개수 | β₂₀₀/β₂ |
|---|---|---|---|
| `pico`, `eap_prochlorococcus`, `eap_synechococcus` | ≤ 1.2 µm | **0** | ~5e−05 |
| `eap_pelagophytes`, `eap_prasinophytes` | 3 µm | 42~45 | 6e−02 |
| `eap_hapto_prymnesiaceae` | 4 µm | 95 | 3.0e−01 |
| `eap_microcystis` | 5 µm | 71 | 1.0e+00 |
| 6 µm 계열 7 종, `nano`(5 µm), `micro`(20 µm) | 5~20 µm | 112~120 | 0.6~8.7 |
| `eap_dinoflagellates`, `eap_raphidophytes` | 24 µm | 121 | 3.3e+01 |

경계는 1.2 µm 와 3 µm 사이다. **이는 자료 결함이 아니라 표현 한계다.** 신규 EAP 자료의 진동 문제는 해결됐다(이웃 각도 사이 변동 최대가 131.8 % → 3.9 %). 남은 것은 전방 첨두이며 이는 물리적 실체다.

### 3.2 무엇을 껐는가

**첫째, 종을 지정하면 멈춘다.** `--ocrt-phyto-group` 을 엽록소가 양수인 구성모델 실행에 주면 종에 관계없이 메시지를 내고 종료한다(코드 2). `main.c` 의 `organic_group_cli_seen` 을 조건으로 쓴다.

```
error: --ocrt-phyto-group is not supported yet (requested '...').
  Reason: the in-water solver expands the phase function to L=200. ...
  Drop the option and rerun: phytoplankton then contributes absorption only,
  and the particle scattering phase is detritus alone.
  Per-species scattering will be restored in a later revision.
```

**둘째, 종을 지정하지 않으면 식물플랑크톤은 흡수만 쓴다.** 산란과 후방산란을 총량에서 빼고 성분값도 0 으로 만든다. 입자 위상은 쇄설물 하나가 된다.

| 파일 | 변경 |
|---|---|
| `src/rt_types.h` | `int organic_phyto_scattering;` 추가, 기본 0 |
| `src/rt_water_rt.h` | 같은 항목 추가 |
| `src/rt_solver.c` | `w_opts.organic_phyto_scattering = cs->organic_phyto_scattering;` |
| `src/rt_water_rt.c` | 총량과 성분에서 `b_phyto`·`bb_phyto` 제거 |
| `src/rt_iop_organic.{h,c}` | `rt_iop_organic_group_is_large()` 와 판정표(측정값 주석 포함) |
| `src/main.c` | 종 지정 시 종료, 도움말 갱신 |

되살릴 때는 `organic_phyto_scattering` 을 1 로 두면 원래 동작으로 돌아간다. 판정표는 어느 종이 안전한지 판단하는 근거로 코드에 남겨 두었다.

### 3.3 부작용 — 기본값이 막힌다

`--ocrt-phyto-group` 의 기본값은 `micro` 다. 따라서 **엽록소를 넣고 종을 지정하지 않으면 문제가 없고, 지정하면 무조건 멈춘다.** 기존 스크립트가 종을 명시했다면 전부 수정해야 한다.

### 3.4 되살리려면 무엇이 필요한가

두 갈래가 있다.

1. **델타 절단을 구성모델 경로에 구현한다.** 절단을 걸면 전방 첨두가 사라져 L=200 으로 충분해진다. 실측으로 확인했다. 절단 적용 시 와편모조류(24 µm)까지 포함해 모든 종이 음수 0 개, 되살림 오차 0.15 % 가 된다.
2. **L 을 올린다.** L=400 이면 절단 없이도 진동이 사라지고 음수가 0 이 된다. 다만 `FIXED_BULK_BETAL_LMAX` 가 소스 상수라 옵션으로 빼는 작업이 따르고, OSOAA 는 `OS_NB` 가 200 고정이라 정합 시험에는 쓸 수 없다.

논문 자료 생산 정확도만 놓고 보면 2 번이 근사가 아니라는 점에서 낫다. 정합 시험까지 고려하면 1 번이 맞다.

---

## 4. 코드 수정 — 절단 플래그가 조용히 무시되던 문제

**증상.** 구성모델 경로에서 `--ocrt-mie-truncation` 을 켠 값과 끈 값이 자릿수까지 같았다(`rrs0minus_I=8.541482e-04`, orders=18 로 동일).

**원인.** 절단을 수행하는 코드는 고정벌크 사전검사 블록에만 있고, 그 블록은 `--iop-mie-phase` 가 있을 때만 돈다. 성분 모멘트를 계산하는 `water_component_phase_moments()` 에는 절단 인자 자체가 없다.

**조치.** 명시적 미구현 오류로 막았다.

```
rt_water_rt: --ocrt-mie-truncation/--ocrt-mie-ss-mode are not implemented for the
OCRT constituent water model.  Use the fixed-bulk path
(--water-model iop --iop-mie-phase F --iop-mie-truncation).
```

**구현하려면 네 가지가 함께 필요하다.**

1. 성분 모멘트 계산에 절단 경로를 붙인다(`rt_aerosol_compute_vector_legendre_gauss_truncated`, 상수 0.85 / 0.92 / 0.1).
2. 산란계수를 `b* = b(1 − A/2)` 로 다시 잡는다.
3. 후방산란 비를 절단된 위상에서 다시 뽑는다.
4. 위상 캐시 키에 절단 여부를 넣는다.

성분이 둘 이상이면 성분별 제거 질량이 서로 달라 단일 규칙이 아예 성립하지 않는다는 점도 남는다. 이번에 식물플랑크톤 산란을 끄면서 성분이 하나가 되었으므로 그 제약은 당분간 걸리지 않는다.

**부수 수정.** 성분 수 하한을 2 에서 1 로 낮췄다(`rt_water_rt.c`). 예전 코드는 "엽록소가 양수면 항상 식물플랑크톤과 쇄설물이 함께 들어온다"를 전제했는데 3절 변경으로 그 전제가 깨졌다.

---

## 5. 검증 규약 변경 두 건

### 5.1 해수 광학두께 — 수송 광학두께 기준으로 바꿨다

**기존 규약.** OSOAA 해수 깊이를 광학두께 6 이 되도록 잡는다.

**문제.** 그 값은 순수 해수·CDOM 처럼 단일산란 반사도가 낮은 계열에서 정해진 것이다. 입자가 들어가면 단일산란 반사도가 0.9 를 넘고 전방산란이 강해 광학두께 6 은 반무한이 아니다. 총부유물 2 g/m³, 443 nm 에서 OSOAA 가 OCRT 보다 **10.25 % 어두웠다.**

**새 규약.**

```
tau * (1 - omega * g) = 6      →      tau = 6 / (1 - omega*g)
```

- 순수 해수·CDOM 은 단일산란 반사도가 작아 자동으로 종전의 6 을 돌려준다. 기존 결과와의 연속성이 유지된다.
- 총부유물 2 g/m³, 443 nm 에서는 tau=24 가 나오고 실측 최적점과 일치한다.
- 적용 후 수면 아래 세기 오차가 8.008 % → 0.191 % 로 42 배 줄었다.

`compare_entry.py` 의 `check_water_optical_depth()` 상한(12)이 이 규약과 충돌하므로 입자 케이스에서는 우회해야 한다.

### 5.2 물 분자 산란은 절대 자르지 않는다

고정벌크 우회로 절단을 걸 때, 물 레일리를 위상에 섞은 뒤 혼합 위상 전체에 절단을 걸면 **레일리 산란까지 76 % 가 제거되어** 결과가 어긋난다. 엽록소 0.3 mg/m³, 443 nm 에서 2.764 % 오차가 났다.

쇄설물만 먼저 자르고 잘린 위상에 물 레일리를 섞은 뒤, 산란계수를 `b = b_w + b_입자(1 − A/2)` 로 직접 주면 **0.132 % 로 21 배 줄어든다.**

관련해서 반대 방향 함정도 있다. 물 레일리를 아예 빼면 후방산란이 절반이 되어 **55 % 어긋난다.** 엽록소 1 mg/m³, 443 nm 에서 순수 해수는 산란의 1.26 % 뿐인데 **후방산란의 52 %** 를 낸다(레일리의 bb/b 가 0.5 이고 입자는 0.006 이기 때문이다).

---

## 6. 검증 결과 — 세 계열 모두 닫힘

수면 아래 반사도 `rrs(0⁻)` 의 세기 성분 기준이다. 정규화는 조합별 OSOAA 최대 크기다.

| 계열 | 조합 수 | 방향 수 | 세기 오차 범위 |
|---|---|---|---|
| CDOM (3단계 × 5밴드) | 15 | 2,535 | 0.04 ~ 0.34 % |
| 총부유물 (2단계 × 3밴드) | 6 | 1,014 | 0.05 ~ 0.49 % |
| 엽록소 (3단계 × 3밴드) | 9 | 1,521 | 0.13 ~ 0.25 % |

**부가 확인**

- 풍속 3·5·7·10 m/s: 세기는 무관, 편광은 단조 증가하나 10 m/s 에서도 수면 아래 최대 0.76 % 이내.
- 태양천정각 0·20·40·60 도: 세기는 각도가 클수록 좋아지고 편광은 반대이나 최대 1.05 % 이내.
- 태양천정각 0 도 축대칭 검사: OCRT 의 U 가 배정밀도 한계(1.06e−16), 방위각 변동이 2.13e−16. **v3 도 통과한다.** 즉 v4 극 한계 수정의 실효는 관측이 천저인 극에 있고 태양이 천정인 극은 v3 에서도 옳았다.
- 수중 절점 64 대 96: CDOM 계열에서 차이가 0.001 % 이내로 수렴해 있다. 입자 케이스는 별도 확인이 필요하다.

**수면 위 `Rrs(0⁺)` 를 읽을 때의 주의.** OSOAA 의 수면 위 값은 흑색 해수 실행을 한 번 더 돌려 빼야 얻어진다. 수출광이 작은 조합에서는 두 큰 값의 차가 작아져 상대오차가 부풀려진다. 실제로 CDOM 1.0 · 412 nm 에서 흑색 성분 비가 0.9 를 넘으며 오차가 13.8 % 로 뛴다. **수면 위 성적은 흑색 차감 비중(`black_frac_*`)을 함께 적지 않으면 쓰면 안 된다.**

---

## 7. 미결 항목

| 항목 | 내용 | 우선도 |
|---|---|---|
| 구성모델 절단 구현 | 4절 참조. 현재 논문 자료 생산 경로에 4~6 % 편향이 그대로 있다 | 높음 |
| 종별 식물플랑크톤 산란 복원 | 3.4절 참조. 위 항목이 되면 함께 풀린다 | 높음 |
| `FIXED_BULK_BETAL_LMAX` 옵션화 | 현재 소스 상수 200 | 중간 |
| 에어로솔 포함 검증 | `.mie` 자산 필요 | 중간 |
| `harness/phase_gate.py` | 유실됨. 이번 세션에서 대체 도구를 새로 짰다 | 낮음 |
| 죽은 덤프 훅 정리 | `OCRT_DUMP_WATER_MIE_MOMENTS`, `OCRT_DUMP_PHASE`, `OCRT_DUMP_FIXEDPHASE` 는 등록만 되어 있고 출력 코드가 없다 | 낮음 |

---

## 8. 회귀 확인 기록

모든 변경은 다음을 통과했다.

- 엽록소를 쓰지 않는 경로(CDOM 단독, 총부유물 단독, 풍속·밴드 조합) 18 건 **바이트 동일**
- v5 검증 프로토콜 `solver_surf`: V4 솔버 기준값 불일치 0 항목, V5 영향 범위 위반 0 조건, V6 천저 극 제약 성립
- v3 기준 빌드 지문 `1135ba699cac`, v5 `22947f48e7ab` 재현 확인(`-march=cascadelake`)

---

## 9. 이 세션에서 배운 검증 원칙

1. **실행 결과값과 종료 코드를 함께 확인한다.** 진단 덤프는 실패 지점보다 앞에서 출력되므로 덤프만 보고 성공을 판정하면 안 된다. 이 세션에서 실제로 한 번 잘못 판정했다.
2. **플래그를 여는 변경은 그 플래그가 실제로 작동하는지 값으로 확인한다.** 관문만 열고 동작이 없으면 조용히 무시된다.
3. **0 에 가까운 양을 자기 자신으로 정규화하지 않는다.** 태양천정각 0 도에서 Q 를 Q 평균으로 정규화해 26 배 변동이라는 허수를 만든 적이 있다. 세기 최대 크기로 정규화하면 2.1e−16 이다.
4. **조합을 섞어 한 번에 정규화하지 않는다.** 조합마다 신호 크기가 200 배까지 차이 나면 큰 조합이 분모를 독점한다.
