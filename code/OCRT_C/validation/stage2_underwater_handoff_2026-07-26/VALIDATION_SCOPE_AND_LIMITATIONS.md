# 제공 수중 정합자료의 검증범위와 제한 — 2026-07-26

## 1. 내부 무결성 검증

제공된 방향별 원자료 `OCRT_OSOAA_underwater_validation_2026-07-26.csv`
(11,154행 × 37열)에서 조합별 지표를 독립 재계산하였다. 재계산 결과는 제공
요약표 `OCRT_OSOAA_underwater_summary_2026-07-26.csv`의 50행 × 26열과 모든
수치열에서 최대 절대차 0으로 일치하였다. 18개 Rrs/rrs I/Q/U 산포도도 원자료에서
독립 재생성하였다.

## 2. 생산 구성모델과 제공 Chl 검증경로의 차이

현재 생산 구성모델은 EAP 식물플랑크톤 산란을 비활성화한다. Chl은 식물플랑크톤
흡수에는 들어가지만 `b_phyto=bb_phyto=0`이며, Chl 연계 입자산란은 detritus 하나가
담당한다. `--ocrt-phyto-group`을 Chl>0과 함께 지정하면 종료코드 2로 실패한다.

제공된 Chl 정합자료는 생산 구성모델을 그대로 최종 RT까지 실행한 결과가 아니다.
`ff_chl2.py`는 구성모델을 한 번 실행하여 scalar IOP를 추출한 뒤, 다음 고정벌크
우회경로로 OCRT를 실행한다.

1. detritus 위상만 먼저 절단한다.
2. 절단된 detritus와 pure-water Rayleigh를 산란계수로 혼합한다.
3. `b = b_w + b_det(1-A/2)`로 산란계수를 재조정한다.
4. `--water-model iop --iop-mie-phase` 경로로 RT를 실행한다.

따라서 제공 Chl rrs I 오차 0.1321–0.2483%는 **고정벌크 우회경로의 검증값**이다.
현재 생산 구성모델은 구성성분별 절단을 fail-loud로 막고 있으므로 이 수치를
생산 구성모델의 직접 종단 정합성으로 해석하면 안 된다.

## 3. 독립 재실행 범위

이번 통합에서는 제공 CSV의 통계·그림을 독립 재계산하였으며, 다음 코드 회귀를
실행하였다.

- 17종 EAP phase API 정규화·Mueller physicality·결정성·오류코드
- 명시적 species+Chl fail-loud
- Chl absorption 유지, EAP b/bb=0, detritus-only phase
- 구성모델 truncation flag fail-loud 및 fixed-bulk truncation 경로 보존
- 비-Chl CDOM/TSM 18조건 stdout/stderr byte identity
- water RAA, FIX123, FIX4 source/operator closure, exact-pole, coupling clamp 회귀

OSOAA 전체 원시 실행은 이번 패키지에 수정된 OSOAA 실행파일·완전한 cumulative
patch·surface-matrix 해시가 모두 포함되어 있지 않으므로 새로 수행하지 않았다.
따라서 제공 CSV는 **내부 일관성이 확인된 전달 기준자료**이며, 완전한 외부 독립
재현 기준선으로 승격하려면 OSOAA provenance를 추가 확보해야 한다.

## 4. EAP 원자료 provenance

패키지에는 생성된 `src/generated/rt_eap_species_data.inc`가 포함되지만, 이 파일의
원자료인 다음 두 파일은 포함되지 않았다.

- `EAP_invivo_means.csv` — SHA-256
  `92051b77cec90744294536dbfe611b65350532b98d7d0f8c32545bb5e4f273a8`
- `501nm_extended_e1701000.mat` — SHA-256
  `50638c22d4596c38da5f34e10d846188faa225d6fb262409f858fce6e9725c02`

따라서 동결된 생성 catalog에서 17종 phase를 재생성할 수는 있으나, 원자료부터
catalog를 다시 만드는 단계는 이번 전달물만으로는 재현할 수 없다.

## 5. 수면 위 Rrs 해석

OSOAA Rrs(0+)는 full-water와 black-water의 두 큰 level-26 값을 차감하여 산출한다.
`black_frac_*`가 큰 조건에서는 상대오차가 수치적으로 증폭된다. Rrs 정확도 수치를
인용할 때에는 반드시 동일 행의 `black_frac_I/Q/U`를 함께 제시해야 한다.
