# OCRT RT 기반 상·하향 투과율 출력 수정

## 판정

기존 공개 출력의 다음 식은 방향별 RT 결과가 아니므로 제거했다.

```text
T_diff_up_hemi  = T_diff_dn_hemi
T_total_up_view = exp(-tau/mu_view) + T_diff_up_hemi
```

수정본은 실제 해수면 상향 BRDF가 atmospheric bottom-source SOS를 통과한
결과로 상향 유효 투과율을 계산한다.

```text
TOA_water_signal_I = I_TOA_total - I_atm_path
T_total_up_view = TOA_water_signal_I / Lu0plus
T_diff_up_view = T_total_up_view - T_dir_up_view
```

하향 확산투과율은 원래부터 대기 SOS의 BOA 확산 radiance를 반구 적분한
RT 값이었다. 이번 수정에서 총 하향값 `T_total_dn_hemi`를 명시적으로 추가했다.

엄밀한 bottom-source pass가 없거나 상향 source가 정의되지 않은 경우에는
근사값을 출력하지 않고 `NaN`, `T_up_rt_valid=0`을 반환한다.

## 출력-only 보장

대표 no-atmosphere, Rayleigh, Rayleigh+aerosol 해양 조건에서 수정 전후의
공통 비투과율 숫자 출력은 모두 동일했다. TOA I/Q/U, Rrs, rrs, Lu, Ed/Eu,
IOP, order 및 convergence를 변경하지 않는다.

## AC LUT 수정

`scripts/ocrt_ac_lut.py`의 `Tup(vza)=Tdn(vza)` 상반성 재구성을 제거했다.
Black/black-Fresnel 계산은 상향 BOA source가 없으므로 `Tup_total=NaN`,
`Tup_rt_valid=0`으로 기록한다. 정확한 상향값은 ocean-coupled 실행의
`T_total_up_view`를 사용해야 한다.

## 별도 발견 이슈

Ocean full-grid driver는 현재 aerosol object를 solver에 전달하지 않는다.
이는 투과율 해석식 문제가 아니라 별도의 aerosol wiring 결함이며, 수정 시
full-grid TOA/Rrs가 변하므로 이번 output-only patch에서 제외했다.

## 주요 파일

- `patches/OCRT_v1.2_RT_DERIVED_TRANSMITTANCE_OUTPUT_CORE_20260721.patch`
- `docs/OUTPUT_RT_DERIVATION_AUDIT_2026-07-21.md`
- `validation/before_after_representative.csv`
- `validation/physical_output_invariance_scatter.png`
- `validation/upward_transmittance_old_approx_vs_rt.png`
- `APPLY_TO_CURRENT_HOTFIX.md`
