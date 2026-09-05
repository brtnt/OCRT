# OCRT AOD 기준파장 입력 규약

## 원칙

사용자는 현재 계산 파장의 AOD를 입력하지 않는다. 하나의 고정 기준파장에서
정의된 AOD를 입력하고, OCRT가 선택한 aerosol `.mie` 파일의 분광 소산계수 비를
이용해 각 계산 파장의 AOD를 내부 산출한다.

```text
AOD(lambda) = AOD(lambda_ref) * Ext(lambda) / Ext(lambda_ref)
```

`.mie` spectral table의 `Nor_Ext_Co`를 사용한다. `Extinct_Co`의 비를 사용해도
동일하지만, 정규화 소산계수가 이 목적에 직접 맞는다.

## 권장 옵션

```bash
--aod-555 0.10
--aod-865 0.05
```

일반 기준파장은 다음처럼 지정한다.

```bash
--aod 0.10 --aod-ref-wavelength 550
```

`--aod` 단독은 호환 옵션이며 기본 기준파장은 555 nm이다.
`--aod-band`는 제거되었고 오류를 반환한다.

## 출력 메타데이터

```text
AOD_ref
AOD_ref_nm
AOD_band
AOD_ext_ratio
```

- `AOD_ref`: 사용자가 입력한 기준 AOD
- `AOD_ref_nm`: 기준파장
- `AOD_band`: 현재 계산파장으로 변환된 AOD
- `AOD_ext_ratio`: `Ext(lambda)/Ext(lambda_ref)`

## 예

M80C에서 `--aod-555 0.1`, 계산파장 443 nm:

```text
AOD_ref=0.1
AOD_ref_nm=555
AOD_ext_ratio=1.0428381863
AOD_band=0.1042838186
```
