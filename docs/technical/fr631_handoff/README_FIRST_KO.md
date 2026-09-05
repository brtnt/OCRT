# OCRT 330–1100 nm FR631 분광자료 인계 세트 (2026-08-20)

이 세트는 2026-08-14 단일 ZIP의 불완전/축소 문제를 대체한다. 198개 canonical Mie 파일을 전수 재생성·검증했으며, 각 파일은 다음 계약을 만족한다.

- bulk spectral table: 330–1100 nm, 1 nm, 771행
- P11/P12/P33 phase table: 각 330–1100 nm, 1 nm, 771열
- scattering angle: FR631 exact grid, 631점, 180°→0°
- 직접 Mie 계산 각도: FR631 전 점; 성긴 각도 계산 후 각도보간 금지
- 소비 보간: 파장·각도 모두 선형 런타임 보간 패치 포함

필수 다운로드는 PART00–PART10이다. PART11은 생성 재현용 source package로 런타임 설치에는 필수가 아니다.

## 추출

모든 필수 ZIP을 같은 상위 폴더에 풀어 동일한 `OCRT_SPECTRAL_DATA_330_1100_FR631_HANDOFF_2026-08-20` 디렉터리에 병합한다. 그 뒤:

```bash
python OCRT_SPECTRAL_DATA_330_1100_FR631_HANDOFF_2026-08-20/05_TOOLS/verify_extracted_bundle.py
```

PASS 확인 후 dry-run:

```bash
python OCRT_SPECTRAL_DATA_330_1100_FR631_HANDOFF_2026-08-20/05_TOOLS/apply_data_patch.py --c-root <OCRT_C_ROOT> --python-root <OCRT_PY_ROOT>
```

실제 설치:

```bash
python OCRT_SPECTRAL_DATA_330_1100_FR631_HANDOFF_2026-08-20/05_TOOLS/apply_data_patch.py --c-root <OCRT_C_ROOT> --python-root <OCRT_PY_ROOT> --apply
```

먼저 `03_INTEGRATION_PATCH`의 C/Python loader 패치를 적용해야 한다. 기존 8192-byte C line buffer로는 771-column phase row를 읽을 수 없다.

## 분할 ZIP 목록

| 파일 | 용도 | 파일 수 | ZIP bytes | SHA-256 |
|---|---|---:|---:|---|
| PART01_WATER_MIE_22.zip | mie | 22 | 166,427,065 | `35e869cf49e910c0c5994159d94605e0b214372b15c312b6c04f76d5b7615bf8` |
| PART02_AEROSOL_SNF_OPAC_16.zip | mie | 16 | 119,213,743 | `2553bac456911974a09e87a3ffcd6da4bb2240762118d1acca49467032af2ba5` |
| PART03_AEROSOL_AHMAD2010_PAPER_RH30_RH50.zip | mie | 20 | 153,597,557 | `c6f59b739119b3184d0b8b71d980b53ae6b35a112b1eaccdccd3eeb624a1d2aa` |
| PART04_AEROSOL_AHMAD2010_PAPER_RH70_RH75.zip | mie | 20 | 153,496,170 | `8f61c27542baa4e88ff61b1bf8dd81e01468203bd3499c0a127bdf70691b08e9` |
| PART05_AEROSOL_AHMAD2010_PAPER_RH80_RH85.zip | mie | 20 | 153,282,507 | `f931d958d02a7f66d1ba095ba4d901ffe419e13cf98c845602b80f9914b11661` |
| PART06_AEROSOL_AHMAD2010_PAPER_RH90_RH95.zip | mie | 20 | 153,155,833 | `c0ef0c5868c57559848d161726bd9f7af68e93ac66202780949332d781df9eb5` |
| PART07_AEROSOL_AHMAD2010_ACCURT_RH30_RH50.zip | mie | 20 | 153,309,899 | `1afb244e30fe64db68b3e06acbdee2b12d480702dc60f32aadfb0c80d037cc21` |
| PART08_AEROSOL_AHMAD2010_ACCURT_RH70_RH75.zip | mie | 20 | 153,552,555 | `1a340588b07ff621fede0283282d96ad0a957c495e3a1236a53e55370d12cc45` |
| PART09_AEROSOL_AHMAD2010_ACCURT_RH80_RH85.zip | mie | 20 | 153,281,972 | `dcbf63bcaa086865f204f033990356f958d626e4de875ab6ebfe1d633cf688bf` |
| PART10_AEROSOL_AHMAD2010_ACCURT_RH90_RH95.zip | mie | 20 | 153,152,650 | `0f0f1c25dd91c882f920c05ff4dcb843e10cd4fc9cb884360b00cdaf4fa204d3` |
| PART11_REPRODUCTION_SOURCES.zip | repro | 324 | 78,520,006 | `70ed8a8c4c324439aba6b118b145f6c4315506f863403d2c03cb4fdef5efa249` |

PART00은 문서·non-Mie·코드패치·검증도구를 담고, PART01–PART10이 198개 Mie를 담는다. PART11은 선택적 재현자료다.
