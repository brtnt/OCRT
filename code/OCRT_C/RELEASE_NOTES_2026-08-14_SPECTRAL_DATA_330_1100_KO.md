# OCRT C 2026-08-14 분광자료 통합

버전: `OCRT-v1.2-2026-08-14-KST-spectral-data-330-1100-component-gated`

- 196개 paired runtime 자료 적용
- 330–1100 nm strict data contract 적용
- PLOPS absorption의 800–1100 explicit zero 사용
- AHN 4종 scalar/vector 확장자료 적용
- atmospheric aerosol 176종 적용
- six-gas 330–1100 xsec 적용
- pure-water/psi range validation
- native coupled LUT solve-once 구조 유지
- 조건부 detritus 330–1100 candidate는 L=200/SOS acceptance 실패로 비활성 archive
- active detritus는 350–850 nm이며 Chl>0 범위 밖 요청은 fail-loud

상세: `docs/spectral_330_1100_20260814/README_FIRST_KO.md`
