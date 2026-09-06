# A2010ver — Ahmad et al. (2010) 80-model Mie family

This set contains 80 OCRT-compatible aerosol Mie tables that follow the Ahmad et al. (2010)
paper. Every file name carries the `A2010ver_` prefix.

Examples:

- `A2010ver_r30f95v01.mie`: RH 30 %, fine-mode volume fraction 0.95
- `A2010ver_r80f20v01.mie`: RH 80 %, fine-mode volume fraction 0.20
- `A2010ver_r95f00v01.mie`: RH 95 %, coarse-only end point

Numerical format (as received):

- 20 wavelengths, 0.350–3.750 µm
- 361 scattering angles, 180° to 0° at 0.5° spacing
- P11, P12 and P33 phase-matrix blocks

Scientific boundary:

- 0.412–0.865 µm is constrained by Ahmad et al. (2010), Tables 3 and 4.
- Outside that interval the package uses the documented constituent-level engineering extension
  from the earlier OCRT paper-centric reconstruction.
- This family is intentionally distinct from the simplified AccuRT effective-mode family.

---

# A2010ver — Ahmad et al. (2010) 80 모델 Mie 계열

이 세트는 Ahmad et al. (2010) 논문을 따르는 OCRT 호환 에어로졸 Mie 표 80개다. 모든 파일 이름에
`A2010ver_` 접두사가 붙는다.

예:

- `A2010ver_r30f95v01.mie`: 상대습도 30 %, 미세모드 부피비 0.95
- `A2010ver_r80f20v01.mie`: 상대습도 80 %, 미세모드 부피비 0.20
- `A2010ver_r95f00v01.mie`: 상대습도 95 %, 조대모드만 있는 끝점

수치 형식(수령 당시):

- 파장 20개, 0.350–3.750 µm
- 산란각 361개, 180° 에서 0° 까지 0.5° 간격
- P11, P12, P33 위상행렬 블록

과학적 범위:

- 0.412–0.865 µm 는 Ahmad et al. (2010) 표 3·4 로 제약된다.
- 그 밖의 구간은 이전 OCRT 논문 기준 재구성에서 문서화한 성분 수준 공학적 확장을 쓴다.
- 이 계열은 단순화된 AccuRT 유효모드 계열과 의도적으로 구분된다.
