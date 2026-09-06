# Water IOP data

- `water_coef_z09_1nm.txt`: pure-water absorption / scattering table.
- `psi_T_rottgers2014_OCRT.txt`: approximate temperature correction; see the header.
- `pico_Synechococcus_EAP.mie`, `nano_Haptophytes_EAP.mie`, `Diatoms_centric_EAP.mie`: OCRT
  phytoplankton P11/P12/P33 (legacy names).
- `Detritus_Stramski2001.mie`: organic detritus P11/P12/P33.
- `aph_ccrr_morel1988_mm01.txt`: canonical CCRR Chl adapter table.
- `aph_bricaud_1998.txt`: legacy file-name compatibility copy, same bytes.
- `source/`: retained normalized source spectrum.

Organic phase files are valid over the documented 350–850 nm range. The CCRR source marks
300–350 nm and 700–1000 nm as extrapolated.

---

# 해수 IOP 자료

- `water_coef_z09_1nm.txt`: 순수해수 흡수/산란 표.
- `psi_T_rottgers2014_OCRT.txt`: 근사 온도 보정. 머리글 참조.
- `pico_Synechococcus_EAP.mie`, `nano_Haptophytes_EAP.mie`, `Diatoms_centric_EAP.mie`: OCRT
  식물플랑크톤 P11/P12/P33(옛 이름).
- `Detritus_Stramski2001.mie`: 유기 쇄설물 P11/P12/P33.
- `aph_ccrr_morel1988_mm01.txt`: 표준 CCRR Chl 어댑터 표.
- `aph_bricaud_1998.txt`: 옛 파일 이름 호환 사본, 내용 동일.
- `source/`: 보존된 정규화 원 스펙트럼.

유기물 위상 파일은 문서화된 350–850 nm 범위에서 유효하다. CCRR 원자료는 300–350 nm 와 700–1000 nm 를
외삽으로 표시한다.
