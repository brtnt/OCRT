# Water IOP data

- `water_coef_z09_1nm.txt`: pure-water absorption/scattering LUT.
- `psi_T_rottgers2014_OCRT.txt`: approximate temperature correction; see header.
- `pico_Synechococcus_EAP.mie`, `nano_Haptophytes_EAP.mie`,
  `Diatoms_centric_EAP.mie`: OCRT phytoplankton P11/P12/P33.
- `Detritus_Stramski2001.mie`: organic detritus P11/P12/P33.
- `aph_ccrr_morel1988_mm01.txt`: canonical CCRR Chl adapter table.
- `aph_bricaud_1998.txt`: legacy filename compatibility copy, same bytes.
- `source/`: retained normalized source spectrum.

Organic phase files are valid for the documented 350–850 nm range. The CCRR source
itself marks 300–350 and 700–1000 nm as extrapolated.
