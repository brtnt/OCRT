# OCRT EAP phytoplankton option guide

## Default behavior

- The OCRT Chl model default is `micro`.
- `micro` maps to canonical `eap_diatoms_centric` and file `EAP_02_Diatoms_centric_D6.mie`.
- Omitting the phase option is therefore the recommended normal-mode usage.
- Explicit `micro` and `eap_diatoms_centric` are both allowed without advanced mode and are byte-identical.

## Advanced selections

All other choices require `OCRT_ADVANCED=1`.

| option | representative file | mode |
|---|---|---|
| `micro` | `EAP_02_Diatoms_centric_D6.mie` | standard default |
| `eap_diatoms_centric` | `EAP_02_Diatoms_centric_D6.mie` | standard equivalent |
| `pico` | `EAP_15_Synechococcus_D1p2.mie` | advanced alias |
| `nano` | `EAP_12_Hapto_Prymnesiaceae_D4.mie` | advanced alias |
| `eap_diatoms_pennate` | `EAP_00_Diatoms_pennate_D6.mie` | advanced |
| `eap_chlorophytes` | `EAP_01_Chlorophytes_D8.mie` | advanced |
| `eap_cryptophytes` | `EAP_03_Cryptophytes_D6.mie` | advanced |
| `eap_cyano_blue` | `EAP_04_Cyano_blue_D6.mie` | advanced |
| `eap_cyano_red` | `EAP_05_Cyano_red_D6.mie` | advanced |
| `eap_dinoflagellates` | `EAP_06_Dinoflagellates_D24.mie` | advanced |
| `eap_eustigmatophytes` | `EAP_07_Eustigmatophytes_D6.mie` | advanced |
| `eap_hapto_pavlovaceae` | `EAP_08_Hapto_Pavlovaceae_D6.mie` | advanced |
| `eap_pelagophytes` | `EAP_09_Pelagophytes_D3.mie` | advanced |
| `eap_prasinophytes` | `EAP_10_Prasinophytes_D3.mie` | advanced |
| `eap_prochlorococcus` | `EAP_11_Prochlorococcus_D0p5.mie` | advanced |
| `eap_hapto_prymnesiaceae` | `EAP_12_Hapto_Prymnesiaceae_D4.mie` | advanced |
| `eap_raphidophytes` | `EAP_13_Raphidophytes_D24.mie` | advanced |
| `eap_rhodophytes` | `EAP_14_Rhodophytes_D6.mie` | advanced |
| `eap_synechococcus` | `EAP_15_Synechococcus_D1p2.mie` | advanced |
| `eap_microcystis` | `EAP_16_Microcystis_D5.mie` | advanced |

## Examples

Normal mode (no phase option required):

```bash
./build/ocrt --surface ocean --water-model ocrt \
  --ocrt-chl 1 --ocrt-tsm 0 --ocrt-adom440 0 ...
```

Advanced C selection:

```bash
OCRT_ADVANCED=1 ./build/ocrt ... \
  --ocrt-phyto-group eap_synechococcus
```

Advanced Python selection:

```bash
OCRT_ADVANCED=1 python ocrt_solve.py single ... \
  --phyto-group eap_synechococcus
```

Windows CMD:

```bat
set OCRT_ADVANCED=1
build\ocrt.exe ... --ocrt-phyto-group eap_synechococcus
```

## Compatibility note

- Default command syntax is unchanged from the previous package.
- Numerical results differ from the obsolete EAP files because the phase data were intentionally regenerated.
- Existing scripts that explicitly use `pico` or `nano` now need `OCRT_ADVANCED=1`.
