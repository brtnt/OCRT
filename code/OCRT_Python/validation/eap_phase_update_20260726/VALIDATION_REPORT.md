# OCRT EAP phase-function update validation

Date: 2026-07-26  
C version: `OCRT-v1.2-2026-07-26-KST-eap-phase-update`  
Python version: `pyOCRT-v1.2-2026-07-26-eap-phase-update`

## 1. Integrated scope

The supplied coated-sphere EAP implementation, frozen 17-species catalog, raw
source tables, generator, and representative runtime `.mie` files were merged
into both OCRT implementations.

The obsolete runtime files were deleted from both packages:

- `pico_Synechococcus_EAP.mie`
- `nano_Haptophytes_EAP.mie`
- `Diatoms_centric_EAP.mie`

The user-facing strings `pico`, `nano`, and `micro` remain accepted as
compatibility aliases, but now resolve to regenerated catalog files:

| Alias | New canonical species | Runtime file |
|---|---|---|
| pico | `eap_synechococcus` | `EAP_15_Synechococcus_D1p2.mie` |
| nano | `eap_hapto_prymnesiaceae` | `EAP_12_Hapto_Prymnesiaceae_D4.mie` |
| micro | `eap_diatoms_centric` | `EAP_02_Diatoms_centric_D6.mie` |

Seventeen canonical `eap_*` names are available in C and Python. The packages
bundle one representative diameter for each species. The included generator
can additionally create the Table-1 set of 70 species/diameter combinations;
those optional 70 files are not pre-bundled.

## 2. Runtime data validation

All 17 representative files passed the following checks:

- 101 spectral wavelengths, 350–850 nm at 5 nm spacing;
- 12 phase wavelengths;
- 361 phase angles;
- finite spectral, P11, P12, and P33 values;
- P11 normalization maximum absolute error: `1.7033e-5`;
- maximum `|P12|/P11`: `0.988732`;
- maximum `|P33|/P11`: `1.000000`;
- 555-nm zeroth Greek/Legendre coefficient error: `0.0` in the Python runtime transform;
- C and Python copies of every `.mie` file are byte-identical.

The four reference `.mie` samples supplied in the update ZIP are byte-identical
to the packaged regenerated files.

## 3. Phase smoothness

For the centric-diatom representative at 443 nm over 90–170 degrees, a local
normalized second-difference metric was evaluated.

| Data | Median roughness | Maximum roughness |
|---|---:|---:|
| Obsolete `Diatoms_centric_EAP.mie` | 2.10875e-2 | 1.04806 |
| New `EAP_02_Diatoms_centric_D6.mie` | 1.43949e-4 | 6.52819e-3 |

The new phase function reduces this metric by approximately 146x at the median
and 161x at the maximum.

## 4. C integration validation

- Release build: PASS.
- Public `rt_eap_mie_phase_compute` API test: PASS.
- 20-name catalog integration test (3 aliases + 17 canonical names): PASS.
- Organic Chl/phase-cache smoke suite: 10/10 PASS.
- Legacy alias vs canonical full RT output:
  - `pico` vs `eap_synechococcus`: stdout byte-identical;
  - `nano` vs `eap_hapto_prymnesiaceae`: stdout byte-identical;
  - `micro` vs `eap_diatoms_centric`: stdout byte-identical.
- Independent `eap_microcystis` full water RT case at 555 nm: PASS,
  `orders=57`, `converged=1`.
- EAP writer rebuilt with deterministic OpenMP radius-bin parallelism.
  A regenerated centric-diatom D6 file is byte-identical to the packaged file.

Core regression results after updating EAP-dependent anchors:

- water internal-reflection m-sign regression: PASS;
- Rrs/rrs I/Q/U output regression: PASS;
- IOP/Kd full-grid CSV regression: PASS;
- phase-wavelength PCHIP regression: PASS;
- external bottom-source bounds/shape regression: PASS;
- public RAA convention regression: PASS.

The IOP/Kd golden values changed intentionally because legacy `micro` now uses
the regenerated centric-diatom EAP data. The updated mixed-water reference has
`Kd0minus=0.3626284561 m^-1`.

## 5. Python integration validation

- Python compileall: PASS.
- Pytest: 11/11 PASS.
- All 20 accepted phyto-group strings load successfully and return finite,
  positive 555-nm total `a`, `b`, and `bb`.
- All 17 `.mie` files pass shape, normalization, and Mueller physicality tests.
- Alias and canonical files/moments are exactly identical.
- Single/full-grid CLI choice lists include all 17 canonical names.

## 6. Generator changes

The supplied physics and catalog were retained. Two implementation-only
optimizations were made to the build-time writer:

1. spectral-table wavelengths use an efficiency-only coated-Mie path and skip
   the unused 361-angle amplitude calculation;
2. independent radius bins are calculated with OpenMP, followed by summation in
   the original increasing-radius order.

The second step preserves deterministic floating-point accumulation. Output was
verified byte-identical to the supplied samples.

## 7. Supporting files

- `eap_species_inventory.csv`
- `alias_equivalence.csv`
- `python_group_iop_555nm.csv`
- `phase_smoothness.csv`
- `validation_summary.json`
- C and Python test logs under each package's
  `validation/eap_phase_update_20260726/` directory.
