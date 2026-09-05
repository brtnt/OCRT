# OCRT 330–1100 nm 최종 acceptance test plan

## 1. 실행 순서

1. payload SHA-256 검사
2. C/Python paired data install
3. source-range and fallback scan
4. C/Python build
5. static data validators
6. component-only IOP tests
7. C/Python matched-input RT
8. overlap regression
9. full 330–1100 spectral sweep
10. final manifest and release note

## 2. Mandatory wavelength set

```text
330, 331, 340, 349, 350, 351, 400, 412, 443, 490, 510, 555,
620, 660, 680, 709, 745, 750, 760, 799, 800, 801, 820, 865,
900, 935, 940, 1000, 1030, 1099, 1100 nm
```

## 3. Mandatory component cases

- pure water at T=0/5/10/20/25/30°C
- Chl-only at Chl=0.1/1/3 mg m-3
- detritus-only
- AHN brown/yellow/red/calcareous individually
- representative aerosol from each family and all-file load test
- each gas individually, all gases, gas-off
- Rayleigh-only
- black Fresnel ocean / glint-decoupled surface

## 4. Output checks

- input IOP breakdown: a, b, bb and component terms
- P11/P12/P33, g, bb/b
- TOA I/Q/U
- Rrs and rrs I/Q/U
- Ed0+/Ed0-/Eu0- direct/diffuse where available
- convergence/order/layer diagnostics

## 5. Pass criteria

- no NaN/Inf
- no unsupported silent clamp
- C/Python file hashes exact
- C/Python numerical parity at roundoff level
- Mie physicality pass
- 20°C pure-water baseline exact
- gas 350–1100 canonical exact
- explicit zero policies exact
- component changes correctly classified as intentional or regression
- all target wavelengths complete without range error
