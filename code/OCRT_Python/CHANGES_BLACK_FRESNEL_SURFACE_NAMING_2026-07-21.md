# Black Fresnel ocean surface naming update

Version: `pyOCRT-v1.2-2026-07-21-black-fresnel-surface-naming`

- Canonical high-level surface mode: `black_fresnel_ocean`.
- The mode means a rough Fresnel air-water interface with zero water-leaving radiance.
- `coxmunk` remains a deprecated compatibility alias.
- The slope-variance model is separate: default `OCRT floor law` uses
  `sigma^2 = 0.003 + 0.00512 * max(0.01, W)`; optional type 0 is
  Nakajima-Tanaka. The default is therefore not described as an exact
  Cox-Munk implementation.
- Numerical algorithms and outputs are unchanged.
