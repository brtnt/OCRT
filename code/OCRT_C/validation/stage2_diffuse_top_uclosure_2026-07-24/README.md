# Stage-2 diffuse-top incoming-U closure validation

This directory records the validation for
`OCRT-v1.2-2026-07-24-KST-stage2-interface-fix123-dtp-uclosure`.

The production change negates only the molecular incoming-U coefficients
`ruI`, `ruQ`, and `ruU` inside `ocrt_add_diffuse_top_primary`. The existing
U-output wrapper sign is retained.

Acceptance evidence:

- reference-free source/operator closure for m=0..4, every incident node and
  I/Q/U basis;
- pure Rayleigh max absolute closure error: 2.776e-17;
- mixed Rayleigh-aerosol max absolute closure error: 5.551e-17;
- 72/72 physical runs and 600/600 matched cells;
- no runtime increase relative to FIX123;
- six independent pure-water Rrs/rrs I/Q/U scatterplots.

The attached other-session harness remains diagnostic because it does not pass
an identical physical water depth to both codes and its 412-nm water SOS run
can stop at order 20. The production 72-run gate uses the prior exact input
contract and raises the water order cap explicitly.
