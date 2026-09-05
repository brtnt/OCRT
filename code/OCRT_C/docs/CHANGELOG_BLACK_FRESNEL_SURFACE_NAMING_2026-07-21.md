# Black Fresnel ocean surface naming update

Version: `OCRT-v1.2-2026-07-21-KST-black-fresnel-surface-naming`

## Decision

The former high-level surface option `coxmunk` mixed two independent concepts:

1. the physical boundary condition: a rough Fresnel air-water interface with a
   perfectly black ocean below it, so no water-leaving radiance is added; and
2. the slope-variance parameterization used by the rough-surface kernel.

The canonical option is now:

```text
--surface black_fresnel_ocean
```

`--surface coxmunk` remains accepted as a deprecated compatibility alias and
produces identical output.

## Slope variance

The default model is named `ocrt-floor`:

```text
sigma^2 = 0.003 + 0.00512 * max(0.01, W)
```

It is CM54-derived but includes the OCRT low-wind floor and is therefore not
described as an unmodified Cox-Munk implementation. The optional alternative is
`nakajima-tanaka`. New string interface:

```text
--sigma-model ocrt-floor
--sigma-model nakajima-tanaka
```

The legacy numeric `--sigma-type` selector remains supported.

## Numerical impact

None. The update changes names, parser aliases, enum names, documentation, and
API labels only. Regression tests confirm byte-identical stdout for the old and
new surface names under the same physical options.
