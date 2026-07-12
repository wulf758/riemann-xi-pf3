# H1314 Independent Uniform Pair Cap Audit

Classification: `h1314_independent_uniform_pair_cap_audit_validated`

## Verdict

Independent audit validates the H943 uniform pair reduction

```text
R_r(t) <= 51/100        for the exact gamma-log family, r >= 5/2,
```

assuming the displayed parameter certificates, which were themselves
independently reduced to exact polynomial positivity. No counterexample found.

## Identities Checked

The dimensionless identities are correct:

```text
K_e = E[(cosh(zX)-1)/(z^2 X)],
K_o = E[(sinh(zX)-zX)/(z^2 X)],
O   = s E[(sinh(zX)-zX)/X^2],
R   = K_e - K_o tanh(O).
```

Global Taylor sign bounds:

```text
K_e - 1/2 <= z^2 E[X^3]/24 + z^4 E[X^5 exp(zX)]/720,
K_o       >= z E[X^2]/6,
O         >= s z^3/6.
```

## Regime `0 <= z <= 1`

With the constant package

```text
s >= 18,  delta2 <= 3/20,  delta3/delta2 <= 7/2,  E[X^5] <= 16/5,
tilted rate <= 21/8,
```

one gets

```text
R - 1/2 <= U(z) = 61 z^2/960 + (16/5) z^4 exp(21z/8)/720 - z tanh(3z^3)/6.
```

A 256-bit Arb monotone-endpoint cover on 1000 rational `z` boxes certifies

```text
U <= 0.002089262033 < 0.01;
```

refinement to 100000 boxes gives `0.002056446373`.

## Regime `z >= 1`

`K_- <= 1/2` and H943-B give

```text
R - 1/2 <= 0.5 exp(E - 2s(sinh z - z)).
```

`E(1) <= 7/4` and `2s(sinh 1 - 1) >= 6` give the junction value

```text
0.5 exp(-17/4) = 0.007132116955 < 0.01.
```

The exponent is decreasing: its positive derivative part is at most
`exp(z/4) + 3/4` while `2s(cosh z - 1) >= 36 exp(z)/8`.

## Parameter Certificates

- For `r >= 3`: all four numerator polynomials have positive coefficients
  after the shift `r = 3 + x`, `P = 60 + y`.
- On `[5/2, 3]`: all bivariate Bernstein coefficients are positive on
  `P in [153/4, 64]`.
- Tilted rate follows from `V8/V7 < q + 7` and `q <= r + 1`.

## Adversarial Scan

```text
15,008,001 points on r in [2.5, 10], z in [1e-4, 1]:
max actual R = 0.50160716353 at r = 2.5, z ~= 0.22972342.
Selected tails through r, z = 1e6: no issue.
```

## Budget Consequences (with H1313)

```text
E[W^2 K]         <= 0.5304 < 1,
E[(W^2+2) W^2 K] <= 2.8458 < 3.
```

No source files were edited by this audit.

## Recommendation Recorded

Promote the uniform cap `R <= 51/100` to a proved gamma-log lemma once the
exact polynomial and Arb verifier is persisted for reproducibility. Keep the
H943-to-Xi transfer and RH explicitly open.

## Provenance

Materialized on 2026-07-09 from AXMEM card
`h1314-independent-uniform-pair-cap-audit`. Numbers transcribed from the
card; verifier scripts not yet persisted.
