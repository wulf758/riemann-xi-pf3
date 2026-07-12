# H1314 Uniform Pair Cap Read-Only Certificate

Classification: `h1314_uniform_pair_cap_readonly_certificate`

## Statement

Candidate rigorous proof of the uniform bound

```text
R_r(t) <= 51/100        for all r >= 5/2, t >= 0,
```

which together with the H1313 moment bounds suffices for H943-C.

## Exact Identities

Under the mixture normalization

```text
nu_m = c_m m^2 / A,
mu = B/A,
X = M/mu,
s = 1/a^2 = A^3/B^2,
z = a t,
```

the pair ratio has the exact dimensionless form

```text
R = K_e - K_o tanh(O),
K_e = E[(cosh(zX)-1)/(z^2 X)],
K_o = E[(sinh(zX)-zX)/(z^2 X)],
O   = s E[(sinh(zX)-zX)/X^2].
```

## Regime `0 <= z <= 1`

Taylor bounds give

```text
R - 1/2 <= U(z) = 61 z^2/960 + (16/5) z^4 exp(21z/8)/720 - z tanh(3z^3)/6.
```

Arb certification, 4096 cells: zero failures, worst majorant

```text
0.0020641952 < 0.01.
```

## Parameter Certificates

On `r in [5/2, 3]`, 4096 Arb boxes on the numerators certify

```text
s >= 18,
delta2 <= 3/20,
delta3 <= (7/2) delta2,
E[X^5] <= 16/5,
```

with zero failures.

For `r >= 3`, exact certificates after the shift `X = r-3`, `Y = P-60`:
coefficient minima

```text
1, 3/20, 1/2, 21/10
```

for `s >= 30`, `delta2`, `delta3`, and `E[X^5] <= 31/10` respectively.

Tilted growth rate `<= 21/8` via

```text
(q+7) T7 - T8 = q (21 q^5 + 280 q^4 + 1050 q^3 + 1204 q^2 + 315 q + 6)
```

and `mu >= r + 3/2`.

## Regime `z >= 1`

```text
R - 1/2 <= 0.5 exp(E - 2s(sinh z - z)),
```

junction value `<= 0.5 exp(-17/4) = 0.00713212 < 0.01`, and the exponent has
negative derivative.

## Budget Consequences (with H1313)

```text
(51/100)(26/25)        = 663/1250   < 1,
(51/100)(7/2 + 52/25)  = 14229/5000 < 3.
```

## Status

Read-only candidate at recording time; subsequently validated by two
adversarial audits (`h1314_independent_uniform_pair_cap_audit.md`,
`h943c_uniform_51over100_route_validated.md`).

## Provenance

Materialized on 2026-07-09 from AXMEM card
`h1314-uniform-pair-cap-readonly-certificate`. `apply_patch` was blocked by an
app usage quota at recording time, so no H1314 file was created then. The
exact/Arb verifier scripts remain to be persisted; the numbers above are
transcribed from the card, not regenerated.
