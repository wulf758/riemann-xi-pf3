# H1314 Uniform Pair Cap: Proved And Independently Audited

Classification: `h1314_uniform_pair_cap_proved_audited_pending_verifier_persistence`

## Consolidated Statement

For the exact gamma-log model:

```text
R_r(t) <= 51/100        for all r >= 5/2, t >= 0.
```

Dimensional identities and normalization were verified by two independent
agents. Together with the H1313 moment bounds this closes both H943-C
gamma-log budgets:

```text
E[W^2 K]         <= 663/1250   < 1,
E[(W^2+2) W^2 K] <= 14229/5000 < 3.
```

An independent analytic proof also closes the strong H946 corrected-pair
bound `R_r(t) <= 1/2 + C_r t^2`; it is stored in
`h1314_h946_corrected_pair_global_proof.md` (AXMEM
`h1314-h946-corrected-pair-global-proof`).

## Certified Inputs

Parameters:

```text
s >= 18,  delta2 <= 3/20,  delta3 <= (7/2) delta2,  E[X^5] <= 16/5.
```

- `[5/2, 3]`: exact bivariate Bernstein coefficients on
  `r in [5/2,3], P in [153/4, 64]`, all positive.
- `r >= 3`: positive coefficients after the shift `r = 3+X`, `P = 60+Y`.
- Tilted growth `E[X^5 exp(zX)] <= 4 exp(21z/8)` via positivity of
  `(q+7) T7 - T8` and `mu >= r + 3/2`.

## Two Regimes

For `0 <= z <= 1`:

```text
R - 1/2 <= U4(z) = 61 z^2/960 + z^4 exp(21z/8)/180 - (z/6) tanh(3z^3),
```

with the five-interval rational proof (`h943c_u4_five_interval_proof.md`)
using `tanh u >= u/sqrt(1+u^2)`; rational maxima

```text
1687/233280, 431/233280, -17/29952, -489/14080, 959/112500,
```

all `< 1/100`.

For `z >= 1`:

```text
R - 1/2 <= 0.5 exp(H(z)),   H(1) <= -17/4,   H' < 0,
```

hence `<= 0.0071322 < 0.01`.

## Adversarial Evidence

Scan of 15,008,001 points: max actual value

```text
R = 0.50160716353 at the corner r = 2.5,
```

no counterexample. See also the wide search `2.5 <= r <= 500`,
`0.001 <= z <= 20` in `h943c_uniform_51over100_route_validated.md`.

## Non-Claims

Does not prove the Xi transfer. Does not prove RH.

## Pending For Publication Status

```text
1. persist the exact polynomial / Arb verifier scripts (reproducibility);
2. explicit continuity statement at z = 0;
3. the convexity statement used for O;
4. rational bounds for pi*exp(r).
```

## Provenance

Materialized on 2026-07-09 from AXMEM cards
`h1314-uniform-pair-cap-proved-audited`,
`h1314-uniform-pair-cap-readonly-certificate`,
`h1314-independent-uniform-pair-cap-audit`,
`h943c-uniform-51over100-route-validated`,
`h943c-u4-five-interval-proof`.
No H1314 files were published at proof time because `apply_patch` was
rejected by an app usage quota; this note is the parent-side integration.
The numbers are transcriptions from the cards, not regenerated outputs; the
verifier scripts must be re-run once persisted to produce canonical JSON
artifacts.
