# H1319 / H920 Verified Piecewise Closure

Classification: `h1319_piecewise_h920_closure_from_verified_dependencies`

## Result

The compact full-Xi assumption used by the earlier conditional assembly is
now discharged.  The canonical commands

```text
python tools/rh_h1319_full_xi_compact_interval_canonical_certificate.py
python tools/rh_h1319_piecewise_h920_conditional_assembly.py --require-unconditional
```

both exit with status `0`.  The second command reports

```text
compact_assumption_discharged = true
all_required_dependencies_pass = true
```

## Curvature chain

The verified full-Xi bound

```text
t(r)^2 kappa_Xi'''(t(r)) >= -3/4
```

is covered without a gap by

```text
r(125) -> 4 -> 5 -> 10 -> 22 -> infinity.
```

The five pieces are:

1. quotient-free Arb boxes on `r(125) <= r <= 4`;
2. score/mean-value Arb boxes on `4 <= r <= 5`;
3. score/mean-value Arb boxes on `5 <= r <= 10`;
4. score/mean-value Arb boxes on `10 <= r <= 22`;
5. an exact analytic bound for `r >= 22`.

The three middle reports contain `360` exact rational boxes, no unresolved
leaf, and exact endpoint joins.  The independent audit uses a 140-digit
mpmath oracle with explicit conversion radii and passes all `180/180`
value/derivative comparisons, the four symbolic identities, and independent
recomputations of the three worst boxes.

## H920 / integer assembly

The already-audited exact arithmetic gives

```text
2 <= n <= 126:  finite H803 certificate;
127 <= n <= 566: compact coefficient a=3/4;
n >= 567: glued coefficient a=141/142.
```

The exact geometric thresholds are

```text
N_geom(3/4) = 15,
N_geom(141/142) = 567.
```

Consequently the H908-A unit barrier holds for every `n >= 127`, and the
exact H803 factorial barrier holds for every `n >= 2`.  The stronger unit
barrier is not claimed for `n=2,...,7`, where the finite audit records that it
is false.

## Canonical evidence

- `research/riemann/h1319_full_xi_compact_interval_certificate.json`
- `research/riemann/h1319_score_mean_value_independent_audit_canonical.json`
- `research/riemann/h1319_full_xi_high_band_r22_analytic_certificate.json`
- `research/riemann/h1319_piecewise_h920_conditional_assembly.json`

The compact certificate pins every component, the mean-value engine, H1272,
the high-precision audit runner, and all three middle reports by SHA-256.

## Scope

This closes the H1319 compact curvature obligation and the H920/H908-to-H803
piecewise obligation represented by these artifacts.  It does **not** close
the remaining H907 all-degree Jensen/total-positivity step and therefore is
not a proof of the Riemann Hypothesis.
