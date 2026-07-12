# H1319 Piecewise H920 Conditional Assembly

Classification: `h1319_piecewise_h920_closure_conditional_on_compact_certificate`

## Result

The off-by-one and transition problem is closed **conditionally**.  Assume:

1. the existing H801/H803/H804 finite certificate supplies the exact H803
   factorial threshold for every integer `2 <= n <= 126`;
2. the full Xi cumulant satisfies

   ```text
   kappa_Xi'''(t) >= -(3/4)/t^2       for 125 <= t <= T_71;
   ```

3. H1317 supplies

   ```text
   kappa_Xi'''(t) >= -(141/142)/t^2   for t >= T_71,
   T_71=(71*pi*exp(71)-(9/4)*71-1)/2.
   ```

Then the H908-A unit barrier holds for every integer `n >= 127`, and the
weaker exact H803 factorial barrier holds for every integer `n >= 2`.

The only new assumption not already represented by a passing dependency is
item 2.  No compatible compact interval certificate exists at the time of this
report, so this is not yet an unconditional closure.

## 1. Exact Cube Lemma

H908 gives

```text
Delta^3 kappa(n-2)
 = integral_[0,1]^3 kappa'''(n-2+r1+r2+r3) dr1 dr2 dr3.
```

The integration cube covers the closed interval

```text
[n-2,n+1].
```

If `kappa'''(t) >= -a/t^2` throughout that cube, then

```text
Delta^3 kappa(n-2) >= -a/(n-2)^2.
```

The exact geometric condition needed to compare this with the unit barrier is

```text
a*n^2 <= (n-2)^2.                              (G_a)
```

Indeed, under `(G_a)`,

```text
Delta^3 kappa(n-2)
 >= -1/n^2
 >= log(1-1/n^2).
```

Therefore

```text
exp(Delta^3 kappa(n-2)) >= 1-1/n^2.
```

The verifier searches the first integer satisfying `(G_a)` with exact rational
arithmetic; no floating-point square root or rounded ceiling is used.

## 2. The Two Exact Geometric Thresholds

For the compact coefficient,

```text
a=3/4:
  (3/4)*14^2 > 12^2,
  (3/4)*15^2 <= 13^2,
  N_geom(3/4)=15.
```

For the H1317 coefficient,

```text
a=141/142:
  (141/142)*566^2 > 564^2,
  (141/142)*567^2 <= 565^2,
  N_geom(141/142)=567.
```

Thus `567` is not a loose rounded threshold: it is the first valid integer for
`141/142`.

## 3. Piecewise Integer Cover

| integer range | cube and pointwise input | conclusion |
| --- | --- | --- |
| `2 <= n <= 126` | existing finite H803 certificate | exact H803 threshold; no unit-barrier claim |
| `127 <= n <= 566` | `[n-2,n+1] subset [125,567] subset [125,T_71]`, use `a=3/4` | H908-A unit barrier, hence H803 |
| `n >= 567` | `[n-2,n+1] subset [565,infinity)`, use the glued global `a=141/142` bound | H908-A unit barrier, hence H803 |

Both endpoints in the compact hypothesis matter.  In particular, the first
analytic cube is

```text
n=127: [125,128],
```

so a statement beginning only at `t>125` would leave an endpoint obligation.

The last cube assigned solely to the compact block is

```text
n=566: [564,567].
```

It lies below `T_71` by a very coarse exact estimate.  From `pi>3` and `e>2`,

```text
T_71 > 2011728121702468861033853/8 > 567.
```

## 4. Why No Mixed Cube Is Lost

On the compact interval, `3/4 < 141/142`, hence

```text
-(3/4)/t^2 >= -(141/142)/t^2.
```

Combining the compact estimate on `[125,T_71]` with H1317 on
`[T_71,infinity)` therefore yields the single weaker estimate

```text
kappa_Xi'''(t) >= -(141/142)/t^2
for every t>=125.                              (P)
```

For `n>=567`, the entire cube starts at `n-2>=565>=125`, and the exact
geometric threshold for `(P)` is already satisfied.  Consequently this one
global bound covers every cube that meets or crosses `T_71`; there is no need
to locate the enormous integer nearest `T_71`, and no floor/ceiling ambiguity
is introduced there.

## 5. Unit Barrier Versus the All-`n` H803 Conclusion

The finite data must not be described as an all-small-`n` unit-barrier proof.
H804 explicitly records failures of the stronger unit barrier at

```text
n=2,3,4,5,6,7.
```

What the finite certificate supplies is the weaker exact H803 threshold.  For
the analytic tail, the unit barrier implies H803 because, with
`C_n=(2n)(2n-1)` and
`theta_n=C_(n-1)C_(n+1)/C_n^2`,

```text
1-theta_n-1/n^2
 = 4(n^2-1)/(n^2(2n-1)^2)
 > 0                                             for n>=2.
```

The correct assembled conclusions are therefore:

```text
H908-A unit barrier:       every n>=127;
exact H803 factorial bound: every n>=2.
```

Claiming the unit barrier for every `n>=2` would be false.

## 6. Verifier And Current Status

Run:

```text
python tools/rh_h1319_piecewise_h920_conditional_assembly.py
```

The verifier:

- checks both geometric thresholds with `fractions.Fraction`;
- checks every endpoint `125`, `126`, `127`, `566`, and `567` used in the
  decomposition;
- verifies the exact coarse lower bound `T_71>567`;
- loads all `125` H804 rows for `2<=n<=126` and checks their inherited H801
  certificate flags;
- executes the H1317 verifier and requires its `all_checks_pass=true` result;
- emits
  `research/riemann/h1319_piecewise_h920_conditional_assembly.json`.

The current machine-readable result is:

```text
conditional_logic_checks_pass = true
compact_assumption_discharged = false
all_required_dependencies_pass = false
```

For a no-cheat CI check, run:

```text
python tools/rh_h1319_piecewise_h920_conditional_assembly.py --require-unconditional
```

That command must exit with status `2` until a compatible compact certificate
is present.

The expected default path is

```text
research/riemann/h1319_full_xi_compact_interval_certificate.json.
```

Its interface must identify the **full Xi** kernel, state the exact closed
`t`-domain `[125,T_71]` and lower bound `-3/4`, and affirm rigorous interval
arithmetic, full domain coverage, saddle-map coverage, and passing local lower
bounds.  The assembly verifier checks this interface only; the compact
verifier and its interval boxes remain the primary numerical proof evidence.

## 7. Scope

Even after the compact certificate is supplied, this assembly closes only the
H803/H908 moment-curvature obligation in the selected route.  It does not prove
the H907 all-degree Jensen/total-positivity implication and therefore does not
prove the Riemann hypothesis.

