# H1315 Gamma-Log Moment-Box Global Assembly

Classification: `h1315_gamma_log_moment_box_global_assembly_proved_from_audited_dependencies`

## Result

For the dominant gamma-log family with `alpha=9/4`, the normalized saddle
moment box

```text
0 < a <= 1,
-a <= E[W] <= 0,
E[W^3] >= -3a
```

holds for every saddle with `r=r_q>=2`. Consequently,

```text
K_alpha'''(q) >= -35/(q^2 r_q).                         (1)
```

This is a global assembly of H941, H1314, H1313, H925, H922a, H927, and
H922. It closes the dominant gamma-log moment-box target. It does not perform
the transfer to the full Xi kernel and does not prove the Riemann Hypothesis.

The autonomous dependency verifier is

```text
tools/rh_h1315_gamma_log_moment_box_global_assembly.py
```

Its canonical output is

```text
research/riemann/h1315_gamma_log_moment_box_global_assembly.json
```

## Setup

With

```text
V_alpha(z) = pi exp(exp(z))-alpha exp(z)-z,
V_alpha'(z_q) = q,
r_q = exp(z_q),
A = V_alpha''(z_q),
B = V_alpha'''(z_q),
W = sqrt(A)(Z-z_q),
a = B/A^(3/2),
```

H925 defines the nonnegative score excess

```text
N(w)=Psi'(w)-w
```

and proves the exact Stein identities

```text
E[W]   = -E[N(W)],
E[W^3] = -E[(W^2+2)N(W)].                              (2)
```

Thus the two lower signs in the moment box follow from

```text
E[N(W)] <= a,
E[(W^2+2)N(W)] <= 3a.                                 (3)
```

## Low Band: `2 <= r <= 5/2`

H941's rigorous score-adaptive cover proves the two H925 lower-sign budgets
on the entire closed interval

```text
2 <= r <= 2.50.
```

This supplies (3), and hence the two lower moment bounds through (2), on the
low band.

## High Band: `r >= 5/2`

H1314 proves and independently audits the uniform corrected-pair cap

```text
K <= 51/100
```

in the paired expectation representation, equivalently the bound denoted
`R_r(t)<=51/100` in H1314. H1313 proves globally on this band that

```text
E[W^2] <= 26/25,
E[W^4] <= 7/2,
E[W^6] <= 21.                                          (4)
```

Only the first two bounds in (4) are needed for the uniform-cap route. Exact
rational arithmetic gives

```text
E[N]/a = E[W^2 K]
       <= (51/100)(26/25)
        = 663/1250
        = 0.5304
        < 1,                                            (5)
```

and

```text
E[(W^2+2)N]/a = E[(W^2+2)W^2 K]
               <= (51/100)(7/2 + 2(26/25))
                = 14229/5000
                = 2.8458
                < 3.                                    (6)
```

Equations (5)-(6) prove (3) for every `r>=5/2`.

The H1315 verifier relaunches H1313's compact Arb/dependency certificate. It
attests the published H1314 statement and arithmetic, but it does not recreate
H1314's polynomial/Arb proof. Consequently this assembly preserves H1314's
published status `proved_and_independently_audited_pending_verifier_persistence`;
it does not silently upgrade that dependency.

## Gap-Free Global Moment Box

The two bands

```text
[2,5/2],
[5/2,infinity)
```

meet at the exact endpoint `5/2`, so there is no uncovered interval. H925 and
the preceding budgets yield

```text
E[W] >= -a,
E[W^3] >= -3a.
```

H922a supplies the global upper signs

```text
E[W] <= 0,
E[W^3] <= 0,
```

and H927 supplies

```text
0<a<=1
```

for every `r>=2`. This proves the full H922 moment box globally.

## Cumulant Consequence

H922 gives the exact normalized identity

```text
K_alpha'''(q)=A^(-3/2) E[(W-EW)^3].
```

The moment box implies

```text
E[(W-EW)^3] >= -5a,
```

so, because `a=B/A^(3/2)`,

```text
K_alpha'''(q) >= -5B/A^3.
```

The coefficient estimate consumed by H922 is

```text
B/A^3 <= 7/(q^2 r_q)
```

for `alpha=9/4` and `r_q>=2`. Multiplying the exact constants `5*7=35`
proves (1).

## Machine-Checked Assembly

The H1315 verifier performs the following checks:

```text
1. rerun H1313 and require all_checks_pass=true;
2. attest the required theorem fragments and hashes for H1314, H941, H925,
   H922a, H927, and H922;
3. compute 663/1250 and 14229/5000 with exact Fraction arithmetic;
4. verify both strict budget inequalities;
5. verify the gap-free band seam at 5/2;
6. assemble the moment box and the exact final constant 35.
```

The canonical JSON must report

```text
all_checks_pass = true.
```

## Exact Scope

Closed by H1315:

```text
the dominant gamma-log normalized moment box for every r_q>=2;
the resulting one-sided bound K_alpha'''(q)>=-35/(q^2 r_q).
```

Not closed by H1315:

```text
effective derivative constants for the full-Xi correction;
finite H920 closure;
all-degree Jensen hyperbolicity or total positivity;
the Riemann Hypothesis.
```
