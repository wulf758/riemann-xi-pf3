# H920 Effective Finite Closure Gate

Classification: `h920_effective_finite_closure_reduction_proved`

## Question

Can the non-effective eventual H908-A result from H919 be converted into a
concrete finite verification target?

## Short Answer

Yes. H919 leaves exactly one effective constant gap.

It is enough to prove an explicit pointwise third-cumulant lower bound

```text
kappa'''(t) >= -a/t^2
```

for all `t >= T`, with some explicit constant `0 < a < 1`. Then H908-A holds
for every integer

```text
n >= N_eff(a,T)
  = ceil(max(3, T+2, 2/(1-sqrt(a)))).
```

Consequently H803 is reduced to a finite exact check for

```text
2 <= n < N_eff(a,T).
```

The existing H803/H804 finite data already covers `2 <= n <= 126`, so any
effective constant extraction giving `N_eff(a,T) <= 127` closes the H803
moment-curvature lemma immediately at the current report level. If
`N_eff(a,T) > 127`, the remaining work is a finite Arb-style extension up to
`N_eff(a,T)-1`.

This does not prove RH. It would close only the H803/H908 moment-curvature
sublemma; H907's global all-degree Jensen/TP assembly would still remain.

## Effective Closure Lemma

Let `kappa(t)=log M_t` be the log-Laplace transform of the Xi log-height
moment measure from H908. Assume there are explicit constants

```text
0 < a < 1,
T >= 0
```

such that

```text
(A_{a,T})    kappa'''(t) >= -a/t^2
```

for every `t >= T`.

Define

```text
N_eff(a,T) = ceil(max(3, T+2, 2/(1-sqrt(a)))).
```

Then for every integer `n >= N_eff(a,T)`,

```text
Delta^3 kappa(n-2) >= log(1 - 1/n^2).
```

Thus the H804 unit barrier holds for every `n >= N_eff(a,T)`.

## Proof

H908 gives the cube-average identity

```text
Delta^3 kappa(n-2)
 = int_[0,1]^3 kappa'''(n-2+r1+r2+r3) dr1 dr2 dr3.
```

If `n >= T+2`, then every point in the integration cube satisfies

```text
t = n-2+r1+r2+r3 >= T.
```

Using `(A_{a,T})` and `t >= n-2`,

```text
kappa'''(t) >= -a/t^2 >= -a/(n-2)^2.
```

Therefore

```text
Delta^3 kappa(n-2) >= -a/(n-2)^2.
```

The condition

```text
n >= 2/(1-sqrt(a))
```

is equivalent to

```text
sqrt(a) n <= n-2,
```

and hence to

```text
a/(n-2)^2 <= 1/n^2.
```

So for `n >= N_eff(a,T)`,

```text
Delta^3 kappa(n-2) >= -1/n^2.
```

Finally, for `0 < x < 1`,

```text
log(1-x) <= -x.
```

Taking `x=1/n^2`,

```text
-1/n^2 >= log(1 - 1/n^2).
```

Combining the inequalities proves

```text
Delta^3 kappa(n-2) >= log(1 - 1/n^2).
```

This is exactly the H908-A unit-barrier condition.

## Consequence For H803

H908 proves

```text
T_{n+1}/T_n = exp(Delta^3 kappa(n-2)).
```

Thus H908-A gives

```text
T_{n+1}/T_n >= 1 - 1/n^2.
```

H804 proves that this stronger unit barrier implies the exact H803 threshold

```text
T_{n+1}/T_n >= C_{n-1}C_{n+1}/C_n^2,
C_n=(2n)(2n-1).
```

Therefore, under `(A_{a,T})`, H803 holds for all `n >= N_eff(a,T)`.

The finite remainder is exactly

```text
2 <= n < N_eff(a,T).
```

## Useful Constants

The geometric part of the threshold is

```text
N_geom(a) = ceil(2/(1-sqrt(a))).
```

Examples:

| `a` | `N_geom(a)` |
| ---: | ---: |
| `3/4` | `15` |
| `2/3` | `11` |
| `1/2` | `7` |
| `1/3` | `5` |
| `1/4` | `4` |

H919 only used `a=3/4` as a convenient non-effective eventual constant. Since
the asymptotic coefficient in H919 is

```text
1/(2U_t) * (1 + O(1/U_t)) -> 0,
```

any fixed `a>0` should eventually be available after making the saddle
estimates effective. The tradeoff is simple:

```text
smaller a  => smaller N_geom(a), but usually larger T.
larger a   => easier analytic constant extraction, but larger N_geom(a).
```

## What Is Still Missing

H920 does not yet provide numerical values for `a` and `T`. The current
H910-H918 reports prove the needed asymptotic shape but do not expose explicit
constants.

So the remaining effective task is:

```text
H921 Explicit Saddle Constants.
Make H910-H918 effective enough to prove kappa'''(t) >= -a/t^2
for a concrete a<1 and all t>=T.
```

Once this is done:

```text
1. compute N_eff(a,T);
2. if N_eff(a,T) <= 127, cite H803/H804 finite data;
3. otherwise extend the finite H803 audit through N_eff(a,T)-1.
```

## No-Cheat Status

Closed here:

```text
non-effective eventual H908-A -> explicit finite closure criterion.
```

Still open:

```text
explicit constants a,T for the Xi saddle cumulant;
finite verification if N_eff(a,T)>127;
global all-degree Jensen/TP-to-RH assembly.
```

Thus H920 turns the H919 gap into a precise constant-extraction problem. It
does not close RH and does not claim total positivity of the de Bruijn-Newman
kernel.
