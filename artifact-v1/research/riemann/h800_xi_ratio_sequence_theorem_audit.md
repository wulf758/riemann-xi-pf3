# H800 Xi Ratio-Sequence Theorem Audit

Classification: `ratio_logconvexity_bottleneck_identified`

## Question

H799 reduces the H798 two-near one-far sparse order-3 family to a ratio-tail criterion. For `a_n=gamma_n/n!`, set `R_n=a_n/a_{n-1}`, `q_n=R_n/R_{n-1}`, `A=R1`, and `B=R2`. The H799 criterion needs:

- positivity of `a_n`;
- `R_n` nonincreasing, equivalently adjacent log-concavity of `a_n`;
- `q_{n+1}>=q_n`, equivalently log-convexity of the ratio sequence `R_n`;
- tail smallness `B>q_m*R_{m-1}*(1+q_m)`, with `m=3` allowed as a finite base case.

## Primary Sources Checked

1. Griffin, Ono, Rolen, Thorner, Tripp, Wagner, `Jensen Polynomials for the Riemann Xi Function`, arXiv:1910.01227, https://arxiv.org/abs/1910.01227

   Relevant facts: the paper defines `psi(z)=sum gamma(j)/j! z^(2j)=xi(1/2+z)`, states positivity of `gamma(n)`, recalls Polya's equivalence between RH and hyperbolicity of all Jensen polynomials, and records that hyperbolicity is known for all `n>=0` and `d<=3`. It also gives effective eventual hyperbolicity for fixed degree.

2. Griffin, Ono, Rolen, Zagier, `Jensen polynomials for the Riemann zeta function and other sequences`, arXiv:1902.07321, https://arxiv.org/abs/1902.07321

   Relevant facts: gives asymptotics and Hermite-polynomial modeling for Jensen polynomials. This supports eventual fixed-degree Turan-type information, not a direct all-index proof of H799's quotient monotonicity.

3. O'Sullivan, `Zeros of Jensen polynomials and asymptotics for the Riemann xi function`, arXiv:2007.13582, https://arxiv.org/abs/2007.13582

   Relevant facts: gives detailed asymptotic expansions for Taylor coefficients of the xi function and related quantities. This is relevant for eventual tail estimates, but it does not directly state the H799 ratio log-convexity condition.

4. Farmer, `Jensen polynomials are not a plausible route to proving the Riemann Hypothesis`, arXiv:2008.07206, https://arxiv.org/abs/2008.07206

   Relevant facts: warns that Jensen-polynomial asymptotics and eventual hyperbolicity should not be treated as evidence strong enough to attack RH directly. This supports keeping H800 as a local lemma audit, not an RH claim.

## What Seems Covered

### Positivity

The xi Jensen paper states positivity of the `gamma(n)` coefficients. Since `a_n=gamma_n/n!`, positivity transfers immediately to `a_n`.

### Adjacent Log-Concavity Of `a_n`

The xi Jensen paper records global hyperbolicity for `d<=3`, and in particular uses the `d=2` Turan inequality as log-concavity of `gamma(n)`. If

```text
gamma_n^2 >= gamma_{n-1}*gamma_{n+1},
```

then

```text
a_n^2/(a_{n-1}*a_{n+1})
= (gamma_n^2/(gamma_{n-1}*gamma_{n+1})) * ((n+1)/n) >= 1.
```

So `R_n=a_n/a_{n-1}` is globally nonincreasing, assuming the cited global `d=2` result.

### Eventual Tail Smallness

The generating function `sum gamma_n z^n/n!` is entire in the relevant variable, so `limsup a_n^(1/n)=0`. If `R_n` is positive and nonincreasing, then `R_n` has a limit; a positive limit would force finite radius of convergence. Hence `R_n -> 0`, so the coarse H799 tail bound `R_{m-1}<R2/2` is eventually true. H799 already verifies the finite window through `m=22`, but this audit has not produced an explicit analytic cutoff.

## Bottleneck Not Covered

H799 also needs

```text
q_{n+1}=R_{n+1}/R_n >= q_n=R_n/R_{n-1}.
```

Equivalently,

```text
R_n^2 <= R_{n-1}*R_{n+1}.
```

In terms of `a_n`, this is

```text
a_n^3*a_{n-2} <= a_{n-1}^3*a_{n+1}.
```

In terms of `gamma_n`, the factorial normalization changes it to

```text
(gamma_n^3*gamma_{n-2})/(gamma_{n-1}^3*gamma_{n+1}) <= n^2/(n^2-1).
```

I did not find this exact ratio-log-convexity theorem in the primary sources checked. Known Jensen/Turan hyperbolicity results give important adjacent and finite-degree constraints, and asymptotic papers can plausibly help eventually, but none of the checked source statements directly closes this H799 condition globally without additional work.

## Route Judgment

H799 survives the audit, but H800 identifies the real bottleneck: prove or invalidate ratio log-convexity of `R_n=a_n/a_{n-1}` for the xi coefficients. This is a cleaner global lemma target than generic sparse-minor mining. It is also not already solved by the standard Jensen asymptotic literature in the form needed here.

## Not Proved

- Global `q_{n+1}>=q_n` for the normalized xi coefficient ratios `R_n=a_n/a_{n-1}`.
- An explicit analytic cutoff for the H799 tail bound.
- The H798 family for all `m`.
- PF3, JP2', PF-infinity, or RH.

## Next Target

`H801 ratio-log-convexity attack`: derive `q_{n+1}>=q_n` from an integral representation, a coefficient asymptotic with explicit remainder, or invalidate it with a higher-precision coefficient audit beyond H784.
