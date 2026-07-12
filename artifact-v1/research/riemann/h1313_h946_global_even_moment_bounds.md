# H1313 Global H946 Even-Moment Bounds

Classification: `h1313_h946_global_even_moment_subproblem_proved`

## Result

For the exact gamma-log normalized saddle with `alpha=9/4`, the three coarse
H946 moment targets hold for every `r>=5/2`:

```text
E_r[W^2] <= 26/25,
E_r[W^4] <= 7/2,
E_r[W^6] <= 21.                              (1)
```

This closes the even-moment subproblem used in H946.  It does **not** prove
the H946 corrected-pair pointwise inequality, H943-C, any Xi-kernel transfer,
or the Riemann Hypothesis.

The dependency and endpoint verifier is

```text
tools/rh_h1313_h946_global_moment_assembly.py
```

and reports `all_checks_pass=true`.

## Canonical Anchor At `s=30`

For the canonical log-Gamma model, let

```text
Y_s ~ Gamma(s,1),
W_s = sqrt(s)*(log(Y_s)-log(s)).
```

Its cumulants are

```text
kappa_1 = sqrt(s)*(digamma(s)-log(s)),
kappa_j = s^(j/2)*(-1)^j*(j-1)!*zeta(j,s),   j>=2.
```

The verifier evaluates these with Arb at the exact integer `s=30` and uses
the complete Bell recurrence

```text
m_n = sum_{j=1}^n binom(n-1,j-1)*kappa_j*m_(n-j),
m_0 = 1.
```

It certifies

```text
M_2(30) = 1.02527798352001903910376321111276544...
          < 26/25,

M_4(30) = 3.29282565411347218644581918283424841...
          < 7/2,

M_6(30) = 18.7088437782569320159946190958899241...
          < 21.
```

The corresponding certified slacks are approximately

```text
0.01472201648,
0.20717434589,
2.29115622174.
```

H1302 proves that every canonical even moment `M_(2k)(s)` is strictly
decreasing in `s`.  Consequently these three inequalities remain true for
every canonical shape `s>=30`.

## Band 1: `5/2 <= r <= 3`

H1311 proves the H1296 covariance inequalities for `k=1,2,3` throughout this
band.  Hence the exact moments `M_2(r),M_4(r),M_6(r)` are nonincreasing in
`r` there.

H1282 supplies the rigorous endpoint anchor at `r=5/2`:

```text
M_2(5/2) <= 1.039525544768851... < 26/25,
M_4(5/2) <= 3.461554174436287... < 7/2,
M_6(5/2) <= 20.990414697940364... < 21.
```

Therefore (1) holds on `5/2<=r<=3`.

Dependencies:

```text
research/riemann/h1282_gamma_log_endpoint_even_moment_certificate_360k.json
research/riemann/h1311_h1296_finite_covariance_certificate_to_3.md
tools/rh_h1306b_finite_covariance_certificate_pilot.py
tools/rh_h1306c_paired_finite_covariance_pilot.py
tools/rh_h1306d_paired_cover_pilot.py
```

## Band 2: `3 <= r <= 49/5`

H1312 proves that the normalized exact radial law of `|W|` is
likelihood-ratio dominated by the canonical radial law with the same local
skew parameter.  Thus every increasing function of `|W|`, in particular
`|W|^2`, `|W|^4`, and `|W|^6`, has no larger expectation under the exact law.

The same exact H1312 certificate proves

```text
s(r)=1/a(r)^2 >= 30
```

on this band.  H1302 and the Arb anchor at `s=30` therefore give

```text
E_exact[W^(2k)]
 <= E_canonical,s(r)[W^(2k)]
 <= E_canonical,30[W^(2k)]
```

for `k=1,2,3`.  The last quantities satisfy (1), so the exact moments do as
well throughout `3<=r<=49/5`.

Dependencies:

```text
research/riemann/h1302_canonical_log_gamma_even_moment_monotonicity.md
research/riemann/h1312_radial_likelihood_ratio_certificate.md
tools/rh_h1312_radial_lr_moment_certificate.py
```

## Band 3: `r >= 49/5`

H1279 independently proves the three bounds (1) on the whole high-`r` range,
using a certified curvature sandwich on `49/5<=r<=100`, H1274 beyond `100`,
and rigorous left and right tails.

Dependency:

```text
research/riemann/h1279_gamma_log_threshold49over5_moments.json
```

## Coverage

The three closed bands are

```text
[5/2,3],
[3,49/5],
[49/5,infinity).
```

Their endpoints agree exactly, so there is no uncovered interval.  Combining
the three band arguments proves (1) for every `r>=5/2`.

## Exact Scope

Proved by H1313:

```text
the three global coarse even-moment inequalities (1).
```

Not proved by H1313:

```text
the corrected-pair inequality R_r(t)<=1/2+C_r*t^2,
the coefficient bound needed to finish that separate H946 route,
H943-C or any later kernel budget,
the gamma-log-to-Xi transfer,
the Riemann Hypothesis.
```
