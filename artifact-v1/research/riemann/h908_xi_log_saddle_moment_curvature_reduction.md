# H908 Xi Log-Saddle Moment-Curvature Reduction

Classification: `h908_xi_log_saddle_third_difference_reduction_proved`

## Question

Can the H803/H804 Xi moment-curvature target be turned into an analytic
log-saddle condition, rather than another finite PF3 or coefficient audit?

## Short Answer

Yes. The H803 target is exactly a third finite-difference lower bound for the
log-Laplace transform of the positive Xi kernel. Equivalently, it is a lower
bound on an averaged tilted third cumulant of the log-height variable.

This is not a proof of RH. It is a cleaner analytic reduction for the selected
H907 total-positivity branch.

## Moment Setup

Use the standard even Xi integral representation in the form

```text
Xi(z)=xi(1/2+z)=int cosh(zu) dPhi(u),
```

with a positive even kernel measure after symmetrization. Equivalently, for the
even Taylor coefficients,

```text
Xi(z)=sum_{n>=0} a_n z^{2n},
M_n=(2n)!a_n=int u^{2n} dPhi(u).
```

Put

```text
x = u^2,
y = log x,
dnu(y) = pushforward of dPhi(u) under y=log(u^2),
M_t = int exp(t y) dnu(y),
kappa(t)=log M_t.
```

Then `M_n` is the moment sequence used in H803.

## H803 Variables

H803 defines

```text
S_n = M_n/M_{n-1},
T_n = S_n/S_{n-1} = M_n M_{n-2}/M_{n-1}^2.
```

H804 showed that the stronger eventual unit barrier

```text
T_{n+1}/T_n >= 1 - 1/n^2
```

implies the exact H803 factorial-threshold inequality, after finite base cases
`n=2,...,7` are checked separately.

## Exact Third-Difference Identity

Since `M_t=exp(kappa(t))`,

```text
T_n = exp(kappa(n)-2kappa(n-1)+kappa(n-2)).
```

Therefore

```text
T_{n+1}/T_n
 = exp(kappa(n+1)-3kappa(n)+3kappa(n-1)-kappa(n-2))
 = exp(Delta^3 kappa(n-2)),
```

where `Delta f(t)=f(t+1)-f(t)`.

Thus H804's unit barrier is equivalent to

```text
Delta^3 kappa(n-2) >= log(1 - 1/n^2).        (U_n)
```

The exact H803 threshold is similarly equivalent to

```text
Delta^3 kappa(n-2)
 >= log(C_{n-1}C_{n+1}/C_n^2),
C_n=(2n)(2n-1).
```

## Cumulant Interpretation

For real `t`, define the tilted probability measure

```text
dnu_t(y)=exp(t y-kappa(t)) dnu(y).
```

Then

```text
kappa'(t)  = E_t[y],
kappa''(t) = Var_t(y) >= 0,
kappa'''(t)=E_t[(y-E_t y)^3].
```

Iterating the identity

```text
Delta f(x)=int_0^1 f'(x+r) dr
```

gives

```text
Delta^3 kappa(x)
 = int_0^1 int_0^1 int_0^1
     kappa'''(x+r_1+r_2+r_3) dr_1 dr_2 dr_3.
```

Consequently, H804's eventual unit barrier follows if one proves, for all large
`n`,

```text
int_[0,1]^3 kappa'''(n-2+r_1+r_2+r_3) dr_1 dr_2 dr_3
 >= log(1 - 1/n^2).
```

A simpler sufficient condition is the pointwise skewness lower bound

```text
kappa'''(t) >= log(1 - 1/n^2)
for all t in [n-2,n+1].
```

Since

```text
log(1 - 1/n^2) = -1/n^2 + O(1/n^4),
```

the analytic target is a lower bound of order `-1/n^2` for the tilted third
central moment of the log-height variable near saddle index `n`.

## Why This Is Better Than Another Grid

H804's finite data suggested

```text
n^2(1 - T_{n+1}/T_n) ~ 0.5 to 0.6
```

in the tested range, while the unit barrier allows constant `1`. H908 explains
what must be proved analytically: the log-saddle skewness may be mildly
negative, but not more negative than order `1/n^2`.

This is a standard asymptotic problem for a positive-kernel moment sequence. It
is Xi-specific because H802 already showed that generic log-concavity and entire
moment growth do not imply the required ratio-log-convexity.

## Remaining Lemma H908-A

The new proof obligation is:

```text
H908-A Log-Saddle Skewness Bound.
For the Xi kernel moment measure dnu, prove that for all sufficiently large n,

  Delta^3 kappa(n-2) >= log(1 - 1/n^2),

or prove the sharper exact H803 threshold directly.
```

A plausible route is saddle-point analysis for

```text
M_t=int exp(t y)dnu(y),
```

including explicit control of the third derivative of `kappa(t)` with error
below the `1/n^2` scale.

## Finite Base

H804 records that the stronger unit barrier fails for `n=2,...,7` but holds in
the tested tail from `n=8` onward. Therefore H908-A only needs an eventual
proof plus separate finite verification of the H803 exact threshold for the
base range.

## Route Decision

H908 promotes the total-positivity branch from finite PF3 evidence to a precise
analytic moment problem:

```text
prove a log-saddle third-cumulant lower bound for the Xi kernel.
```

If this skewness bound cannot be proved or is false, the H907 selected TP branch
should be demoted rather than mined through more local minors.