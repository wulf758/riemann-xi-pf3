# H1310 Mellin Reformulation and Finite Cumulant Hierarchy

Classification: `h1310_mellin_cumulant_hierarchy_open`

This note gives an exact independent reformulation of H1296.  It also records
which generic exponential-family properties are insufficient.  The resulting
cumulant inequalities are numerically supported, not proved.

## 1. Exact Mellin family

Put

```text
q(r)=r(pi exp(r)-alpha),
A=r q'(r),
B=r A'(r),
alpha=9/4.
```

With

```text
x=r exp(w/sqrt(A)),
```

the gamma-log measure becomes

```text
exp(-Psi_r(w))dw
 = C_r x^(q-1) exp(-pi exp(x)+alpha x)dx.
```

Thus, for

```text
Z(q)=int_0^infinity x^(q-1) exp(-pi exp(x)+alpha x)dx,
K(q)=log Z(q),
T=log X,
```

one has

```text
W=sqrt(A)(T-log r).
```

The natural potential of `T` is

```text
U(t)=pi exp(exp(t))-alpha exp(t),
```

and its mode satisfies `U'(log r)=q`.

## 2. Exact derivative identity

Let `M_n=E[W^n]`.  Differentiating both the exponential-family law and the
moving affine normalization gives

```text
r dM_n/dr
 = (nB/(2A))M_n
   +sqrt(A)(M_(n+1)-M_n M_1-nM_(n-1)).              (1)
```

For `n=2k`, (1) is nonpositive exactly when

```text
2k M_(2k-1)-M_(2k+1)+M_1 M_(2k)-k a M_(2k) >= 0,
```

which is H1296 after the Stein identity.  The Mellin change of variables is
therefore exact, but it does not make the original sign automatic.

## 3. A finite cumulant sufficient condition

Define signed, saddle-normalized cumulants by

```text
g_1 = sqrt(A)(log r-K'(q)),
g_j = (-1)^j A^(j/2) K^(j)(q),       j>=2.
```

If `g_1,...,g_6` are positive and nonincreasing in `r`, then
`M_2,M_4,M_6` are nonincreasing.  Indeed, the even Bell-polynomial expansions
contain only positive monomials in the `g_j`.  For example,

```text
M_2 = g_1^2+g_2,

M_4 = g_1^4+6g_1^2g_2+3g_2^2+4g_1g_3+g_4.
```

The same sign mechanism applies to every monomial of `M_6` because its total
cumulant weight is even.

Direct differentiation gives

```text
r g_1' = sqrt(A)(1-g_2+(a/2)g_1),

r g_j' = sqrt(A)((j a/2)g_j-g_(j+1)),    j>=2.      (2)
```

Consequently the following finite hierarchy is sufficient:

```text
g_2 >= 1+(a/2)g_1,
g_(j+1) >= (j a/2)g_j,       2<=j<=6.               (3)
```

A 147-point high-precision scan on `5/2<=r<=49/5` found
`g_1,...,g_7>0` and every margin in (3) positive.  The smallest margins occur
near `r=49/5` and range from approximately `4.9e-6` down to `8e-11`.
This is a numerical scout only.

## 4. Why generic Mellin structure is insufficient

The kernel `exp(qt)` is totally positive and the exact `X` density is strictly
log-concave.  These facts imply monotone likelihood ratio in `q`, but they do
not control moments after recentering at the moving mode and rescaling by its
curvature.

A concrete stress model is

```text
U(t)=100000 exp(t)+exp(200t).
```

At `t=0` it has

```text
q=100200,
A=140000,
a approximately 0.15463.
```

Its Mellin density is strictly log-concave, its family is totally positive,
and every derivative of `U` is positive.  Nevertheless high-precision
quadrature gives positive derivatives of all three standardized moments:

```text
dM_2/dt approximately 3.0005,
dM_4/dt approximately 34.6529,
dM_6/dt approximately 419.567.
```

Therefore total positivity, log-concavity, positivity of all derivatives,
large `q`, and a small cubic parameter do not prove H1296.  A successful proof
must use the exact factorial/Poisson coefficient profile, or prove the special
cumulant hierarchy (3) directly.

## Status

Proved here: the Mellin representation, derivative identity (1), cumulant
derivatives (2), and sufficiency of (3).

Open: positivity and monotonicity of the exact cumulants, hence H1296, H946,
the Xi-kernel transfer, and the Riemann Hypothesis.
