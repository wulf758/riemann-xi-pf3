# H1312 Radial Likelihood-Ratio Certificate On `3 <= r <= 49/5`

Classification: `h1312_exact_radial_lr_domination`

## Result

Let

```text
g_r(t)     = exp(-Psi_r(t)) + exp(-Psi_r(-t)),
g_can,a(t) = exp(-Phi_a(t)) + exp(-Phi_a(-t)),
Phi_a(t)   = (exp(a*t)-1-a*t)/a^2.
```

For the exact gamma-log potential with `alpha=9/4`, and for every

```text
3 <= r <= 49/5,  t >= 0,
```

the radial density ratio is nonincreasing:

```text
d/dt log(g_r(t)/g_can,a(t)) <= 0.
```

Consequently the normalized law of `|W|` under the exact potential is
likelihood-ratio dominated by the canonical log-Gamma radial law.  In
particular, every increasing function of `|W|` has no larger expectation
under the exact law, whenever the expectations exist.

The exact verifier is
`tools/rh_h1312_radial_lr_moment_certificate.py`.

## Dimensionless Reduction

Use the positive exponential decomposition from H1303.  If `M` has law

```text
nu_m = c_m*m^2/A,
```

put

```text
mu = E[M] = B/A,
X  = M/mu,
s  = 1/a^2 = A^3/B^2,
z  = a*t.
```

Then `E[X]=1` and

```text
Psi_r(t) = s E[(exp(zX)-1-zX)/X^2],
Phi_a(t) = s (exp(z)-1-z).
```

Let `H=Psi-Phi` and define

```text
Delta = Phi(t)-Phi(-t) = 2s(sinh(z)-z),
D     = H(t)-H(-t),
p     = 1/(1+exp(Delta)),
p_D   = 1/(1+exp(Delta+D)),
L     = -H'(-t),
ell   = a*L
      = E[(1-exp(-zX))/X] - (1-exp(-z)).
```

The radial log-ratio derivative is exactly

```text
d/dt log(g_r/g_can)
  = -L - p_D*D' + (p-p_D)*Delta'.
```

Jensen gives `L>=0`, `D>=0`, and `D'>=0`.  It is therefore sufficient to
prove

```text
ell >= 2(p-p_D)(cosh(z)-1).                 (1)
```

The two elementary logistic bounds used below are

```text
p-p_D <= D/4,
p-p_D <= exp(-Delta)(1-exp(-D)).            (2)
```

## Exact Parameter Lemmas

Write `P=pi*exp(r)` and let `T_j` be the Touchard polynomial.  Set

```text
V_j = P*T_j(r) - alpha*r,
A   = V_2,
B   = V_3,
G   = A*V_4-B^2.
```

For `delta_j=E[X^j]-1`, the verifier proves exactly on
`3<=r<=49/5` that

```text
s >= 30,
delta_2 >= 1/20,
delta_3 <= (7/2) delta_2,
delta_4 >= 5 delta_2,
E[X^5] <= 30 delta_2.                       (3)
```

The first, third, and fourth inequalities reduce to numerators whose
coefficients are all strictly positive after

```text
X0 = r-3,
Y0 = P-60.
```

The coefficient counts and smallest coefficients are respectively

```text
s>=30:                    28 coefficients, minimum 1,
delta_3/delta_2<=7/2:     36 coefficients, minimum 1/2,
delta_4/delta_2>=5:       60 coefficients, minimum 1.
```

For `delta_2>=1/20` and `E[X^5]<=30 delta_2`, each coefficient of
`Y0^j` is a polynomial in `r`.  Its exact Bernstein coefficients on
`[3,49/5]` are all positive.  The verifier checks 3 rows of Bernstein
coefficients for the first inequality and 6 rows for the second.

Finally `P>60` follows without a floating-point constant: `pi>3` and the
degree-eight Taylor floor gives

```text
exp(3) > 89641/4480,
3*(89641/4480)-60 = 123/4480 > 0.
```

## Tilted Fifth-Moment Bound

For `0<=z<=1`, put `q=r*exp(z/mu)`.  The exact identity

```text
B-(r+3/2)A
  = P*r*(r-1)/2 + alpha*r*(r+1/2) > 0
```

gives `mu>=r+3/2`.  The elementary inequality
`log(1+1/r)>=1/(r+1)` then gives

```text
r <= q <= r+1 <= 54/5.
```

If `V_j(q)=pi*exp(q)*T_j(q)-alpha*q`, then

```text
d/dz log E[X^5 exp(zX)] = V_8(q)/(mu*V_7(q)).
```

The Touchard identity

```text
((q+6)T_7(q)-T_8(q))/q
 = (140-q^2)q^4 + 700q^3 + 903q^2 + 252q + 5
```

is positive because `(54/5)^2<140`.  The `-alpha*q` correction is harmless:
using `pi*exp(q)>60`, the displayed polynomial is at least `5`, while

```text
5*60-alpha*(54/5+5) = 5289/20 > 0.
```

Thus `V_8/V_7<=q+6`.  Since

```text
q+6 <= r+7 <= (5/2)mu,
```

integration and (3) give

```text
E[X^5 exp(zX)] <= 30 delta_2 exp(5z/2).      (4)
```

## Compact Part `0 <= z <= 1`

Taylor inequalities with their correct global remainder signs yield

```text
ell >= L5
    = delta_2*z^3/6
      -delta_3*z^4/24
      +delta_4*z^5/120
      -E[X^5]*z^6/720.
```

Using (3),

```text
L5 >= delta_2*z^3*Q(z),
Q(z) = 1/6 - 7z/48 + z^2/24 - z^3/24 > 0.
```

Similarly,

```text
D <= D_up
  = 2s[delta_3*z^5/120
       +z^7*E[X^5 exp(zX)]/5040]
  <= 2s*delta_2*z^5*K(z),

K(z) = 7/240 + z^2*exp(5z/2)/168.
```

Here the only infinite Taylor tail was absorbed by

```text
sinh(y)-y-y^3/6-y^5/120 <= y^7 exp(y)/7!.
```

Combining (2), `cosh(z)-1<=(cosh(1)-1)z^2`, and
`Delta>=s*z^3/3`, the right side of (1), divided by `L5`, is at most

```text
R(z) = 4(cosh(1)-1) * z*K(z)/Q(z) * m_*(z),
```

where

```text
m_*(z) = sup_{y>=30z^3} y*min(1/4,exp(-y/3)).
```

The maximum in `y` occurs at `3 log 4`.  The verifier uses the rational
split `z=13/25`:

```text
m_*(z) <= 3log(4)/4                         for z<=13/25,
m_*(z) <= 30z^3 exp(-10z^3)                 for z>=13/25.
```

Arb checks 64 rational boxes in each branch.  The worst interval upper
bound returned by the verifier is

```text
R(z) < 0.486278 < 199/250 < 1.
```

This proves (1) on `0<=z<=1` with substantial room.

## Tail `z >= 1`

The function `ell` is increasing because

```text
ell'(z) = E[exp(-zX)]-exp(-z) >= 0.
```

At `z=1`, the compact lower bound and (3) give

```text
ell(1) >= delta_2*Q(1) >= 1/960,
Q(1)=1/48.
```

On the other hand, (2) and `s>=30` give

```text
2(p-p_D)(cosh(z)-1)
 <= 2(cosh(z)-1) exp(-60(sinh(z)-z)).
```

This last expression decreases for `z>=1`.  Arb certifies

```text
value at z=1 < 2.956e-5 < 1/960,
logarithmic derivative at z=1 < -30.42.
```

The derivative remains negative because `coth(z/2)` decreases and
`cosh(z)-1` increases.  This proves (1) on the whole tail.

## Scope

H1312 is a rigorous comparison theorem for the exact gamma-log potential on
`3<=r<=49/5`.  It is not a proof of RH.  The adjacent interval
`5/2<=r<=3` is handled independently by the H1311/H1306d covariance
certificate rather than by this likelihood-ratio argument.
