# H1316 Explicit Xi `m>=2` Tail Transfer

Classification: `h1316_xi_mtail_explicit_derivatives_proved`

## Result

Let

```text
I_(m,beta)(t)
 = integral_0^infinity U^(2t)
   exp(beta*U-pi*m^2*exp(4U)) dU
```

and put

```text
c_5 = 3/(2*pi),
R_tail(t)
 = sum_(m>=2) [
     m^4 I_(m,9)(t)/I_(1,9)(t)
     - c_5 m^2 I_(m,5)(t)/I_(1,9)(t)
   ].
```

For every `t>0` whose dominant gamma-log saddle parameter satisfies
`r=r_(2t)>=5/2`, the infinite series may be differentiated three times and

```text
0 <= R_tail(t) < 1/60,
|R_tail'(t)|   <=  1311/t,
|R_tail''(t)|  <=  6555/t^2,
|R_tail'''(t)| <= 48507/t^3.                    (1)
```

Thus the qualitative `O(exp(-c*t/U_t))` assertion in H918 is not needed for
the H923 transfer gate.  The whole arithmetic tail has explicit constants

```text
B_1,tail = 1311,
B_2,tail = 6555,
B_3,tail = 48507.
```

The constants are deliberately coarse.  Their role is finiteness and full
uniformity, not numerical optimization.

## 1. Sum Over `m` Before Integrating

Set `q=2t`, `rho=4U`, and `Z=log(rho)`.  Under the probability law

```text
dP_q(Z) proportional to exp(q*Z-V_alpha(Z)) dZ,
V_alpha(Z)=pi*exp(exp(Z))-alpha*exp(Z)-Z,
alpha=9/4,
```

the exact tail ratio is

```text
R_tail(t)=E_q[h(rho)],

h(rho)=sum_(m>=2)
  (m^4-c_5*m^2*exp(-rho))
  exp(-pi*(m^2-1)*exp(rho)).                    (2)
```

Every summand in (2) is positive because `m>=2`, `c_5<1/2`, and
`exp(-rho)<=1`.  For absolute-value estimates define the plus-sign majorant

```text
h_bar(rho)=sum_(m>=2)
  (m^4+c_5*m^2*exp(-rho))
  exp(-pi*(m^2-1)*exp(rho)).
```

Write `X=pi*exp(rho)>=pi`, `m=n+2`, and `y=exp(-4X)`.  Then

```text
m^2-4=n*(n+4)>=4n,
(n+2)^4<=16^(n+1),
(n+2)^2<=4^(n+1),
y<=exp(-4*pi)<1/1000.
```

Consequently

```text
h_bar(rho)
 <= exp(-3X) [16/(1-16y)+2/(1-4y)]
 < 19 exp(-3X).                                (3)
```

This is the promised pointwise summation.  In particular it is uniform in
both `m` and the saddle parameter.

## 2. Use The Global Sixth-Moment Bound Once

Let `z_q=log(r)` be the saddle, and use the exact normalizations

```text
A=V_alpha''(z_q),
W=sqrt(A)*(Z-z_q),
S=Z-z_q=W/sqrt(A).
```

The saddle identity and curvature identity are

```text
q=pi*r*exp(r)-alpha*r-1,
A=q*(r+1)+alpha*r^2+r+1 >= q*r.                (4)
```

H1313 supplies, uniformly for `r>=5/2`,

```text
E[W^2]<=26/25,
E[W^4]<=7/2,
E[W^6]<=21.                                    (5)
```

It follows from Cauchy--Schwarz that, for `0<=j<=3`,

```text
mu_j=E[|W|^j] <= 2,                            (6)
```

where `mu_0=1`; for `j=3`, use
`mu_3^2<=E[W^2]E[W^4]<4`.

Put

```text
d=sqrt(A)/r,
G={W>=-d}.
```

On `G`,

```text
rho=exp(Z)>=r*exp(-1/r)>=r-1,
pi*exp(rho)>=pi*exp(r-1)>=q/(e*r).
```

Therefore (3) gives

```text
h_bar <= 19 exp(-3q/(e*r)) on G.               (7)
```

On the complement, `|W|>=d`, so (5) gives, for `j<=3`,

```text
E[|W|^j 1_(G^c)] <= E[W^6]/d^(6-j)
                  <= 21/d^(6-j).              (8)
```

For

```text
n_j=E_q[h(rho) S^j],
D_j=E_q[S^j],
x=q/r,
```

(3), (7), and (8) imply

```text
|n_j|
 <= 19*mu_j*A^(-j/2)*exp(-3x/e)
    +399*r^(6-j)/A^3

 <= q^(-j) [38*x^(j/2)*exp(-3x/e)
             +399*x^(j-3)],                    (9)

|D_j| <= 2*q^(-j)*x^(j/2),   1<=j<=3.          (10)
```

Here (4) was used in the second lines.  Also `x>30`: indeed, for `r>=5/2`,

```text
x=pi*exp(r)-9/4-1/r
 >3*12-9/4-2/5>30.
```

The elementary inequality `exp(5/2)>12` follows already from the Taylor
polynomial through degree six.

For `j=0`, retaining `mu_0=1` in the first line of (9) yields

```text
0<=R_tail(t)
 <=19 exp(-3x/e)+399/x^3
 <19 exp(-30)+399/30^3
 <1/60.                                         (11)
```

## 3. Quotient Derivatives Without Hidden Saddle Derivatives

At the fixed value of `q` under consideration, freeze `z_q` as a constant and
multiply numerator and denominator of (2) by `exp(-q*z_q)`.  This is only a
common exponential factor, so differentiating the quotient gives exactly

```text
R_q' = n_1-n_0 D_1,

R_q''=n_2-2n_1D_1+n_0(2D_1^2-D_2),

R_q'''=n_3-3n_2D_1
       +3n_1(2D_1^2-D_2)
       +n_0(-6D_1^3+6D_1D_2-D_3).              (12)
```

No derivative of the saddle `z_q` occurs: the centering is frozen separately
at each evaluation point before applying the ordinary quotient rule.

For `0<=k<=3`,

```text
x^(k/2) exp(-3x/e)<=1.                          (13)
```

For `k>0`, maximize the left side on `(0,infinity)`; its maximum is
`(k/6)^(k/2)<=1`.  Moreover, whenever
`j+ell_1+...+ell_s=k<=3`, the algebraic exponent from the second term of
(9) is

```text
j-3+(k-j)/2=(j+k)/2-3<=0.
```

Hence each monomial

```text
n_j D_(ell_1)...D_(ell_s)
```

in (12) is at most

```text
437*2^s/q^k,     437=38+399,                    (14)
```

in absolute value.  Counting the coefficients in (12) gives

```text
|R_q'|   <= (1+2)*437/q             = 1311/q,
|R_q''|  <= (1+4+8+2)*437/q^2       = 6555/q^2,
|R_q'''| <= (1+6+24+6+48+24+2)
             *437/q^3               =48507/q^3. (15)
```

Since `q=2t`, the factors `2^j` from changing derivatives exactly cancel the
factors `(2t)^j` in the denominator.  Equations (11) and (15) prove (1).

## 4. Differentiation And Infinite-Sum Audit

The passage from the infinite series to (2), and differentiation up to order
three, are justified by dominated convergence.  The pointwise bound (3)
dominates the series uniformly.  After differentiation the only additional
factor is a power `|Z-z_q|^j`, `j<=3`; it is integrable by (5).  More directly,
at the two ends of the real `Z` axis the gamma-log density has respectively
exponential and double-exponential decay.

Potential failure modes checked explicitly:

1. The `beta=5` term includes `exp(-rho)`; it was retained in (2) and only
   bounded by one in the majorant.
2. The tail is a quotient by `I_(1,9)`, so raw numerator derivatives are not
   sufficient.  The exact quotient identities (12) include every denominator
   term.
3. Centering at the moving saddle does not introduce missing `z_q'` terms,
   because the common exponential is frozen at the evaluation point.
4. The sixth moment, not an assumed sub-Gaussian tail, controls the bad event.
5. The result is uniform for the entire infinite `m`-sum; no fixed-`m`
   asymptotic is used.

## 5. Scope

This closes only the `m>=2` part of the H923 correction.  To obtain explicit
constants for the full

```text
epsilon(t)=-c_5 I_(1,5)(t)/I_(1,9)(t)+R_tail(t),
```

one must combine (1) with a separate explicit estimate for the `m=1`,
`beta=5` ratio.  H1316 by itself does not prove the full Xi transfer, the
Jensen/total-positivity assembly, or the Riemann Hypothesis.

## Verifier

The exact coefficient arithmetic and elementary scalar inequalities are
checked by

```text
python tools/rh_h1316_xi_mtail_explicit_transfer.py
```

which reports `all_checks_pass=true`.

