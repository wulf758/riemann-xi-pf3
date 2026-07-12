# H1317 Beta-5 H923 Derivative Certificate

Classification: `h1317_beta5_h923_derivatives_proved`

## Statement

Let

```text
I_beta(t)=integral_0^infinity
  U^(2t) exp(beta*U-pi*exp(4U)) dU,
R_5(t)=I_5(t)/I_9(t).
```

For every dominant gamma-log saddle with `r>=5/2`,

```text
|R_5'(t)|   <= 207/t,
|R_5''(t)|  <= 1035/t^2,
|R_5'''(t)| <= 7659/t^3.                            (1)
```

These are deliberately coarse absolute bounds.  Together with the H1316
`m>=2` tail certificate, they supply every derivative hypothesis in H923.

## Dominant Law and Frozen-Saddle Quotients

Put `q=2t`, `rho=4U`, `Z=log(rho)`, and `alpha=9/4`.  The dominant law is

```text
dP_q(Z) proportional to exp(qZ-V_alpha(Z)) dZ,
V_alpha(Z)=pi*exp(exp(Z))-alpha*exp(Z)-Z.
```

Let `z_q=log(r)` be its saddle and define

```text
A=V_alpha''(z_q),
W=sqrt(A)(Z-z_q),
S=Z-z_q=W/sqrt(A),
x=q/r.
```

The exact saddle and curvature identities give

```text
q=pi*r*exp(r)-alpha*r-1,
A=q(r+1)+alpha*r^2+r+1 >= q*r,
x=pi*exp(r)-alpha-1/r > 30                         (2)
```

for `r>=5/2`.

In these variables

```text
R_5(q)=E_q[f(Z)],  f(Z)=exp(-exp(Z)).               (3)
```

At a fixed value of `q`, freeze `z_q` and multiply numerator and denominator
of (3) by `exp(-q*z_q)`.  Define

```text
f_0=f(z_q)=exp(-r),
D_j=E[S^j],
n_j=E[(f-f_0)S^j].                                 (4)
```

The constant component `f_0` cancels identically from the quotient
derivatives.  Ordinary quotient differentiation therefore gives, exactly,

```text
R_q' = n_1-n_0 D_1,

R_q''=n_2-2n_1D_1+n_0(2D_1^2-D_2),

R_q'''=n_3-3n_2D_1
       +3n_1(2D_1^2-D_2)
       +n_0(-6D_1^3+6D_1D_2-D_3).                 (5)
```

There is no missing derivative of `z_q`: the centering is frozen at the
evaluation point before the quotient rule is applied.

## Centered Numerator Bound

H1313 proves

```text
E[W^2]<=26/25, E[W^4]<=7/2, E[W^6]<=21.            (6)
```

In particular `E|W|^j<=4` for `1<=j<=4`.

Put

```text
d=sqrt(A)/r,
G={W>=-d}={S>=-1/r}.
```

On the line segment from `z_q` to any `Z` in `G`,

```text
rho=exp(Z)>=r*exp(-1/r)>=r-1>=3/2.
```

Since `rho*exp(-rho)` decreases for `rho>=1`, the mean-value theorem gives

```text
|f(Z)-f_0| <= (r-1)exp(1-r)|S|                     (7)
```

on `G`.  Hence its contribution to `|n_j|` is at most

```text
4(r-1)exp(1-r) A^(-(j+1)/2).                       (8)
```

On `G^c`, use `|f-f_0|<=1`, `|W|>=d`, and (6):

```text
E[|f-f_0||S|^j 1_(G^c)]
 <=A^(-j/2) E[|W|^j 1_(|W|>=d)]
 <=21 A^(-j/2)d^(-(6-j))
 =21r^(6-j)/A^3.                                   (9)
```

Combining (8), (9), and `A>=qr` gives, for `0<=j<=3`,

```text
|n_j| <= q^(-j)[
  4e*x^((j-1)/2)*exp(-r)+21*x^(j-3)].              (10)
```

The denominator moments satisfy, for `1<=j<=3`,

```text
|D_j|<=2A^(-j/2)<=2q^(-j)x^(j/2).                  (11)
```

## Product Bound and Constants

Consider a monomial

```text
n_j D_(ell_1)...D_(ell_s)
```

of total derivative order

```text
k=j+ell_1+...+ell_s<=3.
```

The local power in (10)-(11) is

```text
x^((k-1)/2)exp(-r)<=4.                              (12)
```

Indeed `x<pi*exp(r)`; the worst case is `k=3`, where the left side is
`x*exp(-r)<pi<4`.  The other cases are smaller because `r>=5/2` and `x>30`.

The bad-event power is

```text
x^(j-3+(k-j)/2)=x^((j+k)/2-3)<=1.                  (13)
```

Using `e<3`, every monomial is therefore bounded by

```text
69*2^s/q^k,  69=48+21.                             (14)
```

Counting the coefficients in (5) gives respectively

```text
3, 15, 111.
```

Thus

```text
|R_q'|   <= 3*69/q       =207/q,
|R_q''|  <=15*69/q^2     =1035/q^2,
|R_q'''| <=111*69/q^3    =7659/q^3.                (15)
```

Finally `q=2t`.  The factor `2^j` from `d^j/dt^j` cancels the denominator
`q^j=(2t)^j`, proving (1).

## H923 Assembly

Write

```text
epsilon(t)=-c_5 R_5(t)+R_tail(t),
c_5=3/(2*pi)<1/2.
```

The independently audited H1316 `m`-tail theorem gives

```text
|R_tail'|<=1311/t,
|R_tail''|<=6555/t^2,
|R_tail'''|<=48507/t^3.                            (16)
```

Equations (1) and (16) give the valid integer constants

```text
B_1=1415,
B_2=7073,
B_3=52337.                                         (17)
```

The exact pointwise correction audit gives

```text
-c_5<epsilon(t)<0,
```

so `|epsilon|<1/2`.  H923 now yields

```text
d^3/dt^3 log(1+epsilon(t))
 >= -B_Xi/t^3,

B_Xi=2B_3+12B_1B_2+16B_1^3
    =45,450,578,214.                               (18)
```

In particular,

```text
d^3/dt^3 log(1+epsilon(t)) >= -1/(4t^2)            (19)
```

whenever

```text
t>=4B_Xi=181,802,312,856.                          (20)
```

## Scope

This closes the full explicit H923 transfer gate.  The constants are not
sharp.  More importantly, the gamma-log input H922m starts only at
`q=2t>=exp(100)`, so the resulting Xi theorem still leaves an enormous finite
interval and does not prove H803, the global Jensen/TP assembly, or RH.

