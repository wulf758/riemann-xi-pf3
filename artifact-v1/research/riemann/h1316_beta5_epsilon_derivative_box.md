# H1316 Beta-5 Epsilon Derivative Box

Classification: `h1316_beta5_h923_derivative_box_proved`

## Result

Use the exact notation

```text
q=2t,
R_5(t)=I_5(t)/I_9(t),
c_5=3/(2pi),
epsilon_5(t)=-c_5 R_5(t).
```

For every `q>=exp(100)` (hence `t>=exp(100)/2`),

```text
|R_5'(t)|   <= 207/t,
|R_5''(t)|  <= 1035/t^2,
|R_5'''(t)| <= 7659/t^3.                            (1)
```

Consequently the beta-5 correction alone satisfies H923 with

```text
|epsilon_5|       < 1/2,
|epsilon_5'|      <= (207/2)/t,
|epsilon_5''|     <= (1035/2)/t^2,
|epsilon_5'''|    <= (7659/2)/t^3.                  (2)
```

Thus

```text
|d^3/dt^3 log(1+epsilon_5(t))|
 <= 18389880/t^3.                                  (3)
```

No asymptotic `O(.)` term occurs in (1)-(3).

## Dominant Probability Law

Let `P_q` be the `alpha=9/4` gamma-log law in `z`, and let

```text
z_q=log r,
S=z-z_q,
W=sqrt(A) S,
f(z)=exp(-exp(z)),
f_0=f(z_q)=exp(-r).
```

Then

```text
R_5=E_q[f].                                         (4)
```

Put

```text
D_j=E_q[S^j],
N_j=E_q[f S^j],
n_j_tilde=E_q[(f-f_0)S^j].                          (5)
```

Thus `N_j=f_0 D_j+n_j_tilde`.

## Exact Derivative Formulas And Cancellation

Differentiation in the natural parameter `q` gives

```text
partial_q R_5
 = n_1_tilde-n_0_tilde D_1,

partial_q^2 R_5
 = n_2_tilde-2n_1_tilde D_1
   +n_0_tilde(2D_1^2-D_2),

partial_q^3 R_5
 = n_3_tilde-3n_2_tilde D_1
   +3n_1_tilde(2D_1^2-D_2)
   +n_0_tilde(-6D_1^3+6D_1D_2-D_3).                (6)
```

These are the usual first three derivatives of a ratio of Laplace
transforms. Substituting `N_j=f_0D_j+n_j_tilde` makes every term containing
the constant `f_0` cancel identically. This cancellation is essential: a
bound on `f` without subtracting `f_0` loses the required power of `q`.

## Bounds For The Centered Numerators

For `j=0,1,2,3`, split at

```text
G={W>=-sqrt(A)/r}={S>=-1/r}.                         (7)
```

On the line segment from `z_q` to any `z` in `G`, writing `rho=exp(z)`,

```text
rho>=r exp(-1/r)>=r-1.
```

Since `r>=5/2`, the function `rho exp(-rho)` is decreasing throughout this
range. The mean-value theorem therefore gives

```text
|f-f_0| <= (r-1)exp(1-r)|S|                         (8)
```

on `G`.

H1313 gives `M_2<=26/25`, `M_4<=7/2`, and `M_6<=21`. Hence
`E|W|^ell<=4` for `1<=ell<=4`, and (8) yields

```text
E[|f-f_0||S|^j 1_G]
 <=4(r-1)exp(1-r) A^(-(j+1)/2).                     (9)
```

On `G^c`, `|f-f_0|<=1` and `|W|>=sqrt(A)/r`. Using the sixth moment,

```text
E[|f-f_0||S|^j 1_(G^c)]
 <=21 r^(6-j)/A^3.                                  (10)
```

Let `x=q/r`. Since `A>=qr`, (9)-(10) give

```text
|n_j_tilde|
 <=q^(-j)[4e x^((j-1)/2) exp(-r)+21x^(j-3)].        (11)
```

The same moment bounds give, for `ell=1,2,3`,

```text
|D_ell|<=2A^(-ell/2)
          <=2q^(-ell)x^(ell/2).                     (12)
```

## Constant Count

H922m gives `r>=50` and `r<=log q` for `q>=exp(100)`, so in particular
`x=q/r>30>1`. The saddle equation gives

```text
x exp(-r)=pi-(9/4+1/r)exp(-r)<pi<4.                 (13)
```

In a monomial of total derivative order `k<=3`, the good-set power in
(11)-(12) is

```text
x^((k-1)/2)exp(-r)<=x exp(-r)<4,
```

while the bad-set exponent is

```text
(j+k-6)/2<=0.
```

Using `e<3`, every monomial before its integer coefficient and factors `2`
from (12) is therefore bounded by

```text
(4e*4+21)q^(-k)<69q^(-k).                           (14)
```

The coefficient counts in (6), including one factor `2` per `D`, are

```text
k=1: 1+2=3,
k=2: 1+4+8+2=15,
k=3: 1+6+30+48+24+2=111.                            (15)
```

Equations (14)-(15) prove

```text
|partial_q R_5|   <=207/q,
|partial_q^2 R_5| <=1035/q^2,
|partial_q^3 R_5| <=7659/q^3.                       (16)
```

Since `q=2t`, multiplication by `2^j` converts (16) into (1) with the same
constants.

## H923 Log-Correction Constant

Because `0<R_5<=1` and `c_5<1/2`, (2) holds. Applying H923 with

```text
B_1=207/2, B_2=1035/2, B_3=7659/2
```

gives

```text
B_beta5=2B_3+12B_1B_2+16B_1^3
       =18389880,
```

which proves (3).

## Scope And Assembly Warning

This closes the H923 derivative box for the `m=1`, `beta=5` correction.
It does not by itself bound the `m>=2` tail. In the full correction,

```text
log(1-c_5R_5+R_tail)
```

must not be replaced by `log(1-c_5R_5)+log(1+R_tail)`. One valid assembly is
to factor

```text
1-c_5R_5+R_tail
 =(1-c_5R_5)[1+R_tail/(1-c_5R_5)]
```

and bound the derivatives of the second factor, including its interactions
with `(1-c_5R_5)^(-1)`.

The companion H1316 direct-cumulant proof independently gives
`|g_5'''(t)|<=27/(1250t^2)` and serves as a cross-check, but (3) is the bound
in the exact format requested by H923.
