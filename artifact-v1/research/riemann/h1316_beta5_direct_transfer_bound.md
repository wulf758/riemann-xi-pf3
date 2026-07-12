# H1316 Beta-5 Direct Xi Transfer Bound

Classification: `h1316_beta5_direct_log_correction_bound_proved`

## Result

Let

```text
I_beta(t) = int_0^infty U^(2t) exp(beta U-pi exp(4U)) dU,
R_5(t) = I_5(t)/I_9(t),
c_5 = 3/(2pi),
g_5(t) = log(1-c_5 R_5(t)).
```

For `q=2t>=exp(100)`, let `r=r_q` and `A` be the dominant `alpha=9/4`
gamma-log saddle parameters. Then

```text
|g_5'''(t)| <= 54/(t^2 r^2) <= 27/(1250 t^2).       (1)
```

Thus the `m=1`, `beta=5` term costs at most `zeta_5=27/1250` in the
H921/H920 transfer budget. This direct third-log-derivative estimate bypasses
H923's stronger sufficient requirement of separate `t^-j` bounds for all
three derivatives of `epsilon`.

## Exact Decomposition

The `m=1` part of the Xi moment is exactly

```text
2 pi^2 I_9(t)-3 pi I_5(t)
 =2 pi^2 I_9(t)[1-c_5 R_5(t)].                       (2)
```

Under `x=pi exp(r)`, `r=4U`, and `q=2t`, put

```text
J_alpha(q)=int_pi^infty log(x/pi)^q x^(alpha-1) exp(-x) dx.
```

Since `I_beta=4^(-q-1) pi^(-beta/4) J_(beta/4)`, (2) is equivalently
proportional to

```text
J_(9/4)(q)-(3/2)J_(5/4)(q).                         (3)
```

In the dominant `alpha=9/4` law, (2) multiplies the density by

```text
h(r)=1-c_5 exp(-r).
```

Because `r>0` and `pi>3`, `1/2<h<=1`.

## Normalized Saddle Law

Let `P` be the dominant gamma-log law in `z=log r`, put
`W=sqrt(A)(z-z_q)`, and let `Q` be its corrected version:

```text
dQ=h dP/E_P[h].                                      (4)
```

H1313 proves for every `r>=5/2`:

```text
M_2=E_P[W^2]<=26/25,
M_4=E_P[W^4]<=7/2,
M_6=E_P[W^6]<=21.                                    (5)
```

Moreover

```text
h'(W)=c_5 r exp(-r)/sqrt(A),
Lip(h)<=1/(2sqrt(A)),                                 (6)
```

using `c_5<1/2` and `r exp(-r)<=1/e<1`.

## Covariance Perturbation

For any real random variable `X`, any `L`-Lipschitz `h`, and any
square-integrable `f`, an independent-copy argument and Cauchy-Schwarz give

```text
Var(h(X))<=L^2 Var(X),
|Cov(f(X),h(X))|<=L sqrt(Var(f(X))Var(X)).             (7)
```

Write `m_j=E_P[W^j]`, `n_j=E_Q[W^j]`, and `delta_j=n_j-m_j`.
Equations (4)-(7) imply

```text
|delta_1|<=M_2/sqrt(A),
|delta_2|<=sqrt(M_2 M_4)/sqrt(A),
|delta_3|<=sqrt(M_2 M_6)/sqrt(A).                     (8)
```

Also

```text
|m_1|<=sqrt(M_2),  m_2<=M_2,
|n_1|<=sqrt(2M_2), n_2<=2M_2.                         (9)
```

For `kappa_3(X)=E[X^3]-3E[X]E[X^2]+2E[X]^3`, (8)-(9) yield

```text
|kappa_3(Q)-kappa_3(P)|<=C_3/sqrt(A),

C_3=sqrt(M_2 M_6)
    +3(2M_2^2+M_2 sqrt(M_4))
    +2M_2^2(3+sqrt(2)).                               (10)
```

The rational upper bounds

```text
sqrt(M_2 M_6)<24/5, sqrt(M_4)<15/8, sqrt(2)<10/7
```

give `C_3<467591/17500<27`.

## Third Derivative

The third derivative of a log-Laplace transform is its third cumulant. Hence

```text
d^3/dq^3 log(1-c_5R_5)
 =kappa_(3,Q)(z)-kappa_(3,P)(z)
 =A^(-3/2)[kappa_(3,Q)(W)-kappa_(3,P)(W)],            (11)
```

so its absolute value is at most `27/A^2`. Since `q=2t`, three
differentiations contribute `8`, and `A>=qr=2tr`. Therefore

```text
|g_5'''(t)|<=8*27/A^2<=54/(t^2r^2).                  (12)
```

H922m gives `r>=log(q)/2>=50` for `q>=exp(100)`, proving (1).

## Scope

There is no unquantified remainder in (1). The numerical scout
`tools/rh_h1316_beta5_transfer_scout.py` suggests the true sign is favorable,
but that sign is not used. This note handles only the `m=1`, `beta=5`
correction; the `m>=2` tail still requires its own explicit budget.
