# H1316 Full Gamma-Log to Xi Transfer Assembly

Classification: `h1316_full_xi_h923_transfer_proved_huge_finite_gap_open`

## Result

For the full Xi moment correction, the H923 hypotheses are now explicit. For

```text
t >= T_Xi := exp(100)/2,
```

one has

```text
kappa_Xi'''(t) >= -3/(4t^2).                         (1)
```

The full multiplicative correction has the derivative constants

```text
B_1=1415,
B_2=7073,
B_3=52337,
B_Xi=2B_3+12B_1B_2+16B_1^3
    =45450578214.                                    (2)
```

This completes the effective analytic transfer from the dominant gamma-log
model to the full Xi kernel on `t>=exp(100)/2`. It does **not** close the
finite interval below that threshold: H920 leaves roughly `1.344e43`
integer moment indices, while the existing finite audit stops at `126`.

## Exact Multiplicative Correction

Let

```text
I_(m,beta)(t)
 =int_0^infinity U^(2t)
   exp(beta U-pi m^2 exp(4U)) dU,
c_5=3/(2pi).
```

After suppressing only positive affine-in-`t` factors, which do not affect a
third derivative, the full Xi moment is exactly

```text
M_Xi(t)=2pi^2 I_(1,9)(t)[1+epsilon(t)],

epsilon(t)=-c_5 R_5(t)+R_tail(t),
R_5(t)=I_(1,5)(t)/I_(1,9)(t),

R_tail(t)=sum_(m>=2)[
  m^4 I_(m,9)(t)/I_(1,9)(t)
 -c_5 m^2 I_(m,5)(t)/I_(1,9)(t)].                  (3)
```

All signs and constants in (3) follow directly from the classical positive
Xi kernel.

## Denominator Safety

Under the dominant law put `X=exp(4U)>=1`. Then `epsilon=E[h(X)]`, where

```text
h(X)=-c_5/X+sum_(m>=2)
 (m^4-c_5m^2/X)exp(-pi(m^2-1)X).                    (4)
```

Every summand is positive. Conversely,

```text
sum_(m>=2)m^4 exp(-pi(m^2-1)X)<1/(400X)<c_5/X.      (5)
```

For completeness, the first inequality follows by noting that
`X exp(-aX)` decreases for `X>=1`, and then applying the geometric estimate

```text
sum_(m>=2)m^4 exp(-pi(m^2-1))
 <=16exp(-3pi)/(1-16exp(-4pi))<1/400.
```

The last scalar inequality follows from `pi>3`, the Taylor lower bounds
`exp(9)>7000`, `exp(12)>1000`, and exact rational arithmetic. Also
`1/400<c_5` follows from `pi<4`.

Thus

```text
-c_5/X<h(X)<0,
-c_5<epsilon(t)<0,
1+epsilon(t)>1-c_5>1/2.                             (6)
```

In particular `|epsilon|<1/2` globally; the logarithm has no branch or
small-denominator issue.

## Full Derivative Box

The H1316 beta-5 lemma proves, for `q=2t>=exp(100)`,

```text
|R_5'|<=207/t,
|R_5''|<=1035/t^2,
|R_5'''|<=7659/t^3.                                 (7)
```

The independent H1316 infinite-tail lemma proves already for every saddle
`r>=5/2`:

```text
|R_tail'|<=1311/t,
|R_tail''|<=6555/t^2,
|R_tail'''|<=48507/t^3.                             (8)
```

Since `c_5<1/2`, (3), (7), and (8) give

```text
|epsilon'|   <=(1311+207/2)/t       <1415/t,
|epsilon''|  <=(6555+1035/2)/t^2    <7073/t^2,
|epsilon'''| <=(48507+7659/2)/t^3   <52337/t^3.      (9)
```

This proves the constants in (2). Applying H923 to

```text
G(t)=log(1+epsilon(t))
```

gives the explicit full-Xi correction bound

```text
G'''(t)>=-45450578214/t^3.                          (10)
```

## Combining With The Dominant Model

H922m and H922 give, for `q=2t>=exp(100)`,

```text
kappa_9'''(t)>=-70/(r t^2),                         (11)
```

where the dominant saddle satisfies

```text
q=pi r exp(r)-(9/4)r-1.                             (12)
```

The saddle function is increasing. At `r=94`,

```text
pi*94*exp(94)<exp(100),
```

because `pi<22/7` and the Taylor polynomial through degree seven gives
`exp(6)>2101/7>94*(22/7)`. Hence

```text
q>=exp(100) => r>94.                                (13)
```

Also, for `t>=exp(100)/2`,

```text
45450578214/t<1/188,                                (14)
```

because `exp(100)>2^100>376*45450578214`.

Finally, (10)-(14) yield

```text
kappa_Xi'''(t)
 >=-[70/r+45450578214/t]/t^2
 >-[35/47+1/188]/t^2
 =-3/(4t^2).
```

This proves (1) with the exact H920 constants

```text
a=3/4,
T=exp(100)/2.                                       (15)
```

## What H920 Then Gives

H920's effective closure lemma gives

```text
N_eff
 =ceil(max(3,T+2,2/(1-sqrt(3/4))))
 =ceil(exp(100)/2+2)
 ~=1.3440585709080678*10^43.                        (16)
```

Therefore the H803/H908 moment-curvature inequality is proved analytically
for all integers `n>=N_eff`.

But the existing finite certificate covers only `2<=n<=126`. The uncovered
finite interval is therefore

```text
127 <= n < ceil(exp(100)/2+2).                      (17)
```

This is an enormous finite hole, not a cosmetic one. Closing it by direct
enumeration is presently unrealistic. The next useful move is to sharpen the
dominant gamma-log threshold/constants drastically or prove the discrete
moment-curvature inequality directly on a much larger symbolic range.

## Dependency And Audit Status

- `h1316_beta5_epsilon_derivative_box.md` proves (7).
- `h1316_xi_mtail_explicit_transfer.md` proves (8) and has a standalone
  verifier reporting `all_checks_pass=true`.
- `h1316_xi_transfer_chain_independent_audit.md` independently checks the Xi
  normalization, signs, denominator safety, and identifies why the old H918
  report was not a valid effective proof.
- `h1316_full_xi_transfer_assembly.py` checks all constant additions, the
  H923 polynomial, the `r>94` scalar comparison, and the final `3/4` budget.

This completes the requested transfer. It does not prove the Riemann
Hypothesis, close the huge finite interval (17), or solve the later global
Jensen/total-positivity assembly.
