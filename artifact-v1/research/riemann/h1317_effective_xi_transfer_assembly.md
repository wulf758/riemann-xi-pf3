# H1317 Effective Gamma-Log to Xi Transfer Assembly

Classification: `h1317_effective_xi_transfer_proved_huge_finite_gap_open`

## Result

Let `kappa_Xi(t)` be the logarithm of the positive Xi-kernel moment transform
used in H908, and let `r=r_(2t)` be the dominant `alpha=9/4` gamma-log saddle:

```text
pi*r*exp(r) = 2t + (9/4)r + 1.
```

For every `t` for which `r>=5/2`, the explicit transfer estimates assembled
below give

```text
kappa_Xi'''(t)
 >= -[70/r + 45450578214/t]/t^2.                 (1)
```

In particular, define

```text
T_71 = [71*pi*exp(71) - (9/4)*71 - 1]/2.
```

Then, for every `t>=T_71`,

```text
kappa_Xi'''(t) >= -(141/142)/t^2.                (2)
```

This is an effective transfer to the **full Xi kernel** with a coefficient
strictly smaller than one.  It proves an effective eventual H908-A/H920
pointwise input.  It does not close the finite interval below `T_71`, and it
does not prove the all-degree Jensen/total-positivity assembly or RH.

Numerically,

```text
T_71 ~= 7.6258180846 * 10^32.
```

Thus the transfer succeeds, but the resulting finite gap is enormous.

## 1. Exact Xi Correction

With

```text
I_(m,beta)(t)
 = integral_0^infinity U^(2t)
   exp(beta*U-pi*m^2*exp(4U)) dU,
c_5=3/(2*pi),
```

the full moment is

```text
M_Xi(t)=constant*I_(1,9)(t)*(1+epsilon(t)),

epsilon(t)=-c_5 R_5(t)+R_tail(t),
R_5=I_(1,5)/I_(1,9).
```

The exact normalization, the factor `q=2t`, the factor eight in third
derivatives, and the signs were independently audited in H1316.

Pointwise, with `X=exp(4U)>=1`,

```text
-c_5/X
 < -c_5/X + sum_(m>=2)
      (m^4-c_5*m^2/X) exp(-pi*(m^2-1)X)
 < 0.
```

Indeed,

```text
sum_(m>=2) m^4 exp(-pi*(m^2-1)X)
 <= X^(-1) sum_(m>=2) m^4 exp(-pi*(m^2-1))
 < 1/(400X) < c_5/X.
```

Therefore

```text
-c_5 < epsilon(t) < 0,
1+epsilon(t) > 1-c_5 > 1/2.                    (3)
```

There is no logarithmic branch or denominator problem.

## 2. Explicit `m>=2` Tail

H1316 proves, using the global H1313 moments and summing `m` before
integration,

```text
0 <= R_tail < 1/60,
|R_tail'|   <=  1311/t,
|R_tail''|  <=  6555/t^2,
|R_tail'''| <= 48507/t^3.                      (4)
```

This replaces the invalid non-effective H918 tail argument.

## 3. Explicit Derivatives of `R_5`

Under the dominant gamma-log law, put

```text
rho=exp(Z),
f(Z)=exp(-rho),
W=sqrt(A)*(Z-z_q),
S=W/sqrt(A),
x=q/r>30.
```

For derivatives at a fixed `q`, freeze the saddle center and replace `f` by

```text
g(Z)=f(Z)-f(z_q).
```

The constant disappears from every quotient derivative.  Split at

```text
G={W>=-sqrt(A)/r}.
```

On `G`, `rho>=r-1`.  Since `rho*exp(-rho)` decreases for `rho>=1`,

```text
|g| <= (12r/x)*|W|/sqrt(A).
```

On `G^c`, use `|g|<=1` and H1313's sixth moment.  For

```text
n_j=E[g*S^j],
D_j=E[S^j],
```

the resulting bounds are, for `0<=j<=3`,

```text
|n_j| <= q^(-j)[48*x^((j-3)/2)+21*x^(j-3)],
|D_j| <= 2*q^(-j)*x^(j/2).                     (5)
```

Every monomial of total derivative order `k<=3` in the exact quotient rule
is therefore at most `69*2^s/q^k`, because all residual powers of `x` are
nonpositive.  Counting the quotient coefficients gives, after `q=2t`,

```text
|R_5'|   <=  207/t,
|R_5''|  <= 1035/t^2,
|R_5'''| <= 7659/t^3.                          (6)
```

The separate direct beta-5 argument also proves the much sharper third-log
bound `|d_t^3 log(1-c_5 R_5)|<=54/(t^2r^2)`, but (6) is convenient for the
single H923 correction.

## 4. H923 Assembly

From `c_5<1/2`, (4), and (6), one may take

```text
B_1=1415,
B_2=7073,
B_3=52337
```

in

```text
|epsilon^(j)(t)| <= B_j/t^j.
```

Using (3), H923 gives

```text
|d_t^3 log(1+epsilon)| <= B_Xi/t^3,

B_Xi = 2B_3 + 12B_1B_2 + 16B_1^3
     = 45450578214.                              (7)
```

## 5. Dominant Gamma-Log Term

H1315 and H922 give

```text
kappa_9'''(t) >= -70/(t^2 r).                    (8)
```

Combining (7) and (8) proves (1).

For `r>=71`,

```text
70/r <= 70/71=140/142.
```

Also `T_71>142 B_Xi`; this follows already from `pi>3` and `e>2` by exact
integer arithmetic.  Hence, for `t>=T_71`,

```text
B_Xi/t <= 1/142.
```

Adding the two estimates proves (2).

## 6. What This Gives — And What It Does Not

H920's geometric threshold for `a=141/142` is below `567`.  Since `T_71` is
far larger, H908-A follows for integers approximately

```text
n >= ceil(T_71+2).
```

The existing finite Xi audit reaches only small indices (about `126` in the
current reports).  The unverified finite interval therefore still extends to
roughly `7.6*10^32`.

So the transfer itself is now effective.  What remains is not a hidden
asymptotic remainder, but:

```text
1. sharpen the dominant gamma-log coefficient enough to lower T_71,
   or prove the discrete finite range by another method;
2. prove the all-degree Jensen/total-positivity implication required by H907.
```

Neither item is supplied by H1317.

## 7. Correction to the Earlier Record

H918's written moderate-left-tail quadratic bound is false on its stated
whole interval, so H918/H919 must not be cited as the proof of this transfer.
H1316-H1317 replace that step with explicit quotient estimates.  This
correction does not invalidate the transfer statement; it replaces an invalid
proof with an effective one.

## Verifier

Run

```text
python tools/rh_h1317_effective_xi_transfer_assembly.py
```

The verifier checks all constants, the correction-sign estimate, the saddle
threshold comparison, and the dependency artifacts.
