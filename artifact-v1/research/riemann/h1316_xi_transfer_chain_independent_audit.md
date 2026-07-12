# H1316 Independent Audit of the Gamma-Log to Xi Transfer

Classification: `h1316_transfer_exact_reduction_valid_h918_proof_invalid_beta5_open`

## Verdict

The algebraic transfer from the Xi kernel to the dominant gamma-log model is
correct, including all factors `2`, `4`, and `pi`.  The exact multiplicative
correction and the logarithmic differentiation formula in H923 are also
correct.

However, H918 does **not** prove the claimed uniform transfer.  One of its
left-tail inequalities is false on the interval on which it is asserted, and
its derivative-stability paragraph contains no effective absolute bounds.
Consequently H919 cannot cite H918 as a completed proof.

Two parts of the transfer can now be certified independently:

1. the full correction always lies in `(-3/(2*pi),0)`, so its logarithm has no
   denominator or branch problem;
2. the complete `m>=2` correction and its first three derivatives are bounded
   explicitly by H1316's separately audited `m`-tail lemma.

The only remaining analytic transfer term is the `m=1`, `beta=5` ratio.

## 1. Exact Xi Moment Normalization

Use the de Bruijn-Newman normalization

```text
H_0(w)=integral_0^infinity Phi(U) cos(wU) dU
      =(1/8) xi(1/2+iw/2),

Phi(U)=sum_(m>=1)
  (2*pi^2*m^4*exp(9U)-3*pi*m^2*exp(5U))
  exp(-pi*m^2*exp(4U)).
```

Putting `w=-2iz` gives the exact identity

```text
xi(1/2+z)=8 integral_0^infinity Phi(U) cosh(2zU) dU.
```

Thus, if

```text
xi(1/2+z)=sum_(n>=0) a_n z^(2n),
M_n=(2n)! a_n,
```

then

```text
M_n=8*4^n integral_0^infinity U^(2n) Phi(U) dU.       (1)
```

H908 suppresses the constant `8` and the geometric factor `4^n` when it
writes the positive moment representation.  This is harmless for the chosen
target: it adds an affine function of `t` to `kappa(t)=log M_t`, hence changes
neither `kappa'''` nor `Delta^3 kappa`.

For a real interpolation one may therefore use

```text
M_Xi(t)=8*4^t integral_0^infinity U^(2t) Phi(U) dU,
```

which agrees with (1) at the integers.

## 2. Exact Gamma-Log Change of Variables

Define

```text
I_(m,beta)(t)=integral_0^infinity
  U^(2t) exp(beta U-pi*m^2*exp(4U)) dU.
```

With

```text
q=2t,  a=pi*m^2,  alpha=beta/4,
x=a*exp(4U),
```

one has exactly

```text
I_(m,beta)(t)
 =4^(-q-1) a^(-alpha)
  integral_a^infinity
    (log(x/a))^q x^(alpha-1) exp(-x) dx.              (2)
```

Putting `x=a*exp(r)` and then `r=exp(Z)` gives the potential

```text
V_(a,alpha)(Z)=a*exp(exp(Z))-alpha*exp(Z)-Z,
```

and its saddle satisfies

```text
a*r*exp(r)=q+alpha*r+1.                              (3)
```

For the dominant term, `a=pi`, `alpha=9/4`, and `r=4U` in the original
coordinate.  From (2),

```text
d^3/dt^3 log I_(1,9)(t)=8 K_(9/4)'''(2t).            (4)
```

This verifies the factor `8`, the relation `q=2t`, and the relation `r=4U`
used in H921-H923.

## 3. Exact Correction and Its Sign

Let `c_5=3/(2*pi)`.  Factoring the dominant moment gives

```text
integral U^(2t) Phi(U) dU
 =2*pi^2 I_(1,9)(t) (1+epsilon(t)),

epsilon(t)=-c_5 I_(1,5)(t)/I_(1,9)(t)+R_tail(t),

R_tail(t)=sum_(m>=2) [
  m^4 I_(m,9)(t)/I_(1,9)(t)
  -c_5*m^2 I_(m,5)(t)/I_(1,9)(t)].                  (5)
```

All signs and powers of `m` in H923 agree with the kernel.

There is a useful global fact omitted from H923.  Under the dominant tilted
law, put `X=exp(4U)>=1`.  Then `epsilon(t)=E_t[h(X)]`, where

```text
h(X)=-c_5/X+sum_(m>=2)
  (m^4-c_5*m^2/X) exp(-pi*(m^2-1)X).                (6)
```

Every summand in the series is positive.  On the other hand,

```text
sum_(m>=2) m^4 exp(-pi*(m^2-1)X)
 <= (1/X) sum_(m>=2) m^4 exp(-pi*(m^2-1))
 < 1/(400X) < c_5/X.                               (7)
```

The first inequality follows because `X exp(-aX)` decreases for `X>=1` and
`a>=3*pi`; the last finite constant follows from the `m=2` term and a
geometric ratio bound.  Hence

```text
-c_5/X < h(X) < 0,
-c_5 < epsilon(t) < 0,
1+epsilon(t)>1-c_5>1/2                              (8)
```

for every `t>-1/2`.  Thus H923's condition `|epsilon|<=1/2` is already
proved globally, with strict slack.  No positivity of Xi zeros or other
RH-equivalent input is used.

## 4. Exact Logarithmic Transfer

Put

```text
G(t)=log(1+epsilon(t)).
```

Equations (1) and (5) give

```text
kappa_Xi'''(t)=kappa_9'''(t)+G'''(t),                (9)
```

and direct differentiation gives

```text
G''' = epsilon'''/(1+epsilon)
       -3*epsilon'*epsilon''/(1+epsilon)^2
       +2*(epsilon')^3/(1+epsilon)^3.               (10)
```

H923's bound

```text
|G'''(t)| <=
  (2B_3+12B_1B_2+16B_1^3)/t^3
```

is correct under its derivative hypotheses.  It is a sufficient condition,
not the minimal transfer statement.

For the H803 target, the exact smallest statement at each integer is

```text
Delta^3 G(n-2)
 >= log(Theta_n)-Delta^3 kappa_9(n-2),               (11)

Theta_n=C_(n-1) C_(n+1)/C_n^2,
C_n=(2n)(2n-1).
```

For the stronger H804 unit barrier, replace `log(Theta_n)` in (11) by
`log(1-1/n^2)`.  A direct lower bound for `Delta^3 G`, or merely for the
combination (10), is enough.  Separate absolute estimates for all three
derivatives of `epsilon` are convenient but stronger than necessary.

## 5. Why H918 Is Not a Proof

H918 uses, throughout its alleged moderate left tail,

```text
F_t(U_t)-F_t(U_t-d) >= c_3*(t/U_t)*d^2,
0<=d<=U_t/2,                                        (12)
```

with a fixed positive `c_3`, where

```text
F_t(U)=2t log U+9U-pi exp(4U).
```

Take `d=U_t/2`.  The saddle equation gives

```text
F_t(U_t)-F_t(U_t/2)=2t log 2+O(t/U_t),
```

whereas the right side of (12) is `c_3*t*U_t/4`.  Therefore the ratio of the
two sides without `c_3` is

```text
8 log(2)/U_t+o(1/U_t) -> 0.
```

No fixed `c_3>0` can make (12) true on the stated range.  Direct numerical
values of that ratio are

| `t` | `U_t` | ratio |
| ---: | ---: | ---: |
| `10^30` | `16.11498` | `0.33640` |
| `10^100` | `56.09838` | `0.09821` |
| `10^300` | `170.94907` | `0.03237` |

The intended tail theorem may still be true after splitting the left side
into a genuinely local range and a nonquadratic range.  The written proof
does not establish it.  Moreover, H918's derivative paragraph bounds raw
derivatives without consistently taking absolute values and does not extract
constants for derivatives of the quotient.  H923 correctly treats those
bounds as still open.

## 6. The `m>=2` Part Is Now Closed Separately

The independently audited H1316 `m`-tail lemma proves, for `r>=5/2`,

```text
0<=R_tail(t)<1/60,
|R_tail'(t)|   <= 1311/t,
|R_tail''(t)|  <= 6555/t^2,
|R_tail'''(t)| <= 48507/t^3.                        (13)
```

Its change of variables, frozen-saddle quotient derivatives, product powers,
coefficient counts, and the conversion `q=2t` were checked line by line.  The
standalone verifier reports `all_checks_pass=true`.

Thus (13) rigorously replaces the `m>=2` portion of H918.  It leaves only

```text
R_5(t)=I_(1,5)(t)/I_(1,9)(t)                       (14)
```

and its first three derivatives to complete H923.

## 7. What a Completed Transfer Would Actually Give

Suppose the remaining `R_5` bounds yield a full correction estimate

```text
G'''(t)>=-zeta/t^2.
```

The presently proved coarse gamma-log bound quoted in H923 is

```text
kappa_9'''(t)>=-35/(2t^2 U_t).
```

Even with the ideal value `zeta=0`, H920's coefficient becomes less than one
only when

```text
U_t>35/2, equivalently r=4U_t>70.                  (15)
```

This corresponds to an enormous moment threshold, far beyond the finite Xi
audit through `n=126`.  Therefore completing H923 is worthwhile, but it does
not by itself close H920.  After the transfer, one still needs either a much
sharper effective gamma-log skewness bound, a direct discrete estimate of
(11), or a correspondingly huge finite verification.

## Status Table

| item | audited status |
| --- | --- |
| Xi kernel and moment normalization | exact and valid |
| factors `8`, `q=2t`, `r=4U` | exact and valid |
| signs in `epsilon` | exact and valid |
| global `-c_5<epsilon<0` | proved here |
| H923 logarithmic derivative identity | exact conditional gate |
| H918 uniform transfer proof | invalid as written |
| H1316 `m>=2` transfer | independently verified |
| `m=1`, `beta=5` derivative transfer | open at this audit point |
| H920 finite closure | not implied by transfer alone |

## Non-Circularity

Everything audited here uses only the classical positive Xi kernel, elementary
changes of variables, moment identities, and explicit inequalities.  No
statement about the location of zeta zeros is used.  The problem is missing
effectivity/sharpness, not circular dependence on RH.

