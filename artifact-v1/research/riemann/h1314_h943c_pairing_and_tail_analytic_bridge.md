# H1314/H943-C Pairing And Tail Analytic Bridge

Classification: `h1314_h943c_pairing_tail_bridge_proved`

## Statement

This note supplies the analytic bridge used by the reproducible H1314
uniform-cap certificate. It proves:

```text
R(t) = K_e(t)-K_o(t)tanh(O(t)),
E[g(W)K(W)] = E[g(W)R(|W|)]       for every even g,
```

and, for the exact gamma-log family with `r>=5/2`,

```text
R(t)-1/2
 <=(1/2)exp(r(exp(z/mu)-1)+3z/mu-2s(sinh(z)-z)),
 z=a t>=1.
```

The exponent is at most `-17/4` at `z=1` and is strictly decreasing
thereafter. Consequently

```text
R(t)-1/2 < (1/2)exp(-17/4) < 1/100       (z>=1).
```

No polynomial from the separate strong H946 branch is used.

## 1. Pairing Identity

Let

```text
p_+ = exp(-Psi(t)),       p_- = exp(-Psi(-t)),
K_+ = K(t),               K_- = K(-t),
O   = (Psi(t)-Psi(-t))/2.
```

Then

```text
p_+/p_- = exp(-2O).
```

Define the curvature ratio under the paired radial mass by

```text
R(t)=(K_+p_+ + K_-p_-)/(p_+ + p_-)
    =(K_+exp(-2O)+K_-)/(1+exp(-2O)).                (1)
```

Put

```text
K_e=(K_++K_-)/2,
K_o=(K_+-K_-)/2.
```

Since

```text
tanh(O)=(1-exp(-2O))/(1+exp(-2O)),
```

substitution in (1) gives exactly

```text
R=K_e-K_o tanh(O).                                  (2)
```

The standalone checker simplifies the difference between the two sides to
zero symbolically.

## 2. Even Expectations Are Preserved

Let `Z_P=int_R exp(-Psi(w))dw`, and let `g` be even. Pairing `w=t` and
`w=-t` gives

```text
E[g(W)K(W)]
 =Z_P^(-1) int_0^infinity
   g(t)(K_+p_+ + K_-p_-)dt

 =Z_P^(-1) int_0^infinity
   g(t)R(t)(p_+ + p_-)dt

 =E[g(W)R(|W|)].                                    (3)
```

Thus a pointwise bound on `R(t)` transfers directly to both even H943-C
budgets. There is no lost factor of two: the same paired mass `p_++p_-`
appears in the normalizer on both sides.

## 3. Dimensionless Formulas

Use the exact exponential-mixture normalization

```text
E_nu[X]=1,
mu=B/A,
a=B/A^(3/2),
s=1/a^2=A^3/B^2,
z=a t.
```

The normalized potential is

```text
Psi(t)=s E[(exp(zX)-1-zX)/X^2].                     (4)
```

Since `N(t)=Psi'(t)-t=a t^2K(t)`, differentiation of (4) gives

```text
K_+=E[(exp(zX)-1-zX)/(z^2X)],
K_-=E[(exp(-zX)-1+zX)/(z^2X)].                      (5)
```

Taking the even and odd halves of (5), and taking the odd half of (4), yields

```text
K_e=E[(cosh(zX)-1)/(z^2X)],
K_o=E[(sinh(zX)-zX)/(z^2X)],
O  =s E[(sinh(zX)-zX)/X^2].                         (6)
```

These identities have continuous values at `z=0`.

## 4. Lower Bound For The Odd Potential

For `z>=0`, expand the positive series in (6):

```text
O=s sum_(k>=1)
    z^(2k+1) E[X^(2k-1)]/(2k+1)!.                  (7)
```

Because `X>0`, `E[X]=1`, and `x -> x^n` is convex for every integer `n>=1`,
Jensen gives

```text
E[X^n] >= E[X]^n = 1.
```

Termwise comparison in (7) therefore proves

```text
O >= s(sinh(z)-z) >= s z^3/6.                       (8)
```

This is the exact lower bound used in both the compact and tail regimes.

## 5. The Left Curvature Satisfies `K_-<=1/2`

For `x>=0`, set

```text
d(x)=1-x+x^2/2-exp(-x).
```

Then

```text
d(0)=d'(0)=0,
d''(x)=1-exp(-x)>=0.
```

Hence `d(x)>=0`, equivalently

```text
exp(-x)-1+x <= x^2/2.                               (9)
```

Apply (9) pointwise with `x=zX` in (5):

```text
K_- <= E[X]/2 = 1/2.                                (10)
```

This also follows from the monotonicity of `V'''` in the H943-B integral
representation, but (9) is a direct proof in the mixture coordinates.

## 6. Upper Bound For `K_+`

Let `z_0=log r` denote the saddle coordinate and put

```text
u=t/sqrt(A)=z/mu.
```

H943-B proves the exact differential inequality

```text
0 < d/dy log V'''(y) <= exp(y)+3.                   (11)
```

For `0<=theta<=1`, integrate (11) from `z_0` to `z_0+theta*u`:

```text
log(V'''(z_0+theta*u)/B)
 <= r(exp(theta*u)-1)+3theta*u
 <= r(exp(u)-1)+3u =: E.                            (12)
```

The exact H943 integral formula is

```text
K_+=int_0^1 (1-theta)
     V'''(z_0+theta*u)/B dtheta.
```

Using (12),

```text
K_+ <= exp(E) int_0^1(1-theta)dtheta
     = (1/2)exp(E),
E=r(exp(z/mu)-1)+3z/mu.                             (13)
```

## 7. Pair-Tail Bound

Subtract `1/2` in the paired formula (1):

```text
R-1/2
 =[(K_+-1/2)exp(-2O)+(K_--1/2)]/(1+exp(-2O)).
```

By (10), the second numerator is nonpositive. Also the denominator is at
least one and `K_+-1/2<=K_+`. Therefore (8) and (13) give

```text
R-1/2
 <=K_+exp(-2O)
 <=(1/2)exp(H(z)),                                   (14)

H(z)=r(exp(z/mu)-1)+3z/mu-2s(sinh(z)-z).
```

This proves the pair-tail input that had previously only been recorded as a
displayed formula.

## 8. Exact Junction Bound

The parameter certificate proves

```text
s>=18,
mu>=r+3/2>=4.                                       (15)
```

The elementary integral inequality

```text
log(1+x)=int_0^x dt/(1+t) >= x/(1+x)
```

with `x=1/r` gives

```text
1/mu < 1/(r+1) <= log(1+1/r),
r(exp(1/mu)-1) <= 1.                                (16)
```

Moreover, the positive Taylor series gives exactly

```text
sinh(1)-1
 =1/6+1/120+1/5040+...
 >1/6.                                              (17)
```

Combining (15)-(17),

```text
H(1)
 <=1+3/4-2*18*(1/6)
 =-17/4.                                            (18)
```

The checker regenerates the strict lower partial sum

```text
1/6+1/120=7/40>1/6
```

rather than treating (17) as a hardcoded Boolean.

## 9. The Tail Exponent Is Decreasing

Differentiate:

```text
H'(z)=(r/mu)exp(z/mu)+3/mu-2s(cosh(z)-1).           (19)
```

From (15), for `z>=1`,

```text
(r/mu)exp(z/mu)+3/mu
 <=exp(z/4)+3/4
 <=(7/4)exp(z).                                     (20)
```

It remains to justify the lower bound used on the negative term. Put
`Y=exp(z)`. Since `z>=1` and `e>2`, `Y>=2`. Exact algebra gives

```text
cosh(z)-1-exp(z)/8
 =(3Y-2)(Y-2)/(8Y) >= 0.                            (21)
```

Therefore, using `s>=18`,

```text
2s(cosh(z)-1) >= (9/2)exp(z).                       (22)
```

Equations (19)-(22) yield the explicit strict bound

```text
H'(z) <= -(11/4)exp(z) < 0.                         (23)
```

Thus (18) holds throughout the tail. Finally,

```text
(1/2)exp(-17/4)
 =0.0071321169544996276366... < 1/100,
```

with a 256-bit Arb enclosure regenerated by the checker.

## 10. Standalone Checker

Run

```text
python tools/rh_h1314_h943c_pairing_tail_bridge.py
```

Canonical output:

```text
research/riemann/h1314_h943c_pairing_tail_bridge.json
```

It reports `all_checks_pass=true` and checks:

```text
pairing algebra residual = 0,
paired even-expectation residual = 0,
dimensionless Ke/Ko/O residuals = 0,
K_- scalar derivative identities,
exact sinh partial sum,
exact cosh factorization,
H(1) and H' constants,
parameter-certificate status and hash.
```

## Scope

This closes the analytic pairing and tail bridge for the uniform H1314
certificate. It uses no unrecorded coefficient from the stronger H946
corrected-pair branch. It remains a gamma-log result; it is not by itself an
Xi transfer or a proof of the Riemann Hypothesis.
