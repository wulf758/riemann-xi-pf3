# H1319 Full-Xi Certificate on the First Compact Segment

Classification: `h1319_full_xi_segment_r125_to_4_certified_gap_4_to_71_open`

## Certified result

Let

```text
q(r)=pi*r*exp(r)-(9/4)r-1,
t(r)=q(r)/2,
A(r)=r[pi*exp(r)(r+1)-9/4].
```

The canonical Arb run in
`research/riemann/h1319_full_xi_adaptive_rbox_cover.json` certifies

```text
t(r)^2 kappa_Xi'''(t(r)) >= -3/4
```

on every `r` in `[323/100,4]`.  The interval is partitioned into 192 exact
rational boxes.  All 192 boxes pass and there are no missing or overlapping
leaves.

The unique solution `r(125)` of `q(r)=250` lies in

```text
323/100 < r(125) < 13/4,
```

because the certificate checks `q(323/100)<250<q(13/4)` and
`q'(r)>0` on that bracket.  Therefore the proved interval contains the desired
segment `[r(125),4]`.

The worst certified quotient-free polynomial lower bound is

```text
18999.871163750544 > 0
```

on the first box `[323/100,62093/19200]`.  The corresponding, non-authoritative
quotient diagnostic has lower endpoint

```text
-0.7136080274358392 > -3/4.
```

The SHA-256 digest of the detailed JSON is

```text
7F1B0F3C01FEB9CDCEB8B6B33B65BCF467BC69AAE7DDE55CBE5149FD68826F3C
```

## Exact moment identity

In the saddle coordinate

```text
W=sqrt(A)(z-log r),
rho=r*exp(W/sqrt(A)),
```

the full Xi density is proportional to

```text
exp(-Psi_r(W)) H_r(W),
H_r(W)=sum_(m>=1)
  [m^4-c5*m^2*exp(-rho)]
  exp[-pi*(m^2-1)*exp(rho)],
c5=3/(2*pi).
```

For the fixed translation `c=-1/12`, define

```text
K_j=int_R (W-c)^j exp(-Psi_r(W)) H_r(W) dW,
C=K_3 K_0^2-3 K_1 K_2 K_0+2 K_1^3.
```

Translation invariance of the third cumulant gives exactly

```text
kappa_Xi'''(t(r))=8 C/[A^(3/2) K_0^3].
```

Since `A>0` and `K_0>0`, the target is equivalent to the polynomial inequality

```text
P(r)=8 t(r)^2 C+(3/4)A(r)^(3/2)K_0^3 >= 0.
```

The Arb engine certifies `P(r)>0` directly.  It does not infer the claim from
the displayed floating-point diagnostic.

## Rigorous enclosure method

The cancellation-prone exponent is evaluated through the exact identity

```text
h=W/sqrt(A),
delta=r*expm1(h),
Psi=pi*exp(r)*phi2(delta)
    +(pi*r*exp(r)-(9/4)r)*phi2(h),
phi2(x)=exp(x)-1-x.
```

`phi2` is enclosed by an explicit Taylor polynomial and a geometric remainder;
no subtraction of correlated Arb intervals is used.  On `|W|<=24`, the exact
`m=1` term is integrated by a rigorous second-order midpoint rule.

The positive `m>=2` central correction is bounded with the H1316 pointwise
majorant

```text
h_bar(rho)<19 exp[-3*pi*exp(rho)].
```

For `|W|>=24`, H1316 gives `0<H<1`, so the full-Xi tails are dominated by the
H1271 gamma-log tails.  The exact conservative inputs are

```text
T_0<2e-40,  T_2<2e-37,  T_4<7e-35,
T_1<=T_2/24,  T_3<=T_4/24.
```

Shifted tails use the binomial bound for `|W-c|^j`.  No finite differences are
used anywhere in the certificate.

## Explicit remaining gap

This is not the full compact certificate required by the piecewise H920
assembly.  The endpoint `r=4` corresponds only to

```text
t(4)=2*pi*exp(4)-5 ~= 338.0502940874,
```

whereas the effective H1317 tail begins at

```text
r=71,
T_71=t(71)=[71*pi*exp(71)-(9/4)71-1]/2
     ~=7.6258180846e32.
```

Thus the entire interval

```text
4 < r < 71
```

remains open.  In particular, this artifact must not be relabelled with
`full_domain_covered=true`, and it cannot populate the global compact manifest
for `[125,T_71]`.

The direct translated-moment form becomes poorly conditioned as `r` grows:
the target skewness shrinks while separate odd raw moments still have
order-one interval errors.  The next credible engine is the full-law
score/Stein form.  With

```text
N_full=Psi'-H'/H-W,
a=E[N_full],
b=E[(W^2+2)N_full],
```

Stein integration gives

```text
E[W]=-a,
E[W^3]=-b,
kappa_3(W)=-b+3a E[W^2]-2a^3.
```

Its unnormalized numerator uses only `Z`, `int W^2 p`, `int N_full p`, and
`int (W^2+2)N_full p`; it avoids the unstable subtraction of two odd raw
moments.  The pointwise score certificates already observed at
`r=10,20,50,71` suggest combining this form with Taylor models in `r` rather
than attempting exponentially thin naive r-boxes.

## Artifacts

- canonical entry point:
  `tools/rh_h1319_full_xi_adaptive_rbox_certificate.py`;
- interval engine:
  `tools/rh_h1319_full_xi_adaptive_rbox_cover.py`;
- detailed certificate:
  `research/riemann/h1319_full_xi_adaptive_rbox_cover.json`;
- generated run summary:
  `research/riemann/h1319_full_xi_adaptive_rbox_cover.md`.
