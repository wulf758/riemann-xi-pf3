# H1319 Full-Xi Fixed-r Arb Pilot

Classification: `h1319_full_xi_fixed_r_point_certified_not_r_cover`

## Result

At the rational point

```text
r = 647/200 = 3.235,
t = 124.9637667370406...,
```

the exact full-Xi third cumulant satisfies the rigorous Arb enclosure

```text
-0.6878928068
 < t^2 kappa_Xi'''(t)
 < -0.3935118285.
```

In particular,

```text
t^2 kappa_Xi'''(t) + 3/4 > 0.0621071932.
```

This certifies the proposed `-3/(4t^2)` target at one point just below
`t=125`. It is not yet a proof on an interval of `r`.

## Exact Representation

Under the normalized dominant gamma-log law in `W`, the full Xi kernel is the
positive multiplicative tilt

```text
H(X) = sum_(m>=1)
  (m^4-3m^2/(2X)) exp(-X(m^2-1)),
X = pi exp(r exp(W/sqrt(A))).
```

For

```text
J_j = integral W^j exp(-Psi_r(W)) H(X(W)) dW,
```

the verifier forms

```text
kappa_Xi'''(t)
 = 8 A^(-3/2)
   [J_3/J_0 - 3 J_1 J_2/J_0^2 + 2 J_1^3/J_0^3].
```

No finite difference is used.

## Certificate Shape

- `|W|<=24`: first-order Arb boxes, 50,000 subdivisions.
- Xi arithmetic sum: exact interval evaluation through `m=8`.
- Omitted positive arithmetic tail: `<4.65e-106` uniformly.
- `|W|>=24`: H1272 dominant tail bounds, using the proved pointwise
  inequality `0<H<1`.

The canonical command is

```text
python tools/rh_h1319_full_xi_fixed_r_arb_pilot.py --subdivisions 50000
```

and the generated report is

```text
research/riemann/h1319_full_xi_fixed_r_arb_pilot.json.
```

## Resolution Diagnostic

At 2,000 and 20,000 subdivisions the interval enclosure was too wide to prove
the target. Those runs were inconclusive, not counterexamples. At 50,000
subdivisions the lower margin became strictly positive.

## Remaining Step

Replace the fixed rational value of `r` by adaptive rational boxes covering

```text
[647/200,71].
```

The first-order method is expected to require narrow boxes near the left edge;
Taylor or midpoint remainder forms may be preferable for a full cover.
