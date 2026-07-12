# H1319 Full-Xi Compact Scaled-Cumulant Scout

Classification: `h1319_full_xi_compact_hybrid_dense_scout_not_proof`

## Question

On the compact saddle interval corresponding to

```text
125 <= t,
r(t) <= 71,
```

does the full Xi kernel numerically satisfy

```text
t^2 kappa_Xi'''(t) >= -3/4?                            (1)
```

## Numerical Answer

Yes on every sampled point, with substantial room. The sampled minimum is the
left endpoint

```text
t = 125,
r = 3.2352157911731031290131595418145375261332828802432...
```

where the converged multiprecision values are

```text
t^2 kappa_dominant''' = -0.5414206382046353836753562944014591050110...
t^2 correction'''     = +0.0005271015080698127538602820503536893717...
t^2 kappa_Xi'''       = -0.5408935366965655709214960123511054156393...
```

Thus the sampled minimum has margin

```text
0.2091064633034344290785039876488945843607...
```

above `-3/4`.

This is a finite-window pointwise diagnostic, not a proof of (1).

## Artifacts

Canonical hybrid scout:

```text
tools/rh_h1319_full_xi_compact_hybrid_scout.py
```

Shared stable formulas and preliminary low-`r` engine:

```text
tools/rh_h1319_full_xi_compact_scaled_cumulant_scout.py
```

Canonical output:

```text
research/riemann/h1319_full_xi_compact_scaled_cumulant_scout.json
```

The canonical invocation uses

```text
python tools/rh_h1319_full_xi_compact_hybrid_scout.py --mp-threshold 20
```

## Direct Cumulant Reconstruction

No finite differences are used. For the dominant normalized saddle law, the
scout integrates raw `W` moments through order three and forms

```text
kappa_3(W)=E[W^3]-3E[W]E[W^2]+2E[W]^3.
```

Since `q=2t` and `W=sqrt(A)(Z-z_q)`, the required derivative is

```text
t^2 kappa'''(t)=8 t^2 kappa_3(W)/A^(3/2).              (2)
```

The full Xi kernel is written as a multiplicative correction

```text
H(W)=1+delta(W),
```

where

```text
delta = -3/(2x)
        + sum_{m>=2} (m^4-3m^2/(2x)) exp(-x(m^2-1)),
x = pi exp(r_local).
```

The scout integrates the correction moments separately and applies exact
raw-moment update identities. This retains correction terms too small to be
represented by first forming `1+delta` in floating point.

## Stable Saddle Coordinate

Let

```text
s  = W/sqrt(A),
dr = r(expm1(s)),
P  = pi exp(r).
```

The centered potential is evaluated with the exact cancellation identity

```text
Psi(W)
 = P(expm1(dr)-dr)
   +(P-9/4)r(expm1(s)-s).                              (3)
```

Both `expm1(x)-x` terms are evaluated by a local series near zero. Formula
(3) avoids subtracting quantities of saddle size.

## Sampling And Minimum Zoom

The canonical run contains `427` distinct sample points obtained from:

```text
161 logarithmic points on the whole compact interval;
161 linear points on the low-r region;
explicit anchor points;
4 adaptive zoom rounds around the sampled minimum.
```

Every zoom round selects the left endpoint again. The final zoom bracket has
width below `2e-6` in `r`, and the sampled minimum remains exactly
`r(t=125)`.

After the multiprecision repair described below, the sampled full scaled
cumulant is nondecreasing over the entire ordered grid.

## Hybrid Precision

For `r<20`, stable segmented Gauss-Legendre quadrature uses

```text
float64, W=24, m_cutoff=12, order=32 per segment.
```

For all `71` sampled points with `r>=20`, the moments are recomputed using

```text
mpmath dps=60, W=20, m_cutoff=8, order=24 per segment.
```

This replacement is necessary because the odd normalized cumulant becomes
comparable to binary64 roundoff at large `r`. An all-float64 preliminary run
developed three artificial downward steps beyond `r≈60`; those steps vanish
under direct multiprecision integration. The preliminary high-`r` rows are
therefore noncanonical.

At the engine seam `r=20`, the float64 and multiprecision scaled values differ
by only

```text
7.11e-11.
```

## Convergence At The Minimum

At `r(t=125)`, the scout varies precision, finite window, kernel cutoff, and
quadrature order:

```text
dps:       50, 70, 90
W:         18, 24, 28
m_cutoff:   8, 12, 16
order:     24, 32
```

The coarse `dps=50, W=18, m<=8, order=24` value differs from the reference
`dps=90, W=28, m<=16, order=32` value by approximately

```text
9.9e-33.
```

Changing `W=18` to `W=24` accounts for that difference; subsequent cutoff,
order, and precision changes agree to still greater accuracy.

## High-Endpoint Check

At `r=71`, the reference calculation gives

```text
t^2 kappa_Xi''' = -0.02815822187928669410150891632372665486673...,
```

with margin about `0.7218417781` above `-3/4`. The correction contribution is
positive and of size about `5.35e-64` in this scaling.

The coarse high-end configuration

```text
dps=50, W=18, m<=8, order=20
```

differs from the `dps=90, W=28, m<=16, order=32` reference by only

```text
2.92e-23.
```

## Interpretation

The scan strongly supports all three numerical hypotheses:

```text
1. the target (1) has comfortable room on the compact interval;
2. the worst point is the left endpoint t=125;
3. the full-Xi correction is favorable on the sampled grid.
```

The important margin is about `0.2091`, not a near-zero numerical effect.

## Why This Is Not A Proof

The computation still has the following gaps:

```text
finite W windows rather than proved analytic tails;
sampled r values rather than interval r-box coverage;
floating/multiprecision point arithmetic rather than Arb enclosures;
no theorem excluding a narrow unsampled interior dip;
no certified bound across the float64/mpmath seam.
```

The corresponding rigorous H1319 target is an adaptive Arb cover on
`[r(t=125),71]`, using interval `W` quadrature and analytic left/right tails.

## Scope

This report diagnoses the compact full-Xi third-cumulant inequality only. It
does not close H920, the all-degree Jensen-hyperbolicity step, or the Riemann
Hypothesis.
