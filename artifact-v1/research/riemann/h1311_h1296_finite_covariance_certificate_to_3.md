# H1311 H1296 Finite Covariance Certificate on `5/2 <= r <= 3`

Classification: `h1311_h1296_finite_covariance_certificate_to_3_proved`

This is an Arb computer-assisted certificate for the exact H1296 covariance
criterion for `k=1,2,3` on `5/2<=r<=3`.  It combines the H1305 analytic
component-tail bound with finite paired-score quadrature.

It does not cover `3<r<49/5` and does not prove RH.

## 1. Finite paired-score identity

For a cutoff `L`, put

```text
N_< = sum_{m<L} nu_m G_{lambda_m},
S_< = sum_{m<L} nu_m lambda_m.
```

The finite part of the H1305 margin is

```text
F_{k,L}=Cov(W^(2k),N_<)-k S_< E[W^(2k)].
```

Define

```text
Q_< = N_< -(S_</2) W Psi_r'(W).
```

Stein integration gives

```text
Cov(W^(2k),W Psi_r'(W))=2k E[W^(2k)],
```

and therefore the exact, better-conditioned identity

```text
F_{k,L}=Cov(W^(2k),Q_<).                            (1)
```

The certificate pairs `w=+t` and `w=-t`, evaluates the four required
covariance integrals on `0<=t<=16`, and adds rigorous one-sided tails.

## 2. Rigorous quadrature and tails

Each central cell uses the midpoint value.  An `arb_series` enclosure of the
second derivative on the full cell supplies the standard midpoint remainder

```text
|integral_cell f-h f(mid)| <= h^3 sup|f''|/24.
```

The right tail uses

```text
Psi_r(w)>=w^2/2,       w>=0,
```

and Gaussian-power recurrences.  The left tail uses the convex tangent at
`w=-16`.  For the finite score,

```text
0<=N_<(-t)<=t,
0<=N_<(t)<=N(t),
```

while integration by parts bounds the `W Psi'` contribution.  All omitted
components `m>=L` are covered by the one-term analytic H1305 bound.

## 3. Certified covers

The verified driver is

```text
tools/rh_h1306d_paired_cover_pilot.py
```

which imports the paired implementation

```text
tools/rh_h1306c_paired_finite_covariance_pilot.py
```

and the corrected canonical special-function entrypoint

```text
tools/rh_h1306b_finite_covariance_certificate_pilot.py.
```

The following two runs were executed:

```text
python tools/rh_h1306d_paired_cover_pilot.py \
  --r-lo 5/2 --r-hi 13/5 --step 1/2000 \
  --cutoff 16 --subdivisions 400 --dps 80

python tools/rh_h1306d_paired_cover_pilot.py \
  --r-lo 13/5 --r-hi 3 --step 1/2000 \
  --cutoff 17 --subdivisions 400 --dps 80
```

Both use `W=16`.  The cutoffs satisfy `L>=ceil(4r+5)` throughout their
respective bands, as required by H1305.

### Band `5/2 <= r <= 13/5`

```text
box count: 200
failure count: 0

k=1: worst lower margin 0.00237090031455753846410066
k=2: worst lower margin 0.0513153222183327610201118
k=3: worst lower margin 0.775840630835235307395453
```

All three worst boxes are `[5199/2000,13/5]`.

### Band `13/5 <= r <= 3`

```text
box count: 800
failure count: 0

k=1: worst lower margin 0.000220300394434491459891042
k=2: worst lower margin 0.0207561735347356311363413
k=3: worst lower margin 0.336031646843583548074169
```

All three worst boxes are `[5999/2000,3]`.

Thus every one of the 1000 rational `r`-boxes passes for `k=1,2,3`.

## 4. Important implementation correction

The original file

```text
tools/rh_h1306_finite_covariance_certificate_pilot.py
```

reversed the two arguments of python-flint's Hurwitz zeta in its raw
canonical helper.  It is retained only as an auditable first version and must
not be used directly.  H1306b corrects both the cumulant and trigamma calls;
H1306c and H1306d import that correction before running.  The numbers above
come from the corrected path.

## 5. Consequence and boundary

Combining (1), the finite Arb cover, and the H1305 lower bound for the omitted
tail proves

```text
Cov(W^(2k),N) >= k a E[W^(2k)],   k=1,2,3,
```

for every `5/2<=r<=3`.  Hence the exact moments `M_2,M_4,M_6` are
nonincreasing on this interval.  Together with the H1282 endpoint anchor,
the H946 moment targets hold on `5/2<=r<=3`.

A direct box at `r=4` remains positive, but the naive non-Taylor interval in
`r` is too wide to certify a neighborhood efficiently.  The next route
should use an `r`-Taylor model, the scaled hyperbolic paired score, or the
H1307 radial comparison; brute-force refinement is not the chosen route.

## Status

Proved: H1296 and the resulting even-moment monotonicity on `5/2<=r<=3`.

Still open: the remaining compact interval `3<r<49/5`, H946 globally unless
that gap is closed, the Xi-kernel transfer, and the Riemann Hypothesis.
