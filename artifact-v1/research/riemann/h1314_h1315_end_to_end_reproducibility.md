# H1314-H1315 End-to-End Reproducibility Certificate

Classification: `h1314_h1315_end_to_end_materialized_and_verified`

## Result

The calculations used by the accepted gamma-log-to-Xi transfer chain are now
materialized as executable certificates.  The canonical command is

```text
python tools/rh_h1314_h1315_end_to_end_certificate.py
```

and its canonical report is

```text
research/riemann/h1314_h1315_end_to_end_certificate.json.
```

Two complete runs give the same SHA-256 digest:

```text
5169DA3CB27D2444112A1B26E334AD129533680DCE169713B78421831887904F
```

The final report has `all_checks_pass=true`.

## Regenerated Components

The gate executes every component twice, requires the expected certificate
identifier, checks `all_checks_pass=true`, and compares the two JSON hashes.

```text
H921  main saddle coefficient                       pass
H941  exact low-r cover: 23 artifacts, 724 boxes   pass
H943C five-interval U4 rational certificate         pass
H1314 parameter polynomial certificate              pass
H1314 pairing and tail analytic bridge              pass
H1314 uniform pair cap                              pass
H1315 global moment-box assembly                    pass
```

### H921

The standalone certificate reconstructs the saddle identities and proves

```text
0 < B/A^3 <= 7/(q^2 r),       alpha=9/4, r>=2.
```

### H941

The low-band assembler validates exactly 15 selected H932 artifacts and 8
selected H937 artifacts, for a total of 724 certified boxes.  It checks all
internal seams and proves exact gap-free coverage

```text
[2,5/2].
```

### H1314

The parameter engine reconstructs all Touchard/gamma-log polynomials.  It
checks four full Bernstein tensors on `5/2<=r<=3`, four shifted positive
polynomials on `r>=3`, and the tilted growth constant `21/8`.

The separate U4 certificate regenerates the five exact maxima

```text
1687/233280,
431/233280,
-17/29952,
-489/14080,
959/112500,
```

all strictly below `1/100`.

The analytic bridge records and checks the pair identity, the even-expectation
reduction, the Jensen lower bound for `O`, the `K_minus` and `K_plus` bounds,
and the exact tail constants

```text
H(1) <= -17/4,
H'(z) <= -(11/4) exp(z).
```

Together these prove the executable uniform cap

```text
R_r(t) <= 51/100,       r>=5/2, t>=0.
```

### H1315

The gate consumes the new H941 and H1314 JSON results directly, then compares
the H1314 budgets with H1315:

```text
663/1250   < 1,
14229/5000 < 3.
```

It also consumes the H921 coefficient `7` and verifies the exact H922 factor

```text
5*7=35,
```

giving

```text
K_alpha'''(q) >= -35/(q^2 r),       alpha=9/4, r>=2.
```

## Explicit Link To H1317

The gate executes H1317 twice and checks its expected ID, deterministic
output, full Xi bound and correction constant.  In particular it verifies the
previously implicit change of variables:

```text
q=2t,
d_t^3=8 d_q^3,
8*35/(2^2)=70.
```

Thus the H1317 bound consumes exactly the H1315 constant:

```text
kappa_Xi'''(t)
 >= -(70/r(t) + 45450578214/t)/t^2.
```

## Canonical And Preliminary Entrypoints

Use these canonical files:

```text
tools/rh_h1314_h943c_parameter_polynomial_verifier.py
tools/rh_h1314_h943c_uniform_pair_verifier.py
tools/rh_h1315_gamma_log_moment_box_global_certificate.py
tools/rh_h1314_h1315_end_to_end_certificate.py
```

The sibling H1314 file ending in `_parameter_polynomial_certificate.py` is the
exact symbolic engine; its direct serializer is noncanonical, so the wrapper
ending in `_verifier.py` must be used.  The H1315 file ending in
`_global_assembly.py` is a preserved preliminary run whose text-fragment check
failed; the canonical H1315 entrypoint ends in `_global_certificate.py`.
`rh_h1314_h1315_materialized_chain.py` is the earlier passing integration
gate and is superseded by the stricter end-to-end certificate above.

## Remaining Reproducibility Boundary

The separate, stronger H946 branch

```text
R_r(t) <= 1/2 + C_r t^2
```

still lacks its original complete shifted polynomial because that polynomial
was not stored in the AXMEM cards.  It was not reconstructed by invention.
This branch is redundant for H1315 and for the H1317 Xi transfer, which use
the fully materialized uniform `51/100` route.

This work improves reproducibility.  It does not close the enormous finite
H920 interval, the all-degree Jensen/total-positivity step, or the Riemann
Hypothesis.
