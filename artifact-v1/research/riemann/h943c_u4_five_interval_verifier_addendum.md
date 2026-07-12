# H943-C U4 Five-Interval Verifier Addendum

Classification: `h943c_u4_five_interval_proof_verified`

The proof recorded in `h943c_u4_five_interval_proof.md` now has a standalone
exact verifier:

```text
python tools/rh_h943c_u4_five_interval_certificate.py
```

The script uses `Fraction` arithmetic throughout.  It certifies the five
exponential majorants by Taylor polynomials with geometric remainder bounds,
the hyperbolic coefficient lower bounds

```text
4/9, 4/13, 3/11, 5/32,
```

the monotonicity of each quadratic majorant, and the exact maxima

```text
1687/233280,
431/233280,
-17/29952,
-489/14080,
959/112500.
```

Every maximum is strictly below `1/100`.  The canonical generated report is

```text
research/riemann/h943c_u4_five_interval_certificate.json
```

and reports `all_checks_pass=true`.

The original proof note remains unchanged because the Windows sandbox blocked
an in-place `apply_patch`; this addendum is the canonical persistence link.
