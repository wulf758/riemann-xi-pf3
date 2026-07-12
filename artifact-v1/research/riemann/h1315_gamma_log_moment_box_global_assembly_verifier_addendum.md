# H1315 Verifier Addendum

Classification: `h1315_gamma_log_moment_box_global_assembly_verifier_path_correction`

The canonical H1315 verifier is:

```text
tools/rh_h1315_gamma_log_moment_box_global_certificate.py
```

It writes:

```text
research/riemann/h1315_gamma_log_moment_box_global_assembly.json
```

The earlier file

```text
tools/rh_h1315_gamma_log_moment_box_global_assembly.py
```

is a preliminary non-canonical draft. Its first run correctly returned
`all_checks_pass=false` because one H941 attestation fragment crossed a
Markdown line break; the mathematics and exact budget arithmetic in that run
passed. A Windows sandbox-wrapper failure prevented in-place correction, so a
new complete verifier was materialized with `apply_patch` instead of mutating
the draft through a different write mechanism.

The canonical verifier also avoids binary-float display artifacts in the JSON
by converting every exact `Fraction` through `Decimal`.

This addendum changes no H1314 artifact and does not upgrade H1314's status.
