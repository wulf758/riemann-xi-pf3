# Addendum H1321/H1322 To The H1319 Piecewise Closure

Classification: `h1319_finite_dependency_corrected_by_h1321_h1322`

The original H1319 closure correctly certified the analytic block `n>=127`,
but its finite block `2<=n<=126` cited H804/H801.  H801 explicitly recorded
that its Cauchy/DFT extraction had no rigorous outer-circle alias bound and
that cross-radius agreement was only a stability test.  Therefore the old
finite-block description was stronger than its dependency justified.

H1321 now supplies the missing analytic alias remainder and recomputes all 125
finite H803 inequalities with enlarged Arb coefficient balls.  H1322 is the
strict all-index assembly:

```text
2 <= n <= 126: H1321 alias-repaired direct Arb certificate;
n >= 127:      H1319/H920 full-Xi analytic certificate.
```

Canonical follow-up artifacts:

```text
research/riemann/h1321_h801_alias_repaired_finite_h803.json
research/riemann/h1322_h803_all_n_strict_assembly.json
```

This correction changes neither the H1319 analytic tail proof nor its warning
that H907, PF-infinity, all-degree Jensen hyperbolicity, and RH remain open.

