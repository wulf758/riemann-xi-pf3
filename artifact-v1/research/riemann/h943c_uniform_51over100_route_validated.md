# H943C Uniform 51/100 Route Validated

Classification: `h943c_uniform_51over100_route_validated`

## Verdict

Second adversarial audit: the route is valid and the uniform bound

```text
R_r(t) <= 51/100
```

effectively closes H943-C. No sign, normalization factor, or tail problem was
found.

## Normalizations Checked Exact

```text
E_nu[X] = 1,
lambda_M = a X,
z = a t,
```

with `nu_m = c_m m^2/A`, `mu = B/A`, `X = M/mu`, `s = A^3/B^2`. The announced
formulas for `K_e`, `K_o`, `O` and

```text
R = K_e - K_o tanh(O)
```

are correct.

## Regime `0 <= z <= 1`: Simpler Certificates Suffice

The tight constants `delta2`, `delta3` and `E[X^5] <= 16/5` are unnecessary.
The simpler certified package

```text
s >= 18,
E[X^3] <= 61/40,
E[X^5] <= 4,
E[X^5 exp(zX)] <= E[X^5] exp(21z/8)   for 0 <= z <= 1,
```

yields

```text
R - 1/2 <= U4(z) = 61 z^2/960 + z^4 exp(21z/8)/180 - (z/6) tanh(3z^3) < 1/100.
```

(Note `(16/5)/720 = 1/225` vs `4/720 = 1/180`: the looser `E[X^5] <= 4`
produces the `1/180` coefficient.)

The last inequality has a rational five-interval proof, recorded in
`h943c_u4_five_interval_proof.md`, with successive rational maxima

```text
1687/233280,  431/233280,  -17/29952,  -489/14080,  959/112500,
```

all `< 1/100`.

## Regime `z >= 1`

```text
R - 1/2 <= (1/2) exp(H(z)),
H(z) = r(exp(z/mu) - 1) + 3z/mu - 2s(sinh z - z).
```

`H'(z) < 0` and `H(1) < -4`, hence

```text
R - 1/2 < (1/2) exp(-4) < 1/100.
```

## Coverage

No tail or normalization gap: the pairing covers all `t >= 0`, `z = a t`
covers all `z >= 0`, and the normalizer cancels exactly.

## Wide Numerical Search

```text
2.5 <= r <= 500, 0.001 <= z <= 20:
max observed R = 0.5016071541428996 near r = 2.5, z = 0.23.
No counterexample.
```

## Budget Consequences (with H1313)

```text
(51/100)(26/25)       = 0.5304 < 1,
(51/100)(7/2 + 52/25) = 2.8458 < 3.
```

Both H943-C budgets are closed.

## Provenance

Materialized on 2026-07-09 from AXMEM card
`h943c-uniform-51over100-route-validated`. Numbers transcribed from the card;
the standalone verifier remains to be persisted.
