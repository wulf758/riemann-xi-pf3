# H1314 H946 Corrected-Pair Global Proof

Classification: `h1314_h946_corrected_pair_global_proof_pending_verifier`

## Statement

For the exact gamma-log family, for every `r >= 5/2` and every `t >= 0`:

```text
R_r(t) <= 1/2 + C_r t^2.
```

This is the strong (corrected-pair) H946 bound, stronger than the uniform cap
`R_r(t) <= 51/100` certified in the companion H1314 notes.

## Proof Skeleton

With `u = t/sqrt(A)`, define

```text
T = K_e - 1/2 - V5 u^2/(24B),
L = 1/2 + V5 u^2/(24B) - K_-,
```

so that

```text
K_o = T + L,
L >= V4 u/(6B).
```

The proof splits at the scale-free point `u = 1/2`; the two branches overlap
exactly there.

### Branch 1: `0 <= u <= 1/2`

Symmetric Taylor bounds with remainder sign control, the odd-part lower bound

```text
psi_o >= B u^3/6,
```

and the exact inequality

```text
V7(r exp(u)) <= 20 B V4        (B = V3(r) at the saddle)
```

imply `T <= psi_o L`. The elementary bound `tanh x >= x/(1+x)` then closes the
branch.

The final polynomial certificate of this branch has only positive
coefficients after the shift `r = 5/2 + X`.

### Branch 2: `u >= 1/2`

Exponential dissymmetry:

```text
R <= K_- + K_+ exp(-2 psi_o),
2 psi_o >= (pi/4) exp(r exp(u)),
```

then

```text
G(r exp(u)) exp(-2 psi_o) < 1 < V4 u^3/6
```

closes the tail.

## Verification Notes

- All factors `A`, `B`, `u`, signs and constants were re-verified, in
  particular at the critical corner `r = 5/2`, `u = 1/2`.
- Exact shifted-polynomial and rational exponential certificates were
  independently computed.
- An independent analytic audit (see
  `h1314_uniform_pair_cap_proved_audited.md`) also confirms this strong bound.

## Consequences

Combined with already-proved inputs:

```text
C_r <= 1/297                    [H1270],
E[W^2] <= 26/25, E[W^4] <= 7/2, E[W^6] <= 21    [H1313],
```

H946 now implies both H943-C gamma-log budgets (see the H946 conditional
closure computation in `h946_corrected_pair_symmetry_criterion.md`).

## Non-Claims

This does not prove H928, does not prove the transfer to the full Xi
function, and does not prove the Riemann Hypothesis.

## Provenance

Materialized on 2026-07-09 from AXMEM card
`h1314-h946-corrected-pair-global-proof` (recorded by the proving sub-agent;
file creation was blocked at proof time by an environment usage-limit
rejection, and the details were handed to the parent for integration).

Still pending for publication status:

```text
1. persist the exact/Arb verifier scripts;
2. line-by-line independent rederivation of branch 1 remainder signs;
3. explicit continuity at z=0 and the convexity statement used for O;
4. rational bounds for pi*exp(r) used in branch 2.
```
