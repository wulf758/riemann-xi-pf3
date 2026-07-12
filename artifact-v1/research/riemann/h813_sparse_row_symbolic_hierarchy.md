# H813 Sparse-Row Symbolic Hierarchy

Classification: `h813_sparse_row_symbolic_hierarchy_partial`

## Method

- Represent a_k/a_0 after removing the common positive R1^k factor as a monomial in q_j=R_j/R_{j-1}.
- Expand each order-3 Toeplitz determinant exactly as an integer polynomial in q_j.
- Remove the common positive monomial and integer gcd; sign is unchanged.
- Search simple certificates: nonnegative coefficients, factors 1-q_i with nonnegative quotient, adjacent equality vanishing, and exact D3 polynomial matches.

## Summary

```json
{
  "active_sparse_row_pairs": 2310,
  "zero_polynomials": 0,
  "nonnegative_coefficients": 896,
  "q_box_lower_bound": 973,
  "one_minus_factor_with_nonnegative_quotient": 0,
  "matches_translated_contiguous_d3": 0,
  "mixed_unclassified": 441
}
```

## H812 Boundary Families

### rows [0, 1, 3] cols [3, 4, 8]

Family: `row_gaps=(1,2), col_gaps=(1,4)`

Polynomial: `1 - q4 - q3*q4^2*q5^2*q6^2*q7 + q3*q4^3*q5^3*q6^3*q7^2*q8 + q2*q3^2*q4^3*q5^2*q6^2*q7 - q2*q3^2*q4^3*q5^3*q6^3*q7^2*q8`

Certificate:

```json
{
  "type": "mixed",
  "one_minus_factor_candidates": [],
  "adjacent_monotone_equality_vanishes": []
}
```

### rows [0, 1, 3] cols [3, 5, 7]

Family: `row_gaps=(1,2), col_gaps=(2,2)`

Polynomial: `1 - q4*q5 - q4*q5^2*q6 + q4^2*q5^3*q6^2*q7 + q2*q3^2*q4^3*q5^3*q6 - q2*q3^2*q4^3*q5^3*q6^2*q7`

Certificate:

```json
{
  "type": "mixed",
  "one_minus_factor_candidates": [],
  "adjacent_monotone_equality_vanishes": []
}
```

### rows [0, 1, 4] cols [4, 5, 7]

Family: `row_gaps=(1,3), col_gaps=(1,2)`

Polynomial: `1 - q5 - q3*q4^2*q5^2*q6 + q3*q4^2*q5^3*q6^2*q7 + q2*q3^2*q4^3*q5^3*q6 - q2*q3^2*q4^3*q5^3*q6^2*q7`

Certificate:

```json
{
  "type": "mixed",
  "one_minus_factor_candidates": [],
  "adjacent_monotone_equality_vanishes": []
}
```

### rows [0, 1, 5] cols [4, 5, 7]

Family: `row_gaps=(1,4), col_gaps=(1,2)`

Polynomial: `1 - q5 - q2*q3^2*q4^2*q5^2*q6 + q2*q3^2*q4^2*q5^3*q6^2*q7`

Certificate:

```json
{
  "type": "mixed",
  "one_minus_factor_candidates": [],
  "adjacent_monotone_equality_vanishes": []
}
```

### rows [0, 1, 6] cols [4, 6, 7]

Family: `row_gaps=(1,5), col_gaps=(2,1)`

Polynomial: `1 - q5*q6 - q2*q3*q4*q5*q6 + q2*q3*q4*q5^2*q6^2*q7`

Certificate:

```json
{
  "type": "mixed",
  "one_minus_factor_candidates": [],
  "adjacent_monotone_equality_vanishes": []
}
```

### rows [0, 2, 3] cols [3, 5, 8]

Family: `row_gaps=(2,1), col_gaps=(2,3)`

Polynomial: `1 - q4*q5*q6 - q3*q4^2*q5 + q3*q4^3*q5^3*q6^3*q7^2*q8 + q2*q3^2*q4^3*q5^2*q6 - q2*q3^2*q4^3*q5^3*q6^3*q7^2*q8`

Certificate:

```json
{
  "type": "mixed",
  "one_minus_factor_candidates": [],
  "adjacent_monotone_equality_vanishes": []
}
```

### rows [0, 2, 3] cols [3, 6, 7]

Family: `row_gaps=(2,1), col_gaps=(3,1)`

Polynomial: `1 - q5 - q3*q4^2*q5^2*q6 + q3*q4^2*q5^3*q6^2*q7 + q2*q3^2*q4^3*q5^3*q6 - q2*q3^2*q4^3*q5^3*q6^2*q7`

Certificate:

```json
{
  "type": "mixed",
  "one_minus_factor_candidates": [],
  "adjacent_monotone_equality_vanishes": []
}
```

### rows [0, 2, 4] cols [4, 5, 8]

Family: `row_gaps=(2,2), col_gaps=(1,3)`

Polynomial: `1 - q4*q5 - q3*q4^2*q5^2*q6 + q3*q4^3*q5^4*q6^3*q7^2*q8 + q2*q3^2*q4^3*q5^3*q6 - q2*q3^2*q4^3*q5^4*q6^3*q7^2*q8`

Certificate:

```json
{
  "type": "mixed",
  "one_minus_factor_candidates": [],
  "adjacent_monotone_equality_vanishes": []
}
```

### rows [0, 2, 4] cols [4, 6, 7]

Family: `row_gaps=(2,2), col_gaps=(2,1)`

Polynomial: `1 - q4*q5 - q4*q5^2*q6 + q4^2*q5^3*q6^2*q7 + q2*q3^2*q4^3*q5^3*q6 - q2*q3^2*q4^3*q5^3*q6^2*q7`

Certificate:

```json
{
  "type": "mixed",
  "one_minus_factor_candidates": [],
  "adjacent_monotone_equality_vanishes": []
}
```

### rows [0, 2, 5] cols [5, 6, 7]

Family: `row_gaps=(2,3), col_gaps=(1,1)`

Polynomial: `1 - q5*q6 - q3*q4*q5 + q3*q4*q5^2*q6^2*q7 + q2*q3^2*q4^2*q5^2*q6 - q2*q3^2*q4^2*q5^2*q6^2*q7`

Certificate:

```json
{
  "type": "mixed",
  "one_minus_factor_candidates": [],
  "adjacent_monotone_equality_vanishes": []
}
```

### rows [0, 2, 6] cols [5, 6, 7]

Family: `row_gaps=(2,4), col_gaps=(1,1)`

Polynomial: `1 - q5*q6 - q2*q3*q4*q5 + q2*q3*q4*q5^2*q6^2*q7`

Certificate:

```json
{
  "type": "mixed",
  "one_minus_factor_candidates": [],
  "adjacent_monotone_equality_vanishes": []
}
```

### rows [0, 3, 4] cols [4, 6, 8]

Family: `row_gaps=(3,1), col_gaps=(2,2)`

Polynomial: `1 - q4*q5 - q3*q4^2*q5^2*q6 + q3*q4^3*q5^4*q6^3*q7^2*q8 + q2*q3^2*q4^3*q5^3*q6 - q2*q3^2*q4^3*q5^4*q6^3*q7^2*q8`

Certificate:

```json
{
  "type": "mixed",
  "one_minus_factor_candidates": [],
  "adjacent_monotone_equality_vanishes": []
}
```

## Unclassified Examples

- rows [0, 1, 3] cols [1, 4, 5]: `1 - q3*q4 - q2*q3*q4 + q2*q3^2*q4^2*q5`; certificate `mixed`
- rows [0, 1, 3] cols [1, 4, 6]: `1 - q3*q4^2*q5 - q2*q3*q4 + q2*q3^2*q4^3*q5^2*q6`; certificate `mixed`
- rows [0, 1, 3] cols [1, 4, 7]: `1 - q3*q4^2*q5^2*q6 - q2*q3*q4 + q2*q3^2*q4^3*q5^3*q6^2*q7`; certificate `mixed`
- rows [0, 1, 3] cols [1, 4, 8]: `1 - q3*q4^2*q5^2*q6^2*q7 - q2*q3*q4 + q2*q3^2*q4^3*q5^3*q6^3*q7^2*q8`; certificate `mixed`
- rows [0, 1, 3] cols [1, 5, 6]: `1 - q4*q5 - q2*q3*q4*q5 + q2*q3*q4^2*q5^2*q6`; certificate `mixed`
- rows [0, 1, 3] cols [1, 5, 7]: `1 - q4*q5^2*q6 - q2*q3*q4*q5 + q2*q3*q4^2*q5^3*q6^2*q7`; certificate `mixed`
- rows [0, 1, 3] cols [1, 5, 8]: `1 - q4*q5^2*q6^2*q7 - q2*q3*q4*q5 + q2*q3*q4^2*q5^3*q6^3*q7^2*q8`; certificate `mixed`
- rows [0, 1, 3] cols [1, 6, 7]: `1 - q5*q6 - q2*q3*q4*q5*q6 + q2*q3*q4*q5^2*q6^2*q7`; certificate `mixed`
- rows [0, 1, 3] cols [1, 6, 8]: `1 - q5*q6^2*q7 - q2*q3*q4*q5*q6 + q2*q3*q4*q5^2*q6^3*q7^2*q8`; certificate `mixed`
- rows [0, 1, 3] cols [1, 7, 8]: `1 - q6*q7 - q2*q3*q4*q5*q6*q7 + q2*q3*q4*q5*q6^2*q7^2*q8`; certificate `mixed`
- rows [0, 1, 3] cols [2, 3, 4]: `1 - q3 - q2*q3 + q2*q3^2*q4`; certificate `mixed`
- rows [0, 1, 3] cols [2, 3, 5]: `1 - q3 - q2*q3^2*q4 + q2*q3^3*q4^2*q5`; certificate `mixed`
- rows [0, 1, 3] cols [2, 3, 6]: `1 - q3 - q2*q3^2*q4^2*q5 + q2*q3^3*q4^3*q5^2*q6`; certificate `mixed`
- rows [0, 1, 3] cols [2, 3, 7]: `1 - q3 - q2*q3^2*q4^2*q5^2*q6 + q2*q3^3*q4^3*q5^3*q6^2*q7`; certificate `mixed`
- rows [0, 1, 3] cols [2, 3, 8]: `1 - q3 - q2*q3^2*q4^2*q5^2*q6^2*q7 + q2*q3^3*q4^3*q5^3*q6^3*q7^2*q8`; certificate `mixed`
- rows [0, 1, 3] cols [2, 4, 5]: `1 - 2*q3*q4 + q3^2*q4^2*q5`; certificate `mixed`
- rows [0, 1, 3] cols [2, 4, 6]: `1 - q3*q4 - q3*q4^2*q5 + q3^2*q4^3*q5^2*q6`; certificate `mixed`
- rows [0, 1, 3] cols [2, 4, 7]: `1 - q3*q4 - q3*q4^2*q5^2*q6 + q3^2*q4^3*q5^3*q6^2*q7`; certificate `mixed`
- rows [0, 1, 3] cols [2, 4, 8]: `1 - q3*q4 - q3*q4^2*q5^2*q6^2*q7 + q3^2*q4^3*q5^3*q6^3*q7^2*q8`; certificate `mixed`
- rows [0, 1, 3] cols [2, 5, 6]: `1 - q4*q5 - q3*q4*q5 + q3*q4^2*q5^2*q6`; certificate `mixed`

## Alternate Ideas Kept Alive

- `Xi kernel saddle-point unit-deficit` (H804/H805): The clean target T_{n+1}/T_n >= 1-1/n^2 remains analytic, not grid-mined.
- `Li coefficient interval certification` (H212-H217): Blocked by rigorous outer-circle zeta'/zeta bounds, but could revive with stronger interval tooling.
- `Endpoint/Nyman-Beurling profile lemmas` (H574/H611): Formal candidates survive after several invalidations; they are independent of the PF3 corridor.

## Decision

The H812 sparse-row boundary has simple exact factors, but the active sparse-row universe is not yet reduced to elementary coefficient/factor certificates. The route needs a stronger symbolic identity, likely a Plucker/Dodgson or compound-matrix inequality using H811 minors.

## Next Target

`H814 Plucker-Dodgson sparse-row reduction or analytic pivot`: Either express the mixed unclassified sparse-row polynomials using H811-covered consecutive-row minors through a determinant identity, or pivot to the Xi-kernel unit-deficit lemma instead of more grid search.
