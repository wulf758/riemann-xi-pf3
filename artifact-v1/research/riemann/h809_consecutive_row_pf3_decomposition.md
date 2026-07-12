# H809 Consecutive-Row PF3 Decomposition

Classification: `consecutive_row_identity_decomposition_verified_but_pf2_d3_proof_invalidated`

## Statement

For one-sided Toeplitz order-3 minors with consecutive rows, the minors split into exact ratio/geometry identities. The earlier claim that PF2 plus translated contiguous D3 nonnegativity implies all sparse-column minors is false; extra anchored minors are needed.

## Families

| family | columns | formula | source of nonnegativity |
| --- | --- | --- | --- |
| `first_col_0` | `(0,j,k), 2<=j<k` | `det = a0*a_{j-2}*a_{k-2}*(R_{j-1}-R_{k-1})` | PF2 gives R_{j-1}>=R_{k-1}. |
| `first_col_1` | `(1,j,k), 2<=j<k` | `det = a0*a_{j-2}*a_{k-2}*(G_j-G_k), G_c=R_{c-1}(R1-R_c)` | Requires anchored fixed-left minors with cols (1,c,c+1); translated contiguous D3 alone is not enough. |
| `first_col_ge2` | `(i,j,k), 2<=i<j<k` | `det = a_{i-2}*a_{j-2}*a_{k-2}*orient((X_i,Y_i),(X_j,Y_j),(X_k,Y_k))` | X_c=R_{c-1} is monotone by PF2; contiguous D3 gives local convexity, hence all triple orientations. |

## Identity Checks

```json
{
  "trials": 100,
  "max_col": 10,
  "failure_count": 0,
  "failures": []
}
```

## Invalidating Counterexample

```json
{
  "source": "H796",
  "sequence": [
    1,
    7,
    42,
    252,
    1512,
    7560
  ],
  "ratios": [
    7,
    6,
    6,
    6,
    5
  ],
  "translated_contiguous_d3": {
    "D0": "1",
    "D1": "7",
    "D2": "0",
    "D3": "0"
  },
  "negative_sparse_consecutive_row_minor": {
    "rows": [
      0,
      1,
      2
    ],
    "cols": [
      1,
      2,
      5
    ],
    "determinant": "-1260"
  },
  "lesson": "G-monotonicity is controlled by anchored minors (1,c,c+1), not by translated contiguous D3."
}
```

## Shape Coverage

```json
{
  "max_index": 8,
  "counts": {
    "consecutive_rows_total": 588,
    "consecutive_rows_covered": 588,
    "consecutive_rows_zero_by_negative_shift": 378,
    "consecutive_rows_formula_covered": 210,
    "other_rows_uncovered": 6468,
    "total_order3_minors": 7056
  },
  "examples": {
    "covered_formula": {
      "rows": [
        0,
        1,
        2
      ],
      "cols": [
        0,
        1,
        2
      ]
    },
    "uncovered_other_rows": {
      "rows": [
        0,
        1,
        3
      ],
      "cols": [
        0,
        1,
        2
      ]
    },
    "covered_zero_shift": {
      "rows": [
        1,
        2,
        3
      ],
      "cols": [
        0,
        1,
        2
      ]
    }
  }
}
```

## Decision

All tested exact identities passed, but the proposed proof from PF2 plus translated contiguous D3 is invalidated by H796. The first-column-1 family needs monotonicity of G_c=R_{c-1}(R1-R_c); adjacent differences of G_c are fixed-left minors with columns (1,c,c+1), not the translated contiguous D3 family (c-1,c,c+1). H809 is therefore an identity decomposition and a map of the missing anchored-minor hierarchy, not a proof slice.

## Limitations

- This is an identity decomposition, not a proof from the H807 assumptions.
- PF2 plus translated contiguous D3 is insufficient, by H796.
- A viable route needs an anchored-minor hierarchy such as fixed-left minors (1,c,c+1), or a different Xi-specific constraint.

## Next Target

`H810 anchored-minor hierarchy audit`: Try to transform sparse-row/consecutive-column and fully sparse order-3 Toeplitz minors into a dual convexity or monotone-slope problem.
