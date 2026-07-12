# H816 D3 Cone Sparse-Row Coupling

Classification: `h816_d3_cone_direct_certificate_no_hit`

## Candidate Certificate

P(q) = Q(q) + sum c_i M_i(q) D3_{s_i}(q), where Q has nonnegative coefficients, c_i>0, M_i are monomials, and D3_s are translated contiguous order-3 Toeplitz minors normalized by positive factors.

D3 shape:

`D3_k = 1 - 2 q_k + q_k^2 q_{k+1} + q_{k-1} q_k^2 - q_{k-1} q_k^2 q_{k+1} for k>=3, with the k=2 endpoint D3_2=1-2q2+q2^2 q3.`

## Counts

```json
{
  "h815_remaining_input": 356,
  "covered_by_direct_d3_cone": 0,
  "uncovered": 356
}
```

## Covered Families


## Covered Examples

## Uncovered Families

- `row_gaps=(1,2), col_gaps=(1,2)`: 10
- `row_gaps=(1,2), col_gaps=(2,1)`: 10
- `row_gaps=(1,2), col_gaps=(1,1)`: 10
- `row_gaps=(2,1), col_gaps=(2,1)`: 10
- `row_gaps=(2,1), col_gaps=(1,1)`: 10
- `row_gaps=(2,2), col_gaps=(1,1)`: 10
- `row_gaps=(1,2), col_gaps=(3,2)`: 6
- `row_gaps=(1,2), col_gaps=(1,3)`: 6
- `row_gaps=(1,2), col_gaps=(2,2)`: 6
- `row_gaps=(1,2), col_gaps=(3,1)`: 6
- `row_gaps=(1,3), col_gaps=(3,1)`: 6
- `row_gaps=(1,3), col_gaps=(1,2)`: 6
- `row_gaps=(1,3), col_gaps=(2,1)`: 6
- `row_gaps=(1,3), col_gaps=(1,1)`: 6
- `row_gaps=(2,1), col_gaps=(2,2)`: 6
- `row_gaps=(2,1), col_gaps=(3,1)`: 6
- `row_gaps=(2,1), col_gaps=(1,2)`: 6
- `row_gaps=(2,2), col_gaps=(1,2)`: 6
- `row_gaps=(2,2), col_gaps=(2,1)`: 6
- `row_gaps=(2,3), col_gaps=(2,1)`: 6

## Decision

The direct monomial-D3 cone did not certify any H815 remainder form. This argues against another local coefficient-greedy refinement; move to Plucker-Dodgson identities or the Xi-kernel analytic pivot.

## Next Target

`H817 Plucker-Dodgson sparse-row identity`: Search for determinant identities expressing remaining sparse-row minors through H811-covered consecutive-row minors, PF2 minors, and D3 inputs. If direct identities are circular, pivot to the H804 Xi-kernel unit-deficit route.
