# H815 Four-Term Sparse Bridge

Classification: `h815_four_term_sparse_bridge_partial_with_failures`

## Lemma

For 0<q2<=q3<=...<=1 with q2<=1/2, certify normalized determinant polynomials P=1-A-B+Z by either: (i) B contains q2, Z=B*W, and W>=A^beta for some beta<=2; or (ii) the stronger subcase B=A*q2*U, U<=1, Z=B*W, and W^2>=A^3.

Scalar reduction: `case (i): P>=1/2-A+(1/2)A^beta>=1/2-A+(1/2)A^2>=0; case (ii): P>=1-A-(A/2)*(1-A^(3/2))>=0.`

## Counts

```json
{
  "h813_mixed_unclassified_input": 441,
  "covered_by_four_term_bridge": 85,
  "four_term_failures": 200,
  "not_four_term_remaining": 156
}
```

## Covered Families

- `row_gaps=(1,2), col_gaps=(1,1)`: 5
- `row_gaps=(1,2), col_gaps=(3,1)`: 4
- `row_gaps=(1,3), col_gaps=(2,1)`: 4
- `row_gaps=(1,3), col_gaps=(1,1)`: 4
- `row_gaps=(1,2), col_gaps=(4,1)`: 3
- `row_gaps=(1,3), col_gaps=(4,1)`: 3
- `row_gaps=(1,3), col_gaps=(2,2)`: 3
- `row_gaps=(1,4), col_gaps=(3,1)`: 3
- `row_gaps=(1,4), col_gaps=(2,1)`: 3
- `row_gaps=(1,4), col_gaps=(1,1)`: 3
- `row_gaps=(2,2), col_gaps=(3,1)`: 3
- `row_gaps=(1,2), col_gaps=(5,1)`: 2
- `row_gaps=(1,3), col_gaps=(5,1)`: 2
- `row_gaps=(1,4), col_gaps=(5,1)`: 2
- `row_gaps=(1,4), col_gaps=(3,2)`: 2
- `row_gaps=(1,4), col_gaps=(2,2)`: 2
- `row_gaps=(1,5), col_gaps=(4,1)`: 2
- `row_gaps=(1,5), col_gaps=(3,1)`: 2
- `row_gaps=(1,5), col_gaps=(2,1)`: 2
- `row_gaps=(1,5), col_gaps=(1,1)`: 2

## Covered Examples

### rows [0, 1, 3] cols [1, 4, 5]

Polynomial: `1 - q3*q4 - q2*q3*q4 + q2*q3^2*q4^2*q5`

Certificate:

```json
{
  "type": "half_compensated_four_term_bridge",
  "A": "q3*q4",
  "B": "q2*q3*q4",
  "Z": "q2*q3^2*q4^2*q5",
  "W_Z_over_B": "q3*q4*q5",
  "beta": "2/1",
  "proof": "P=1-A-B+Z=1-A-B*(1-W). Since B contains q2, B<=1/2. Since W>=A^beta with beta<=2 under monotone q_i, P>=1-A-(1/2)*(1-A^beta)>=1/2-A+(1/2)*A^2>=0."
}
```

### rows [0, 1, 3] cols [1, 5, 6]

Polynomial: `1 - q4*q5 - q2*q3*q4*q5 + q2*q3*q4^2*q5^2*q6`

Certificate:

```json
{
  "type": "half_compensated_four_term_bridge",
  "A": "q4*q5",
  "B": "q2*q3*q4*q5",
  "Z": "q2*q3*q4^2*q5^2*q6",
  "W_Z_over_B": "q4*q5*q6",
  "beta": "2/1",
  "proof": "P=1-A-B+Z=1-A-B*(1-W). Since B contains q2, B<=1/2. Since W>=A^beta with beta<=2 under monotone q_i, P>=1-A-(1/2)*(1-A^beta)>=1/2-A+(1/2)*A^2>=0."
}
```

### rows [0, 1, 3] cols [1, 6, 7]

Polynomial: `1 - q5*q6 - q2*q3*q4*q5*q6 + q2*q3*q4*q5^2*q6^2*q7`

Certificate:

```json
{
  "type": "half_compensated_four_term_bridge",
  "A": "q5*q6",
  "B": "q2*q3*q4*q5*q6",
  "Z": "q2*q3*q4*q5^2*q6^2*q7",
  "W_Z_over_B": "q5*q6*q7",
  "beta": "2/1",
  "proof": "P=1-A-B+Z=1-A-B*(1-W). Since B contains q2, B<=1/2. Since W>=A^beta with beta<=2 under monotone q_i, P>=1-A-(1/2)*(1-A^beta)>=1/2-A+(1/2)*A^2>=0."
}
```

### rows [0, 1, 3] cols [1, 7, 8]

Polynomial: `1 - q6*q7 - q2*q3*q4*q5*q6*q7 + q2*q3*q4*q5*q6^2*q7^2*q8`

Certificate:

```json
{
  "type": "half_compensated_four_term_bridge",
  "A": "q6*q7",
  "B": "q2*q3*q4*q5*q6*q7",
  "Z": "q2*q3*q4*q5*q6^2*q7^2*q8",
  "W_Z_over_B": "q6*q7*q8",
  "beta": "2/1",
  "proof": "P=1-A-B+Z=1-A-B*(1-W). Since B contains q2, B<=1/2. Since W>=A^beta with beta<=2 under monotone q_i, P>=1-A-(1/2)*(1-A^beta)>=1/2-A+(1/2)*A^2>=0."
}
```

### rows [0, 1, 3] cols [2, 3, 4]

Polynomial: `1 - q3 - q2*q3 + q2*q3^2*q4`

Certificate:

```json
{
  "type": "half_compensated_four_term_bridge",
  "A": "q3",
  "B": "q2*q3",
  "Z": "q2*q3^2*q4",
  "W_Z_over_B": "q3*q4",
  "beta": "2/1",
  "proof": "P=1-A-B+Z=1-A-B*(1-W). Since B contains q2, B<=1/2. Since W>=A^beta with beta<=2 under monotone q_i, P>=1-A-(1/2)*(1-A^beta)>=1/2-A+(1/2)*A^2>=0."
}
```

### rows [0, 1, 4] cols [1, 5, 6]

Polynomial: `1 - q3*q4*q5 - q2*q3*q4*q5 + q2*q3^2*q4^2*q5^2*q6`

Certificate:

```json
{
  "type": "half_compensated_four_term_bridge",
  "A": "q3*q4*q5",
  "B": "q2*q3*q4*q5",
  "Z": "q2*q3^2*q4^2*q5^2*q6",
  "W_Z_over_B": "q3*q4*q5*q6",
  "beta": "2/1",
  "proof": "P=1-A-B+Z=1-A-B*(1-W). Since B contains q2, B<=1/2. Since W>=A^beta with beta<=2 under monotone q_i, P>=1-A-(1/2)*(1-A^beta)>=1/2-A+(1/2)*A^2>=0."
}
```

### rows [0, 1, 4] cols [1, 6, 7]

Polynomial: `1 - q4*q5*q6 - q2*q3*q4*q5*q6 + q2*q3*q4^2*q5^2*q6^2*q7`

Certificate:

```json
{
  "type": "half_compensated_four_term_bridge",
  "A": "q4*q5*q6",
  "B": "q2*q3*q4*q5*q6",
  "Z": "q2*q3*q4^2*q5^2*q6^2*q7",
  "W_Z_over_B": "q4*q5*q6*q7",
  "beta": "2/1",
  "proof": "P=1-A-B+Z=1-A-B*(1-W). Since B contains q2, B<=1/2. Since W>=A^beta with beta<=2 under monotone q_i, P>=1-A-(1/2)*(1-A^beta)>=1/2-A+(1/2)*A^2>=0."
}
```

### rows [0, 1, 4] cols [1, 7, 8]

Polynomial: `1 - q5*q6*q7 - q2*q3*q4*q5*q6*q7 + q2*q3*q4*q5^2*q6^2*q7^2*q8`

Certificate:

```json
{
  "type": "half_compensated_four_term_bridge",
  "A": "q5*q6*q7",
  "B": "q2*q3*q4*q5*q6*q7",
  "Z": "q2*q3*q4*q5^2*q6^2*q7^2*q8",
  "W_Z_over_B": "q5*q6*q7*q8",
  "beta": "2/1",
  "proof": "P=1-A-B+Z=1-A-B*(1-W). Since B contains q2, B<=1/2. Since W>=A^beta with beta<=2 under monotone q_i, P>=1-A-(1/2)*(1-A^beta)>=1/2-A+(1/2)*A^2>=0."
}
```

## Remaining Families

- `row_gaps=(1,2), col_gaps=(1,2)`: 10
- `row_gaps=(2,1), col_gaps=(2,1)`: 10
- `row_gaps=(2,2), col_gaps=(1,1)`: 10
- `row_gaps=(1,2), col_gaps=(2,1)`: 10
- `row_gaps=(1,2), col_gaps=(1,1)`: 10
- `row_gaps=(2,1), col_gaps=(1,1)`: 10
- `row_gaps=(1,2), col_gaps=(3,2)`: 6
- `row_gaps=(1,2), col_gaps=(1,3)`: 6
- `row_gaps=(1,2), col_gaps=(2,2)`: 6
- `row_gaps=(1,2), col_gaps=(3,1)`: 6
- `row_gaps=(1,3), col_gaps=(1,2)`: 6
- `row_gaps=(1,3), col_gaps=(2,1)`: 6
- `row_gaps=(2,1), col_gaps=(2,2)`: 6
- `row_gaps=(2,1), col_gaps=(3,1)`: 6
- `row_gaps=(2,2), col_gaps=(1,2)`: 6
- `row_gaps=(2,2), col_gaps=(2,1)`: 6
- `row_gaps=(2,3), col_gaps=(2,1)`: 6
- `row_gaps=(2,3), col_gaps=(1,1)`: 6
- `row_gaps=(3,1), col_gaps=(2,1)`: 6
- `row_gaps=(3,2), col_gaps=(1,1)`: 6

## Decision

The four-term bridge proves many H813 mixed forms but not every four-term case. Inspect the failed four-term cases before promoting this to a lemma.

## Next Target

`H816 six-term boundary/D3 coupling`: Use the translated contiguous D3 input to certify the remaining higher-term mixed sparse-row polynomials, or pivot to a Plucker-Dodgson identity if direct product bounds stall.
