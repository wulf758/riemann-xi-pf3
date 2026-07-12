# H814 Sparse-Row Lemma Candidates After H813

Classification: `h814_sparse_row_lemma_candidates_extracted`

## Current State

H811 proves the conditional order-3 Toeplitz theorem for consecutive rows under:

- `a_n>0`
- `R_n=a_n/a_{n-1}` nonincreasing
- `q_n=R_n/R_{n-1}` nondecreasing
- `q_2<=1/2`
- translated contiguous order-3 Toeplitz minors nonnegative

H812 found no exact sparse-row counterexample on the full `d=16`, length-7 rational grid under those inputs.

H813 then expanded all active sparse-row order-3 minors in exact `q_n` variables. It certified:

- `896` active sparse-row forms by nonnegative coefficients;
- `973` more by the elementary `q_2` budget box bound;
- `441` active forms still mixed/unclassified.

So the next move is not a denser grid. The remaining problem is a symbolic sparse-row reduction.

## Lemma Candidates

### Lemma A: q2-Budget Monomial Bound

Let `0<q_2<=1/2` and `0<q_i<=1` for all `i>=3`. For any finite polynomial

`P(q)=c + P_+(q) - sum_j b_j M_j(q)`,

where `c>=0`, `P_+` has nonnegative coefficients, `b_j>0`, and `M_j` are monomials, a sufficient condition for `P(q)>=0` is

`c >= sum_j b_j 2^{-v_2(M_j)}`,

where `v_2(M_j)` is the exponent of `q_2` in `M_j`.

Status: proved as an elementary certificate and implemented in H813. This is useful but not enough.

### Lemma B: Four-Term Sparse Bridge

For monotone ratios

`0<q_2<=q_3<=...<=q_m<=1`, `q_2<=1/2`,

prove the nonnegativity of the sparse-row four-term family

`1 - A(q) - q_2 B(q) + q_2 C(q) >= 0`,

where `A,B,C` are interval-product monomials produced by rows `(0,1,r)` and columns `(1,j,k)`, and `C` is the determinant compensation monomial. The first unresolved examples are:

- `1 - q_3 q_4 - q_2 q_3 q_4 + q_2 q_3^2 q_4^2 q_5`
- `1 - q_4 q_5 - q_2 q_3 q_4 q_5 + q_2 q_3 q_4^2 q_5^2 q_6`
- `1 - q_5 q_6 - q_2 q_3 q_4 q_5 q_6 + q_2 q_3 q_4 q_5^2 q_6^2 q_7`

This is the nearest symbolic target because it accounts for most low-complexity unclassified examples.

### Lemma C: Six-Term Boundary/D3 Coupling

For monotone `q_n` with `q_2<=1/2`, prove the six-term boundary family

`1 - A(q) - B(q) + B'(q) + q_2 C(q) - q_2 C'(q) >= 0`

whenever the corresponding translated contiguous D3 minors are nonnegative.

Representative unresolved forms include:

- `1 - q_4 - q_3 q_4 + q_3 q_4^2 q_5 + q_2 q_3^2 q_4^2 - q_2 q_3^2 q_4^2 q_5`
- `1 - q_4 - q_3 q_4^2 q_5 + q_3 q_4^3 q_5^2 q_6 + q_2 q_3^2 q_4^3 q_5 - q_2 q_3^2 q_4^3 q_5^2 q_6`

This likely needs the translated D3 input; the q-box bound alone cannot see it.

### Lemma D: Plucker-Dodgson Sparse-Row Reduction

Let `T=(a_{j-i})_{i,j>=0}` be the one-sided Toeplitz matrix. Under H811 inputs, show that every order-3 minor with nonconsecutive rows can be expressed through a nonnegative Plucker/Dodgson relation involving:

- H811-covered consecutive-row order-3 minors;
- order-2 Toeplitz minors from PF2;
- translated contiguous order-3 minors.

This is the cleanest route to a proof of PF3 from H811 plus a finite symbolic identity family.

### Lemma E: Xi-Kernel Unit-Deficit Pivot

If Lemmas B-D stall, pivot away from sparse-row algebra and return to the Xi-specific analytic target from H804:

`T_{n+1}/T_n >= 1 - 1/n^2`

for the Xi kernel moment-curvature quotient, eventually with a finite base check.

This route is independent of the sparse-row determinant algebra and remains underexploited.

## Decision

H813 did not close PF3, but it changed the shape of the problem: the hard part is now a small family of mixed sparse-row polynomials, not arbitrary minor mining.

Next work should attack Lemma B first, then Lemma C or D. If that becomes circular or RH-equivalent, pivot to Lemma E rather than increasing the grid.
