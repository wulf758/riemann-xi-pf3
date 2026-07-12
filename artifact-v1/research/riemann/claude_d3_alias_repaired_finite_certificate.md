# Alias-Repaired Finite Translated-D3 Certificate

Classification: `xi_translated_d3_alias_repaired_finite_certified`

## Result

The alias-repaired Xi coefficient balls rigorously certify translated D3 for every integer 2<=n<=126. This closes the finite seam only; an analytic or independently certified tail for n>=127 is still required.

The certified inequality is

```text
D_n = (1-q_(n+1))^2 - q_(n+1)^2 (1-q_n)(1-q_(n+2)) > 0,
q_n = a_n a_(n-2) / a_(n-1)^2.
```

- certified integer range: `[2, 126]`
- row count: `125`
- worst lower bound: `[3.500735924525638225876966076684099048181659519265812448360346059928989252545840256423996567025827292492290024648440111684301152785307660224377475619393780434610820958012364654632032464815594949524824036525634727694201652133299447322034734715906799298707420053174178624492044682731e-6 +/- 1.63e-286]` at `n=126`
- determinant cross-evaluation worst lower bound: `3.500735924525638225875386395948707667480194697116556080982812547400706733313757723635717411525547504425048828125000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e-6`

## Proof Contract

The H1321 root-of-unity coefficient balls are recomputed. Each is enlarged by

```text
Xi(iR) R^(-2n) (r/R)^N / (1-(r/R)^N),  r=11, R=12, N=8192,
```

before any quotient or determinant is formed. Every row is then evaluated both
from the quotient formula and from the contiguous 3x3 Toeplitz determinant
divided by `a_n^3`.

## Checks

- `python_flint_loaded`: `True`
- `reference_radius_is_11`: `True`
- `outer_radius_is_12`: `True`
- `sample_count_is_8192`: `True`
- `n_max_reaches_128`: `True`
- `outer_xi_real_part_positive`: `True`
- `outer_xi_imaginary_part_contains_zero`: `True`
- `alias_ratio_strictly_between_zero_and_one`: `True`
- `all_alias_bounds_positive`: `True`
- `all_repaired_coefficients_positive`: `True`
- `pf2_ratios_strictly_decreasing_through_R128`: `True`
- `q2_strictly_below_one_half`: `True`
- `q_strictly_increasing_through_q128`: `True`
- `row_indices_exactly_2_through_126`: `True`
- `all_d3_q_lower_bounds_positive`: `True`
- `all_d3_determinant_lower_bounds_positive`: `True`
- `all_cross_evaluation_differences_contain_zero`: `True`
- `coefficient_mutation_canary_breaks_certificate`: `True`
- `sign_canary_breaks_certificate`: `True`

## Canaries

- coefficient mutation killed the target row: `True`
- sign/comparator mutation killed the target row: `True`

## Scope

- Does not prove translated D3 for n>=127.
- Does not prove full PF3 (only H811 consecutive-row consequences are relevant).
- Does not prove PF-infinity.
- Does not prove all-degree Jensen hyperbolicity.
- Does not prove the Riemann Hypothesis.
