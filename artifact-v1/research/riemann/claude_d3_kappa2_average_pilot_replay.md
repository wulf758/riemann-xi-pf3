# Rigorous Kappa2-Average D3 Pilot

Classification: `xi_d3_kappa2_average_pilot_certified`

## Result

The exact kappa2-average route certifies D_127>0 and overlaps H1321 on q_127,q_128. This is a rigorous pilot row, not a uniform tail proof.

The exact bridge used is

```text
A_n = [Lambda(2n)-2Lambda(2n-2)+Lambda(2n-4)]/4
    = int_[0,1]^2 kappa2(2(n-2)+2(x+y)) dx dy,
q_n = exp(4 A_n) * (2n-2)(2n-3)/((2n)(2n-1)).
```

- target row: `D_127`
- kappa2-average lower bound: `[3.42360862005536875001261459229022687101476980176795312950e-6 +/- 2.06e-63]`
- direct moment-ratio lower bound: `[3.42360862005536874999906541057753358931911822524503539585e-6 +/- 7.66e-64]`
- overlap with H1321 for q127/q128: `True`

## Checks

- `five_even_theta_moments_computed`: `True`
- `all_moment_lower_bounds_positive`: `True`
- `all_q_cross_evaluations_overlap`: `True`
- `d3_kappa2_average_lower_bound_positive`: `True`
- `d3_moment_ratio_lower_bound_positive`: `True`
- `d3_cross_evaluations_overlap`: `True`
- `all_h1321_overlap_checks_pass`: `True`
- `moment_mutation_canary_breaks_certificate`: `True`
- `sign_canary_breaks_certificate`: `True`

## Scope

This is one certified row beyond the H1321 finite seam. It does not
supply a uniform tail lemma for all `n>=127`, PF-infinity, Jensen
hyperbolicity in all degrees, or RH.
