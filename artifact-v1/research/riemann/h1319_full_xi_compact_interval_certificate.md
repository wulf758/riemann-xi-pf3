# H1319 Full-Xi Compact Interval Certificate

Classification: `h1319_full_xi_compact_interval_closed`

## Result

The five certified pieces form a closed, gap-free chain

`r(125) -> 4 -> 5 -> 10 -> 22 -> infinity`.

Since the saddle map is strictly increasing, this proves

`t^2*kappa_Xi'''(t) >= -3/4` for every `125 <= t <= T_71`.

All manifest checks pass: `True`.

## Components

- `r(125) <= r <= 4`: `quotient_free_arb`, `research\riemann\h1319_full_xi_segment_interval_certificate.json`, PASS=`True`.
- `4 <= r <= 5`: `score_mean_value_arb`, `research\riemann\h1319_full_xi_score_mean_value_r_cover_4_5.json`, PASS=`True`.
- `5 <= r <= 10`: `score_mean_value_arb`, `research\riemann\h1319_full_xi_score_mean_value_r_cover_5_10.json`, PASS=`True`.
- `10 <= r <= 22`: `score_mean_value_arb`, `research\riemann\h1319_full_xi_score_mean_value_r10_to_22_arb.json`, PASS=`True`.
- `22 <= r <= infinity`: `exact_analytic_high_band`, `research\riemann\h1319_full_xi_high_band_r22_analytic_certificate.json`, PASS=`True`.

Each dependency is SHA-256 pinned in the JSON. The three numerical
middle pieces use Arb midpoint-plus-interval-derivative mean-value
enclosures; they do not infer interval coverage from sampled points.
The high piece is an exact analytic inequality from `r=22` onward.

## Scope

This closes the compact pointwise-curvature obligation consumed by
H920. It is not, on its own, a proof of the Riemann Hypothesis.

## Independent audit gate

Canonical audit: `research/riemann/h1319_score_mean_value_independent_audit_canonical.json`. PASS=`True`. The gate pins the engine, H1272, the 360 mean-value boxes, and the 140-dps precision-safe audit runner by SHA-256.
