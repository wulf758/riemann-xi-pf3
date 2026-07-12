# H943C U4 Five-Interval Rational Proof

Classification: `h943c_u4_five_interval_proof`

## Statement

Rational proof that

```text
U4(z) = 61 z^2/960 + z^4 exp(21z/8)/180 - (z/6) tanh(3z^3) < 1/100
```

for `0 <= z <= 1`.

## Method

Use the elementary lower bound

```text
tanh u >= u / sqrt(1 + u^2),
```

substitute `y = z^2`, and split `[0, 1]` at

```text
1/3,  1/2,  3/4,  4/5.
```

On each of the five subintervals, a quadratic rational upper bound dominates
`U4`; the exact respective upper maxima are

```text
[0, 1/3]   :  1687/233280
[1/3, 1/2] :   431/233280
[1/2, 3/4] :   -17/29952
[3/4, 4/5] :  -489/14080
[4/5, 1]   :   959/112500
```

All five are `< 1/100`. The quadratic upper bounds are decreasing on each
non-initial `y` interval.

## Notes Recorded

- The small-`z` scalar inequality needs no Arb certificate and no
  `E[X^5] <= 16/5` bound; `E[X^5] <= 4` suffices.
- These five rational cases are to be included in the formal H943-C closure
  note.

## Provenance

Materialized on 2026-07-09 from AXMEM card `h943c-u4-five-interval-proof`.
Rational maxima transcribed exactly from the card.
