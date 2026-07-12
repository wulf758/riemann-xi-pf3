# H807 Essential-Minor And Hard-Solver Audit

Classification: `essential_minor_theorem_not_found_hard_counterexample_not_found`

## Candidate

Positive Toeplitz sequence + PF2 + nondecreasing q_n=R_n/R_{n-1} + all translated contiguous order-3 Toeplitz minors nonnegative implies PF3.

## Source Audit

| source | finding | impact |
| --- | --- | --- |
| [Gasca and Pena / Neville-elimination criterion as summarized in Launois-Lenagan-Rigal](https://arxiv.org/abs/1207.3613) | Finite matrix total positivity can be tested by special initial minors, but the criterion asks for a broad family of initial minors, not just PF2, ratio-log-convexity, and translated contiguous D3 minors. | Does not close H807 candidate; points to missing essential minors. |
| [Katkova and Vishnyakova, A sufficient condition for a sequence to be a Pólya frequency sequence](https://arxiv.org/abs/math/0504183) | Strong adjacent dominance conditions can imply PF_m, but for m=3 the available sufficient condition is far stronger than the Xi tail data, where q_n tends to 1. | Useful benchmark, too strong for current Xi route. |
| [Pólya/de Bruijn-Newman kernel PF-order warning](https://arxiv.org/abs/2602.20313) | The continuous de Bruijn-Newman kernel is reported not to be PF5. This warns against continuous-kernel total-positivity shortcuts, but does not directly address the discrete Xi coefficient sequence. | Prevents overclaiming kernel TP; not an invalidation of JP2'. |

## Hard Search

Parameters:

```json
{
  "length": 7,
  "denominator": 512,
  "attempts": 120000,
  "seconds": 25.0,
  "batch_size": 4096,
  "seed": 807,
  "float_tolerance": 0.0
}
```

Stats:

```json
{
  "checked": 120000,
  "contiguous_d3_nonnegative_float_candidates": 87042,
  "elapsed_seconds": 3.4970982000086224,
  "best_float_margin": -2.335028279454195e-28,
  "best_float_margin_rows": [
    0,
    1,
    6
  ],
  "best_float_margin_cols": [
    6,
    7,
    8
  ],
  "best_q_ints": [
    189,
    197,
    220,
    332,
    494,
    512,
    512
  ],
  "best_q_exact_diagnosis": {
    "q_values": [
      {
        "fraction": "189/512",
        "float": 0.369140625
      },
      {
        "fraction": "197/512",
        "float": 0.384765625
      },
      {
        "fraction": "55/128",
        "float": 0.4296875
      },
      {
        "fraction": "83/128",
        "float": 0.6484375
      },
      {
        "fraction": "247/256",
        "float": 0.96484375
      },
      {
        "fraction": "1",
        "float": 1.0
      },
      {
        "fraction": "1",
        "float": 1.0
      }
    ],
    "contiguous_d3_negative_count": 0,
    "first_contiguous_d3_negative": null,
    "all_order3_negative_count": 0,
    "minimum_order3_minor": {
      "rows": [
        0,
        1,
        2
      ],
      "cols": [
        0,
        7,
        8
      ],
      "determinant": {
        "fraction": "0",
        "float": 0.0
      }
    }
  }
}
```

No certified hard counterexample was found in this bounded run.


## Decision

H807 did not find a primary theorem proving PF2+q+contiguous-D3 => PF3, and the faster bounded search did not find a certified hard counterexample. Treat the candidate as open but unsupported; next progress needs either a formal order-3 symbolic proof attempt or a stronger nonlinear solver.

## Limitations

- The source audit is a targeted primary-source scout, not a complete literature review.
- The hard search is bounded and randomized over rational grids.
- No-hit evidence is not proof of the candidate statement.
- The candidate statement is about PF3 only, not PF-infinity or RH.

## Next Target

`H808 symbolic PF3 cone or nonlinear certificate`: For order 3, symbolically reduce sparse Toeplitz minors under PF2+q+D3, or use a stronger nonlinear solver to find a certified hard counterexample.
