# H819 Adjacent D3 Barrier

Classification: `h819_adjacent_two_d3_barrier_grid_counterexample`

## Question

Do the 20 H818 no-single-shield four-term failures reduce to an adjacent two-D3 barrier using only D3_start_2 and D3_start_3?

## H818 Leftovers

```json
{
  "h818_stats": {
    "grid_total": 170544,
    "checked": 170544,
    "q2_le_half_candidates": 167112,
    "h815_four_term_failures": 200,
    "failures_with_monotone_negative": 102,
    "failures_with_negative_surviving_all_d3": 0,
    "failures_with_single_d3_shield": 82,
    "single_d3_shield_counts": {
      "1": 0,
      "2": 82,
      "3": 0,
      "4": 0,
      "5": 0,
      "6": 0
    }
  },
  "leftover_count": 20,
  "unique_count": 6,
  "leftover_family_counts": [
    [
      "row_gaps=(1,2), col_gaps=(3,1)",
      3
    ],
    [
      "row_gaps=(1,3), col_gaps=(2,1)",
      3
    ],
    [
      "row_gaps=(1,2), col_gaps=(3,2)",
      2
    ],
    [
      "row_gaps=(1,3), col_gaps=(2,2)",
      2
    ],
    [
      "row_gaps=(2,2), col_gaps=(3,1)",
      2
    ],
    [
      "row_gaps=(2,3), col_gaps=(2,1)",
      2
    ],
    [
      "row_gaps=(1,3), col_gaps=(2,3)",
      1
    ],
    [
      "row_gaps=(1,3), col_gaps=(3,2)",
      1
    ],
    [
      "row_gaps=(2,2), col_gaps=(3,2)",
      1
    ],
    [
      "row_gaps=(2,3), col_gaps=(2,2)",
      1
    ],
    [
      "row_gaps=(2,3), col_gaps=(3,1)",
      1
    ],
    [
      "row_gaps=(3,2), col_gaps=(3,1)",
      1
    ]
  ]
}
```

## Unique Polynomials

- P0: `1 - q4*q5 - q3*q4*q5 + q3*q4^2*q5^2*q6`
  - multiplicity: `6`
  - active D3 starts among negatives: `[2, 3]`
- P1: `1 - q4*q5^2*q6 - q3*q4*q5 + q3*q4^2*q5^3*q6^2*q7`
  - multiplicity: `4`
  - active D3 starts among negatives: `[2, 3]`
- P2: `1 - q4*q5 - q3*q4^2*q5^2*q6 + q3*q4^3*q5^3*q6^2*q7`
  - multiplicity: `4`
  - active D3 starts among negatives: `[2, 3]`
- P3: `1 - q4*q5 - q3*q4^2*q5^3*q6^2*q7 + q3*q4^3*q5^4*q6^3*q7^2*q8`
  - multiplicity: `2`
  - active D3 starts among negatives: `[2, 3]`
- P4: `1 - q4*q5*q6 - q4*q5^2*q6^2*q7 + q4^2*q5^3*q6^3*q7^2*q8`
  - multiplicity: `2`
  - active D3 starts among negatives: `[2, 3]`
- P5: `1 - q4*q5^2*q6 - q3*q4^2*q5^2*q6 + q3*q4^3*q5^4*q6^3*q7^2*q8`
  - multiplicity: `2`
  - active D3 starts among negatives: `[2, 3]`

## Scan Summary

```json
{
  "classification": "h819_adjacent_two_d3_barrier_grid_counterexample",
  "decision": "The adjacent pair D3_2,D3_3 is not enough on the tested grids. Pivot to a wider Plucker-Dodgson identity or Xi-specific input.",
  "pair_has_negative": true,
  "d3_2_only_has_negative": true,
  "d3_3_only_has_negative": true,
  "all_d3_has_negative": false,
  "worst_pair_case": {
    "denominator": 32,
    "polynomial_key": "4",
    "polynomial": "1 - q4*q5*q6 - q4*q5^2*q6^2*q7 + q4^2*q5^3*q6^3*q7^2*q8",
    "best_float_value": -0.0065553787778835915,
    "best_exact": {
      "fraction": "-236182411370793/36028797018963968",
      "float": -0.006555378777883619
    },
    "best_q_values": [
      {
        "name": "q2",
        "fraction": "1/32",
        "float": 0.03125
      },
      {
        "name": "q3",
        "fraction": "11/16",
        "float": 0.6875
      },
      {
        "name": "q4",
        "fraction": "29/32",
        "float": 0.90625
      },
      {
        "name": "q5",
        "fraction": "31/32",
        "float": 0.96875
      },
      {
        "name": "q6",
        "fraction": "31/32",
        "float": 0.96875
      },
      {
        "name": "q7",
        "fraction": "31/32",
        "float": 0.96875
      },
      {
        "name": "q8",
        "fraction": "31/32",
        "float": 0.96875
      }
    ]
  }
}
```

## Compact Scan Table

```json
[
  {
    "denominator": 16,
    "conditions": {
      "monotone_q2_half": 167112,
      "d3_2_only": 164947,
      "d3_3_only": 161521,
      "d3_2_and_d3_3": 159486,
      "all_d3": 113542
    },
    "polynomials": [
      {
        "key": 0,
        "polynomial": "1 - q4*q5 - q3*q4*q5 + q3*q4^2*q5^2*q6",
        "negative_counts": {
          "monotone_q2_half": 72,
          "d3_2_only": 18,
          "d3_3_only": 24,
          "d3_2_and_d3_3": 0,
          "all_d3": 0
        },
        "best_d3_2_and_d3_3": {
          "fraction": "0",
          "float": 0.0
        }
      },
      {
        "key": 1,
        "polynomial": "1 - q4*q5^2*q6 - q3*q4*q5 + q3*q4^2*q5^3*q6^2*q7",
        "negative_counts": {
          "monotone_q2_half": 48,
          "d3_2_only": 12,
          "d3_3_only": 16,
          "d3_2_and_d3_3": 0,
          "all_d3": 0
        },
        "best_d3_2_and_d3_3": {
          "fraction": "0",
          "float": 0.0
        }
      },
      {
        "key": 2,
        "polynomial": "1 - q4*q5 - q3*q4^2*q5^2*q6 + q3*q4^3*q5^3*q6^2*q7",
        "negative_counts": {
          "monotone_q2_half": 64,
          "d3_2_only": 28,
          "d3_3_only": 16,
          "d3_2_and_d3_3": 0,
          "all_d3": 0
        },
        "best_d3_2_and_d3_3": {
          "fraction": "0",
          "float": 0.0
        }
      },
      {
        "key": 3,
        "polynomial": "1 - q4*q5 - q3*q4^2*q5^3*q6^2*q7 + q3*q4^3*q5^4*q6^3*q7^2*q8",
        "negative_counts": {
          "monotone_q2_half": 40,
          "d3_2_only": 14,
          "d3_3_only": 16,
          "d3_2_and_d3_3": 0,
          "all_d3": 0
        },
        "best_d3_2_and_d3_3": {
          "fraction": "0",
          "float": 0.0
        }
      },
      {
        "key": 4,
        "polynomial": "1 - q4*q5*q6 - q4*q5^2*q6^2*q7 + q4^2*q5^3*q6^3*q7^2*q8",
        "negative_counts": {
          "monotone_q2_half": 92,
          "d3_2_only": 74,
          "d3_3_only": 8,
          "d3_2_and_d3_3": 0,
          "all_d3": 0
        },
        "best_d3_2_and_d3_3": {
          "fraction": "0",
          "float": 0.0
        }
      },
      {
        "key": 5,
        "polynomial": "1 - q4*q5^2*q6 - q3*q4^2*q5^2*q6 + q3*q4^3*q5^4*q6^3*q7^2*q8",
        "negative_counts": {
          "monotone_q2_half": 40,
          "d3_2_only": 6,
          "d3_3_only": 16,
          "d3_2_and_d3_3": 0,
          "all_d3": 0
        },
        "best_d3_2_and_d3_3": {
          "fraction": "0",
          "float": 0.0
        }
      }
    ]
  },
  {
    "denominator": 24,
    "conditions": {
      "monotone_q2_half": 2003976,
      "d3_2_only": 1985609,
      "d3_3_only": 1950005,
      "d3_2_and_d3_3": 1932775,
      "all_d3": 1350040
    },
    "polynomials": [
      {
        "key": 0,
        "polynomial": "1 - q4*q5 - q3*q4*q5 + q3*q4^2*q5^2*q6",
        "negative_counts": {
          "monotone_q2_half": 504,
          "d3_2_only": 108,
          "d3_3_only": 180,
          "d3_2_and_d3_3": 0,
          "all_d3": 0
        },
        "best_d3_2_and_d3_3": {
          "fraction": "0",
          "float": 0.0
        }
      },
      {
        "key": 1,
        "polynomial": "1 - q4*q5^2*q6 - q3*q4*q5 + q3*q4^2*q5^3*q6^2*q7",
        "negative_counts": {
          "monotone_q2_half": 264,
          "d3_2_only": 48,
          "d3_3_only": 120,
          "d3_2_and_d3_3": 0,
          "all_d3": 0
        },
        "best_d3_2_and_d3_3": {
          "fraction": "0",
          "float": 0.0
        }
      },
      {
        "key": 2,
        "polynomial": "1 - q4*q5 - q3*q4^2*q5^2*q6 + q3*q4^3*q5^3*q6^2*q7",
        "negative_counts": {
          "monotone_q2_half": 348,
          "d3_2_only": 120,
          "d3_3_only": 120,
          "d3_2_and_d3_3": 0,
          "all_d3": 0
        },
        "best_d3_2_and_d3_3": {
          "fraction": "0",
          "float": 0.0
        }
      },
      {
        "key": 3,
        "polynomial": "1 - q4*q5 - q3*q4^2*q5^3*q6^2*q7 + q3*q4^3*q5^4*q6^3*q7^2*q8",
        "negative_counts": {
          "monotone_q2_half": 240,
          "d3_2_only": 108,
          "d3_3_only": 84,
          "d3_2_and_d3_3": 0,
          "all_d3": 0
        },
        "best_d3_2_and_d3_3": {
          "fraction": "0",
          "float": 0.0
        }
      },
      {
        "key": 4,
        "polynomial": "1 - q4*q5*q6 - q4*q5^2*q6^2*q7 + q4^2*q5^3*q6^3*q7^2*q8",
        "negative_counts": {
          "monotone_q2_half": 618,
          "d3_2_only": 506,
          "d3_3_only": 60,
          "d3_2_and_d3_3": 0,
          "all_d3": 0
        },
        "best_d3_2_and_d3_3": {
          "fraction": "0",
          "float": 0.0
        }
      },
      {
        "key": 5,
        "polynomial": "1 - q4*q5^2*q6 - q3*q4^2*q5^2*q6 + q3*q4^3*q5^4*q6^3*q7^2*q8",
        "negative_counts": {
          "monotone_q2_half": 204,
          "d3_2_only": 72,
          "d3_3_only": 84,
          "d3_2_and_d3_3": 0,
          "all_d3": 0
        },
        "best_d3_2_and_d3_3": {
          "fraction": "0",
          "float": 0.0
        }
      }
    ]
  },
  {
    "denominator": 32,
    "conditions": {
      "monotone_q2_half": 12449712,
      "d3_2_only": 12355603,
      "d3_3_only": 12170932,
      "d3_2_and_d3_3": 12080581,
      "all_d3": 8606983
    },
    "polynomials": [
      {
        "key": 0,
        "polynomial": "1 - q4*q5 - q3*q4*q5 + q3*q4^2*q5^2*q6",
        "negative_counts": {
          "monotone_q2_half": 1616,
          "d3_2_only": 297,
          "d3_3_only": 688,
          "d3_2_and_d3_3": 0,
          "all_d3": 0
        },
        "best_d3_2_and_d3_3": {
          "fraction": "0",
          "float": 0.0
        }
      },
      {
        "key": 1,
        "polynomial": "1 - q4*q5^2*q6 - q3*q4*q5 + q3*q4^2*q5^3*q6^2*q7",
        "negative_counts": {
          "monotone_q2_half": 1024,
          "d3_2_only": 189,
          "d3_3_only": 464,
          "d3_2_and_d3_3": 0,
          "all_d3": 0
        },
        "best_d3_2_and_d3_3": {
          "fraction": "0",
          "float": 0.0
        }
      },
      {
        "key": 2,
        "polynomial": "1 - q4*q5 - q3*q4^2*q5^2*q6 + q3*q4^3*q5^3*q6^2*q7",
        "negative_counts": {
          "monotone_q2_half": 1280,
          "d3_2_only": 505,
          "d3_3_only": 368,
          "d3_2_and_d3_3": 0,
          "all_d3": 0
        },
        "best_d3_2_and_d3_3": {
          "fraction": "0",
          "float": 0.0
        }
      },
      {
        "key": 3,
        "polynomial": "1 - q4*q5 - q3*q4^2*q5^3*q6^2*q7 + q3*q4^3*q5^4*q6^3*q7^2*q8",
        "negative_counts": {
          "monotone_q2_half": 928,
          "d3_2_only": 377,
          "d3_3_only": 256,
          "d3_2_and_d3_3": 0,
          "all_d3": 0
        },
        "best_d3_2_and_d3_3": {
          "fraction": "0",
          "float": 0.0
        }
      },
      {
        "key": 4,
        "polynomial": "1 - q4*q5*q6 - q4*q5^2*q6^2*q7 + q4^2*q5^3*q6^3*q7^2*q8",
        "negative_counts": {
          "monotone_q2_half": 2176,
          "d3_2_only": 1794,
          "d3_3_only": 240,
          "d3_2_and_d3_3": 62,
          "all_d3": 0
        },
        "best_d3_2_and_d3_3": {
          "fraction": "-236182411370793/36028797018963968",
          "float": -0.006555378777883619
        }
      },
      {
        "key": 5,
        "polynomial": "1 - q4*q5^2*q6 - q3*q4^2*q5^2*q6 + q3*q4^3*q5^4*q6^3*q7^2*q8",
        "negative_counts": {
          "monotone_q2_half": 816,
          "d3_2_only": 298,
          "d3_3_only": 272,
          "d3_2_and_d3_3": 0,
          "all_d3": 0
        },
        "best_d3_2_and_d3_3": {
          "fraction": "0",
          "float": 0.0
        }
      }
    ]
  }
]
```

## Underexploited Next Ideas

- `Adjacent two-D3 jump-barrier lemma` (primary): H818 leftovers only need starts 2 and 3 on the grid; H816 did not test this scalar two-barrier form.
- `Dodgson condensation identity for sparse rows` (secondary): If the scalar barrier proof stalls, express sparse minors through adjacent contiguous minors and PF2 factors.
- `Xi-specific unit-deficit input` (fallback): H817 showed monotone q alone is false; an analytic input from Xi coefficients may avoid proving a general q-cone theorem.
- `Cross-route Nyman-Beurling endpoint bridge` (separate branch): Older H520/H662 routes still contain formal candidates, especially endpoint Abel cancellation, not used by the PF3 corridor.

## Decision

The adjacent pair D3_2,D3_3 is not enough on the tested grids. Pivot to a wider Plucker-Dodgson identity or Xi-specific input.

## Next Target

`H820 formal adjacent D3 lemma`: Try to prove the six unique normalized four-term forms from monotonicity, q2<=1/2, D3_2>=0, and D3_3>=0, or produce a continuous counterexample.
