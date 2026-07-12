# H818 D3 Shield Audit

Classification: `h818_mixed_d3_shield_requires_global_identity`

## Question

For H815 four-term failures that are negative under monotonicity alone, are the negative grid points eliminated by one translated D3 constraint or only by a multi-D3 interaction?

## Stats

```json
{
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
}
```

## Examples

```json
{
  "failure_index": 6,
  "rows": [
    0,
    1,
    3
  ],
  "cols": [
    2,
    3,
    5
  ],
  "family": "row_gaps=(1,2), col_gaps=(1,2)",
  "polynomial": "1 - q3 - q2*q3^2*q4 + q2*q3^3*q4^2*q5",
  "negative_count": 8,
  "negative_with_all_d3_nonnegative": 0,
  "single_d3_nonnegative_survivor_counts": {
    "1": 8,
    "2": 0,
    "3": 8,
    "4": 8,
    "5": 8,
    "6": 8
  },
  "best_negative": {
    "float_value": -0.010020226240158081,
    "exact_value": {
      "fraction": "-336223/33554432",
      "float": -0.010020226240158081
    },
    "q_values": [
      {
        "name": "q2",
        "fraction": "1/2",
        "float": 0.5
      },
      {
        "name": "q3",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q4",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q5",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q6",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q7",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q8",
        "fraction": "15/16",
        "float": 0.9375
      }
    ],
    "d3_values": {
      "D3_start_1": {
        "fraction": "15/64",
        "float": 0.234375,
        "negative": false
      },
      "D3_start_2": {
        "fraction": "-193/8192",
        "float": -0.0235595703125,
        "negative": true
      },
      "D3_start_3": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_4": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_5": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_6": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      }
    },
    "violated_d3_starts": [
      2
    ]
  },
  "single_d3_shields": [
    2
  ],
  "has_single_d3_shield": true,
  "all_grid_negatives_blocked_by_all_d3": true
}
```

```json
{
  "failure_index": 7,
  "rows": [
    0,
    1,
    3
  ],
  "cols": [
    2,
    3,
    6
  ],
  "family": "row_gaps=(1,2), col_gaps=(1,3)",
  "polynomial": "1 - q3 - q2*q3^2*q4^2*q5 + q2*q3^3*q4^3*q5^2*q6",
  "negative_count": 7,
  "negative_with_all_d3_nonnegative": 0,
  "single_d3_nonnegative_survivor_counts": {
    "1": 7,
    "2": 0,
    "3": 7,
    "4": 7,
    "5": 7,
    "6": 7
  },
  "best_negative": {
    "float_value": -0.019885963651177008,
    "exact_value": {
      "fraction": "-2733106033/137438953472",
      "float": -0.019885963651177008
    },
    "q_values": [
      {
        "name": "q2",
        "fraction": "1/2",
        "float": 0.5
      },
      {
        "name": "q3",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q4",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q5",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q6",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q7",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q8",
        "fraction": "15/16",
        "float": 0.9375
      }
    ],
    "d3_values": {
      "D3_start_1": {
        "fraction": "15/64",
        "float": 0.234375,
        "negative": false
      },
      "D3_start_2": {
        "fraction": "-193/8192",
        "float": -0.0235595703125,
        "negative": true
      },
      "D3_start_3": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_4": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_5": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_6": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      }
    },
    "violated_d3_starts": [
      2
    ]
  },
  "single_d3_shields": [
    2
  ],
  "has_single_d3_shield": true,
  "all_grid_negatives_blocked_by_all_d3": true
}
```

```json
{
  "failure_index": 8,
  "rows": [
    0,
    1,
    3
  ],
  "cols": [
    2,
    3,
    7
  ],
  "family": "row_gaps=(1,2), col_gaps=(1,4)",
  "polynomial": "1 - q3 - q2*q3^2*q4^2*q5^2*q6 + q2*q3^3*q4^3*q5^3*q6^2*q7",
  "negative_count": 8,
  "negative_with_all_d3_nonnegative": 0,
  "single_d3_nonnegative_survivor_counts": {
    "1": 8,
    "2": 0,
    "3": 8,
    "4": 8,
    "5": 8,
    "6": 8
  },
  "best_negative": {
    "float_value": -0.025274591345235464,
    "exact_value": {
      "fraction": "-14228330020543/562949953421312",
      "float": -0.025274591345235464
    },
    "q_values": [
      {
        "name": "q2",
        "fraction": "1/2",
        "float": 0.5
      },
      {
        "name": "q3",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q4",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q5",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q6",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q7",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q8",
        "fraction": "15/16",
        "float": 0.9375
      }
    ],
    "d3_values": {
      "D3_start_1": {
        "fraction": "15/64",
        "float": 0.234375,
        "negative": false
      },
      "D3_start_2": {
        "fraction": "-193/8192",
        "float": -0.0235595703125,
        "negative": true
      },
      "D3_start_3": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_4": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_5": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_6": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      }
    },
    "violated_d3_starts": [
      2
    ]
  },
  "single_d3_shields": [
    2
  ],
  "has_single_d3_shield": true,
  "all_grid_negatives_blocked_by_all_d3": true
}
```

```json
{
  "failure_index": 9,
  "rows": [
    0,
    1,
    3
  ],
  "cols": [
    2,
    3,
    8
  ],
  "family": "row_gaps=(1,2), col_gaps=(1,5)",
  "polynomial": "1 - q3 - q2*q3^2*q4^2*q5^2*q6^2*q7 + q2*q3^3*q4^3*q5^3*q6^3*q7^2*q8",
  "negative_count": 6,
  "negative_with_all_d3_nonnegative": 0,
  "single_d3_nonnegative_survivor_counts": {
    "1": 6,
    "2": 0,
    "3": 6,
    "4": 6,
    "5": 6,
    "6": 6
  },
  "best_negative": {
    "float_value": -0.027306050451698205,
    "exact_value": {
      "fraction": "-62963465543284753/2305843009213693952",
      "float": -0.027306050451698212
    },
    "q_values": [
      {
        "name": "q2",
        "fraction": "1/2",
        "float": 0.5
      },
      {
        "name": "q3",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q4",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q5",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q6",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q7",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q8",
        "fraction": "15/16",
        "float": 0.9375
      }
    ],
    "d3_values": {
      "D3_start_1": {
        "fraction": "15/64",
        "float": 0.234375,
        "negative": false
      },
      "D3_start_2": {
        "fraction": "-193/8192",
        "float": -0.0235595703125,
        "negative": true
      },
      "D3_start_3": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_4": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_5": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_6": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      }
    },
    "violated_d3_starts": [
      2
    ]
  },
  "single_d3_shields": [
    2
  ],
  "has_single_d3_shield": true,
  "all_grid_negatives_blocked_by_all_d3": true
}
```

```json
{
  "failure_index": 10,
  "rows": [
    0,
    1,
    3
  ],
  "cols": [
    2,
    4,
    6
  ],
  "family": "row_gaps=(1,2), col_gaps=(2,2)",
  "polynomial": "1 - q3*q4 - q3*q4^2*q5 + q3^2*q4^3*q5^2*q6",
  "negative_count": 104,
  "negative_with_all_d3_nonnegative": 0,
  "single_d3_nonnegative_survivor_counts": {
    "1": 104,
    "2": 0,
    "3": 80,
    "4": 104,
    "5": 104,
    "6": 96
  },
  "best_negative": {
    "float_value": -0.05466297245584428,
    "exact_value": {
      "fraction": "-234775679/4294967296",
      "float": -0.05466297245584428
    },
    "q_values": [
      {
        "name": "q2",
        "fraction": "1/16",
        "float": 0.0625
      },
      {
        "name": "q3",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q4",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q5",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q6",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q7",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q8",
        "fraction": "15/16",
        "float": 0.9375
      }
    ],
    "d3_values": {
      "D3_start_1": {
        "fraction": "3599/4096",
        "float": 0.878662109375,
        "negative": false
      },
      "D3_start_2": {
        "fraction": "-3119/65536",
        "float": -0.0475921630859375,
        "negative": true
      },
      "D3_start_3": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_4": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_5": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_6": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      }
    },
    "violated_d3_starts": [
      2
    ]
  },
  "single_d3_shields": [
    2
  ],
  "has_single_d3_shield": true,
  "all_grid_negatives_blocked_by_all_d3": true
}
```

```json
{
  "failure_index": 11,
  "rows": [
    0,
    1,
    3
  ],
  "cols": [
    2,
    4,
    7
  ],
  "family": "row_gaps=(1,2), col_gaps=(2,3)",
  "polynomial": "1 - q3*q4 - q3*q4^2*q5^2*q6 + q3^2*q4^3*q5^3*q6^2*q7",
  "negative_count": 48,
  "negative_with_all_d3_nonnegative": 0,
  "single_d3_nonnegative_survivor_counts": {
    "1": 48,
    "2": 0,
    "3": 32,
    "4": 48,
    "5": 48,
    "6": 48
  },
  "best_negative": {
    "float_value": -0.06615871153650232,
    "exact_value": {
      "fraction": "-1163876361809/17592186044416",
      "float": -0.06615871153650232
    },
    "q_values": [
      {
        "name": "q2",
        "fraction": "1/16",
        "float": 0.0625
      },
      {
        "name": "q3",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q4",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q5",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q6",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q7",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q8",
        "fraction": "15/16",
        "float": 0.9375
      }
    ],
    "d3_values": {
      "D3_start_1": {
        "fraction": "3599/4096",
        "float": 0.878662109375,
        "negative": false
      },
      "D3_start_2": {
        "fraction": "-3119/65536",
        "float": -0.0475921630859375,
        "negative": true
      },
      "D3_start_3": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_4": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_5": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_6": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      }
    },
    "violated_d3_starts": [
      2
    ]
  },
  "single_d3_shields": [
    2
  ],
  "has_single_d3_shield": true,
  "all_grid_negatives_blocked_by_all_d3": true
}
```

```json
{
  "failure_index": 12,
  "rows": [
    0,
    1,
    3
  ],
  "cols": [
    2,
    4,
    8
  ],
  "family": "row_gaps=(1,2), col_gaps=(2,4)",
  "polynomial": "1 - q3*q4 - q3*q4^2*q5^2*q6^2*q7 + q3^2*q4^3*q5^3*q6^3*q7^2*q8",
  "negative_count": 48,
  "negative_with_all_d3_nonnegative": 0,
  "single_d3_nonnegative_survivor_counts": {
    "1": 48,
    "2": 0,
    "3": 32,
    "4": 48,
    "5": 48,
    "6": 48
  },
  "best_negative": {
    "float_value": -0.07049249096362287,
    "exact_value": {
      "fraction": "-5079519296579039/72057594037927936",
      "float": -0.07049249096362285
    },
    "q_values": [
      {
        "name": "q2",
        "fraction": "1/16",
        "float": 0.0625
      },
      {
        "name": "q3",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q4",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q5",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q6",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q7",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q8",
        "fraction": "15/16",
        "float": 0.9375
      }
    ],
    "d3_values": {
      "D3_start_1": {
        "fraction": "3599/4096",
        "float": 0.878662109375,
        "negative": false
      },
      "D3_start_2": {
        "fraction": "-3119/65536",
        "float": -0.0475921630859375,
        "negative": true
      },
      "D3_start_3": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_4": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_5": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_6": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      }
    },
    "violated_d3_starts": [
      2
    ]
  },
  "single_d3_shields": [
    2
  ],
  "has_single_d3_shield": true,
  "all_grid_negatives_blocked_by_all_d3": true
}
```

```json
{
  "failure_index": 13,
  "rows": [
    0,
    1,
    3
  ],
  "cols": [
    2,
    5,
    6
  ],
  "family": "row_gaps=(1,2), col_gaps=(3,1)",
  "polynomial": "1 - q4*q5 - q3*q4*q5 + q3*q4^2*q5^2*q6",
  "negative_count": 72,
  "negative_with_all_d3_nonnegative": 0,
  "single_d3_nonnegative_survivor_counts": {
    "1": 72,
    "2": 18,
    "3": 24,
    "4": 72,
    "5": 72,
    "6": 72
  },
  "best_negative": {
    "float_value": -0.023946702480316162,
    "exact_value": {
      "fraction": "-401759/16777216",
      "float": -0.023946702480316162
    },
    "q_values": [
      {
        "name": "q2",
        "fraction": "1/16",
        "float": 0.0625
      },
      {
        "name": "q3",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q4",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q5",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q6",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q7",
        "fraction": "15/16",
        "float": 0.9375
      },
      {
        "name": "q8",
        "fraction": "15/16",
        "float": 0.9375
      }
    ],
    "d3_values": {
      "D3_start_1": {
        "fraction": "3599/4096",
        "float": 0.878662109375,
        "negative": false
      },
      "D3_start_2": {
        "fraction": "-3119/65536",
        "float": -0.0475921630859375,
        "negative": true
      },
      "D3_start_3": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_4": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_5": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      },
      "D3_start_6": {
        "fraction": "31/65536",
        "float": 0.0004730224609375,
        "negative": false
      }
    },
    "violated_d3_starts": [
      2
    ]
  },
  "single_d3_shields": [],
  "has_single_d3_shield": false,
  "all_grid_negatives_blocked_by_all_d3": true
}
```

## Decision

All monotone negatives are blocked by D3 collectively, but not always by a single D3. This points to a Plucker-Dodgson or multi-D3 determinant identity.

## Next Target

`H819 D3 jump-barrier lemma or Plucker-Dodgson`: If single-D3 shields dominate, derive a local jump-barrier inequality from D3_k>=0. Otherwise search determinant identities combining several D3 constraints.
