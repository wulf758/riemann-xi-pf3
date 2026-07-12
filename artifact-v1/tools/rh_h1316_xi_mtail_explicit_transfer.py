#!/usr/bin/env python3
"""Mechanical checks for the H1316 explicit Xi m>=2 tail constants.

The mathematical inputs M2<=26/25, M4<=7/2, M6<=21 are proved by H1313.
This script checks the exact geometric-majorant arithmetic, the elementary
threshold estimates, and the quotient-rule coefficient count used in H1316.
"""

from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path


def main() -> None:
    # At y<=1/1000 the two geometric sums in the pointwise majorant are
    # bounded by these exact rational numbers.
    y = Fraction(1, 1000)
    fourth_power_sum = Fraction(16, 1) / (1 - 16 * y)
    square_sum_with_c5 = Fraction(2, 1) / (1 - 4 * y)
    pointwise_constant = fourth_power_sum + square_sum_with_c5

    # exp(5/2)>12 from its degree-six Taylor polynomial, all exact.
    x0 = Fraction(5, 2)
    taylor6 = sum((x0**n) / math.factorial(n) for n in range(7))
    saddle_ratio_lower = 3 * 12 - Fraction(9, 4) - Fraction(2, 5)

    # M2, M4 imply E|W|^3<2; M6 is the Markov input.
    mu3_squared = Fraction(26, 25) * Fraction(7, 2)
    bad_event_coefficient = 19 * 21

    # Quotient-rule coefficient count.  Each denominator moment contributes
    # a factor at most 2; the lists are the resulting multiples of C=437.
    base = 38 + bad_event_coefficient
    first_multipliers = [1, 2]
    second_multipliers = [1, 4, 8, 2]
    third_multipliers = [1, 6, 24, 6, 48, 24, 2]
    derivative_constants = {
        "B1_tail": base * sum(first_multipliers),
        "B2_tail": base * sum(second_multipliers),
        "B3_tail": base * sum(third_multipliers),
    }

    r0_upper = 19 * math.exp(-30.0) + 399 / (30.0**3)
    scalar_maxima = {}
    for k in range(4):
        if k == 0:
            maximum = 1.0
        else:
            # max_x x^(k/2) exp(-3x/e) = (k/6)^(k/2).
            maximum = (k / 6.0) ** (k / 2.0)
        scalar_maxima[str(k)] = maximum

    root = Path(__file__).resolve().parents[1]
    dependency = root / "research" / "riemann" / "h1313_h946_global_even_moment_bounds.md"
    dependency_text = dependency.read_text(encoding="utf-8")
    dependency_checks = {
        "file_exists": dependency.exists(),
        "M2_bound_present": "E_r[W^2] <= 26/25" in dependency_text,
        "M4_bound_present": "E_r[W^4] <= 7/2" in dependency_text,
        "M6_bound_present": "E_r[W^6] <= 21" in dependency_text,
        "range_present": "every `r>=5/2`" in dependency_text,
    }

    checks = {
        "pointwise_constant_lt_19": pointwise_constant < 19,
        "taylor6_exp_5_over_2_gt_12": taylor6 > 12,
        "saddle_ratio_lower_gt_30": saddle_ratio_lower > 30,
        "mu3_lt_2": mu3_squared < 4,
        "bad_event_coefficient_is_399": bad_event_coefficient == 399,
        "base_is_437": base == 437,
        "R_tail_lt_1_over_60": r0_upper < 1 / 60,
        "scalar_maxima_le_1": all(v <= 1 for v in scalar_maxima.values()),
        "B1_is_1311": derivative_constants["B1_tail"] == 1311,
        "B2_is_6555": derivative_constants["B2_tail"] == 6555,
        "B3_is_48507": derivative_constants["B3_tail"] == 48507,
        "dependency_checks_pass": all(dependency_checks.values()),
    }

    result = {
        "id": "h1316_xi_mtail_explicit_transfer",
        "pointwise_geometric_constant": {
            "exact": str(pointwise_constant),
            "float": float(pointwise_constant),
            "target": 19,
        },
        "taylor6_exp_5_over_2": {
            "exact": str(taylor6),
            "float": float(taylor6),
        },
        "saddle_ratio_lower": {
            "exact": str(saddle_ratio_lower),
            "float": float(saddle_ratio_lower),
        },
        "R_tail_scalar_upper": r0_upper,
        "scalar_maxima": scalar_maxima,
        "derivative_constants": derivative_constants,
        "dependency_checks": dependency_checks,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    print(json.dumps(result, indent=2, sort_keys=True))

    if not result["all_checks_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
