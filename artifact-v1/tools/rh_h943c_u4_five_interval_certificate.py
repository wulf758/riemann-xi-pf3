#!/usr/bin/env python3
"""Exact certificate for the H943-C five-interval U4 bound.

The proof uses only rational arithmetic.  Exponential upper bounds are
certified by a Taylor polynomial plus a geometric bound for the omitted tail.
"""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path


Q = Fraction


def qstr(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def exp_upper_certificate(x: Q, bound: Q) -> dict[str, object]:
    """Prove exp(x) < bound with exact Taylor/geometric-tail arithmetic."""

    partial = Q(1)
    term = Q(1)
    for degree in range(1, 257):
        term *= x / degree
        partial += term
        first_omitted = term * x / (degree + 1)
        tail_ratio = x / (degree + 2)
        if tail_ratio >= 1:
            continue
        upper = partial + first_omitted / (1 - tail_ratio)
        if upper < bound:
            return {
                "x": qstr(x),
                "bound": qstr(bound),
                "taylor_degree": degree,
                "exact_upper": qstr(upper),
                "tail_ratio": qstr(tail_ratio),
                "check": True,
            }
    raise AssertionError(f"no Taylor certificate found for exp({x}) < {bound}")


def main() -> None:
    documented_exp_bounds = [
        (Q(7, 8), Q(657, 272), Q(5, 2)),
        (Q(21, 16), Q(123, 32), Q(4)),
        (Q(63, 32), Q(986581, 133120), Q(8)),
        (Q(21, 10), Q(32507, 3800), Q(9)),
        (Q(21, 8), Q(2278981, 163840), Q(14)),
    ]
    exp_checks = []
    for x, documented, coarse in documented_exp_bounds:
        cert = exp_upper_certificate(x, documented)
        cert["documented_lt_coarse"] = documented < coarse
        cert["coarse_bound"] = qstr(coarse)
        exp_checks.append(cert)

    # (z_left, z_right, coarse exp bound, hyperbolic coefficient lower bound,
    #  exact upper maximum recorded by the original H943-C proof).
    intervals = [
        (Q(0), Q(1, 3), Q(5, 2), Q(0), Q(1687, 233280)),
        (Q(1, 3), Q(1, 2), Q(4), Q(4, 9), Q(431, 233280)),
        (Q(1, 2), Q(3, 4), Q(8), Q(4, 13), Q(-17, 29952)),
        (Q(3, 4), Q(4, 5), Q(9), Q(3, 11), Q(-489, 14080)),
        (Q(4, 5), Q(1), Q(14), Q(5, 32), Q(959, 112500)),
    ]

    interval_checks = []
    for index, (z_left, z_right, exp_coarse, c_hyp, expected) in enumerate(
        intervals, start=1
    ):
        y_left = z_left**2
        y_right = z_right**2
        quadratic = exp_coarse / 180 - c_hyp
        derivative_left = Q(61, 960) + 2 * quadratic * y_left
        derivative_right = Q(61, 960) + 2 * quadratic * y_right

        if index == 1:
            monotonicity = derivative_left >= 0 and derivative_right >= 0
            maximizing_y = y_right
            endpoint_rule = "increasing_on_interval"
            hyperbolic_check = True  # The negative term is simply discarded.
        else:
            monotonicity = quadratic < 0 and derivative_left <= 0
            maximizing_y = y_left
            endpoint_rule = "decreasing_on_interval"
            # C <= 1/(2*sqrt(1+9*z_right^6)), checked after squaring.
            hyperbolic_check = (2 * c_hyp) ** 2 * (1 + 9 * z_right**6) <= 1

        maximum = Q(61, 960) * maximizing_y + quadratic * maximizing_y**2
        checks = {
            "hyperbolic_coefficient": hyperbolic_check,
            "quadratic_monotonicity": monotonicity,
            "maximum_matches_record": maximum == expected,
            "maximum_lt_1_over_100": maximum < Q(1, 100),
        }
        interval_checks.append(
            {
                "index": index,
                "z_interval": [qstr(z_left), qstr(z_right)],
                "y_interval": [qstr(y_left), qstr(y_right)],
                "exp_coarse_bound": qstr(exp_coarse),
                "hyperbolic_coefficient_lower": qstr(c_hyp),
                "quadratic_y2_coefficient": qstr(quadratic),
                "derivative_at_left": qstr(derivative_left),
                "derivative_at_right": qstr(derivative_right),
                "endpoint_rule": endpoint_rule,
                "maximum": qstr(maximum),
                "checks": checks,
            }
        )

    coverage = (
        intervals[0][0] == 0
        and intervals[-1][1] == 1
        and all(intervals[i][1] == intervals[i + 1][0] for i in range(4))
    )
    all_checks = (
        coverage
        and all(bool(item["check"] and item["documented_lt_coarse"]) for item in exp_checks)
        and all(all(row["checks"].values()) for row in interval_checks)
    )

    report = {
        "id": "h943c_u4_five_interval_certificate",
        "statement": (
            "U4(z)=61*z^2/960+z^4*exp(21*z/8)/180-"
            "z*tanh(3*z^3)/6 < 1/100 on 0<=z<=1"
        ),
        "analytic_reduction": {
            "tanh_bound": "tanh(u)>=u/sqrt(1+u^2) for u>=0",
            "proof": (
                "after squaring, the inequality reduces to sinh(u)^2>=u^2; "
                "this follows from sinh(u)>=u"
            ),
            "substitution": "y=z^2",
        },
        "exponential_certificates": exp_checks,
        "interval_certificates": interval_checks,
        "coverage_check": coverage,
        "all_checks_pass": all_checks,
        "scope": "exact scalar H943-C U4 certificate; parameter polynomials are separate",
    }

    output = Path(__file__).resolve().parents[1] / "research" / "riemann" / (
        "h943c_u4_five_interval_certificate.json"
    )
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    if not all_checks:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
