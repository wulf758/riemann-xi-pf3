#!/usr/bin/env python3
"""Quadratic-in-r Arb Taylor pilot for one H1319 full-Xi box.

This prototype repairs the main dependency loss of a direct r-ball.  The
central m=1 moments are expanded to order two at the rational midpoint in r;
their third Taylor coefficients are bounded on the whole r-box.  The W
midpoint remainder, H1316 arithmetic tail, and H1271 outer tail are then
added as uniform interval errors.

The result is rigorous for the configured single box, but this file is not an
adaptive cover or a whole-range certificate.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

from flint import arb, arb_series, ctx

import rh_h1319_full_xi_adaptive_rbox_cover as base


DEFAULT_OUT = Path(
    "research/riemann/h1319_full_xi_r_taylor_box_pilot.json"
)


def stable_expm1_general(value: arb | arb_series) -> arb | arb_series:
    if not isinstance(value, arb_series):
        return value.expm1()
    raw = value.exp() - 1
    coefficients = [value[0].expm1()]
    coefficients.extend(raw[index] for index in range(1, ctx.cap))
    return arb_series(coefficients)


def stable_rem2_general(value: arb | arb_series) -> arb | arb_series:
    if not isinstance(value, arb_series):
        return base.stable_rem2_arb(value)
    raw = value.exp() - 1 - value
    coefficients = [base.stable_rem2_arb(value[0])]
    coefficients.extend(raw[index] for index in range(1, ctx.cap))
    return arb_series(coefficients)


class TaylorFullXiBoxModel(base.FullXiBoxModel):
    """H1319 m=1 model whose stable Psi supports r-series to degree three."""

    def psi(self, w: arb | arb_series) -> arb | arb_series:
        h = w / self.sqrt_A
        delta = self.r * stable_expm1_general(h)
        return (
            self.pi * self.exp_r * stable_rem2_general(delta)
            + (self.pi * self.r * self.exp_r - self.alpha * self.r)
            * stable_rem2_general(h)
        )


def series_variable(constant: arb) -> arb_series:
    return arb_series([constant, arb(1), arb(0), arb(0)])


def moment_taylor_enclosures(
    r_lo: Fraction,
    r_hi: Fraction,
    width: Fraction,
    subdivisions: int,
    center: Fraction,
) -> tuple[list[arb], dict[str, list[arb] | arb]]:
    if not r_lo < r_hi:
        raise ValueError("r_lo must be strictly less than r_hi")
    midpoint_r = (r_lo + r_hi) / 2
    radius_r = (r_hi - r_lo) / 2
    x = base.ball(-radius_r, radius_r)

    model_center = TaylorFullXiBoxModel(series_variable(base.qarb(midpoint_r)))
    model_range = TaylorFullXiBoxModel(series_variable(base.ball(r_lo, r_hi)))
    model_box = TaylorFullXiBoxModel(base.ball(r_lo, r_hi))

    step = width / subdivisions
    step_a = base.qarb(step)
    w_remainder_scale = step_a**3 / 12
    center_a = base.qarb(center)
    coefficients = [[arb(0) for _ in range(3)] for _ in range(4)]
    third_coefficient_abs = [arb(0) for _ in range(4)]
    w_midpoint_error = [arb(0) for _ in range(4)]

    for index in range(subdivisions):
        w_lo = index * step
        w_hi = w_lo + step
        w_midpoint = (w_lo + w_hi) / 2

        at_center = base.paired_main_values(
            model_center, base.qarb(w_midpoint), center_a
        )
        on_r_box = base.paired_main_values(
            model_range, base.qarb(w_midpoint), center_a
        )
        w_series = series_variable(base.ball(w_lo, w_hi))
        on_w_box = base.paired_main_values(model_box, w_series, center_a)

        for order in range(4):
            for degree in range(3):
                coefficients[order][degree] += step_a * at_center[order][degree]
            third_coefficient_abs[order] += step_a * base.abs_upper(
                on_r_box[order][3]
            )
            w_midpoint_error[order] += w_remainder_scale * base.abs_upper(
                on_w_box[order][2]
            )

    radius_r_a = base.qarb(radius_r)
    r_taylor_error = [
        bound * radius_r_a**3 for bound in third_coefficient_abs
    ]
    arithmetic_tail, h_bar = base.central_arithmetic_tail_bounds(
        model_box, width, center
    )
    outer_tail = base.outer_translated_tail_bounds(center)

    moments: list[arb] = []
    for order in range(4):
        polynomial = (
            coefficients[order][0]
            + coefficients[order][1] * x
            + coefficients[order][2] * x * x
        )
        main_error = w_midpoint_error[order] + r_taylor_error[order]
        polynomial += base.symmetric(main_error)
        tail_error = arithmetic_tail[order] + outer_tail[order]
        if order in (0, 2):
            polynomial += base.zero_to(tail_error)
        else:
            polynomial += base.symmetric(tail_error)
        moments.append(polynomial)

    diagnostics: dict[str, list[arb] | arb] = {
        "coefficient_0": [row[0] for row in coefficients],
        "coefficient_1": [row[1] for row in coefficients],
        "coefficient_2": [row[2] for row in coefficients],
        "third_coefficient_abs": third_coefficient_abs,
        "r_taylor_error": r_taylor_error,
        "w_midpoint_error": w_midpoint_error,
        "arithmetic_tail": arithmetic_tail,
        "outer_tail": outer_tail,
        "h1316_hbar": h_bar,
        "t_box": model_box.t,
        "A_box": model_box.A,
    }
    return moments, diagnostics


def payload_list(values: list[arb], digits: int) -> list[dict[str, str]]:
    return [base.interval_payload(value, digits) for value in values]


def certificate(
    r_lo: Fraction,
    r_hi: Fraction,
    width: Fraction,
    subdivisions: int,
    center: Fraction,
    beta: Fraction,
    digits: int,
) -> dict[str, Any]:
    dependencies = base.validate_dependencies()
    moments, diagnostics = moment_taylor_enclosures(
        r_lo, r_hi, width, subdivisions, center
    )
    k0, k1, k2, k3 = moments
    cumulant_numerator = k3 * k0 * k0 - 3 * k1 * k2 * k0 + 2 * k1**3
    t_box = diagnostics["t_box"]
    a_box = diagnostics["A_box"]
    assert isinstance(t_box, arb) and isinstance(a_box, arb)
    denominator = a_box * a_box.sqrt() * k0**3
    scaled = 8 * t_box**2 * cumulant_numerator / denominator
    margin = scaled + base.qarb(beta)
    checks = {
        "dependencies_pass": dependencies["all_pass"],
        "mass_positive": bool(k0.lower() > 0),
        "denominator_positive": bool(denominator.lower() > 0),
        "target_margin_positive": bool(margin.lower() > 0),
    }
    return {
        "id": "h1319_full_xi_r_taylor_box_pilot",
        "classification": "h1319_single_r_box_taylor_certificate_not_cover",
        "config": {
            "r_lo": str(r_lo),
            "r_hi": str(r_hi),
            "r_width": str(r_hi - r_lo),
            "W": str(width),
            "subdivisions": subdivisions,
            "center": str(center),
            "beta": str(beta),
            "series_degree_in_r": 2,
            "r_remainder_derivative_degree": 3,
            "finite_differences_used": False,
        },
        "moments": payload_list(moments, digits),
        "cumulant_numerator": base.interval_payload(cumulant_numerator, digits),
        "scaled_t2_kappa3": base.interval_payload(scaled, digits),
        "margin": base.interval_payload(margin, digits),
        "errors": {
            "r_taylor": payload_list(
                diagnostics["r_taylor_error"], digits  # type: ignore[arg-type]
            ),
            "W_midpoint": payload_list(
                diagnostics["w_midpoint_error"], digits  # type: ignore[arg-type]
            ),
            "central_arithmetic_tail": payload_list(
                diagnostics["arithmetic_tail"], digits  # type: ignore[arg-type]
            ),
            "outer_tail": payload_list(
                diagnostics["outer_tail"], digits  # type: ignore[arg-type]
            ),
        },
        "h1316_hbar": base.interval_payload(
            diagnostics["h1316_hbar"], digits  # type: ignore[arg-type]
        ),
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "scope": "one exact rational r-box only; adaptive assembly remains open",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r-lo", default="10")
    parser.add_argument("--r-hi", default="101/10")
    parser.add_argument("--W", default="24")
    parser.add_argument("--subdivisions", type=int, default=300)
    parser.add_argument("--center", default="0")
    parser.add_argument("--beta", default="3/4")
    parser.add_argument("--dps", type=int, default=90)
    parser.add_argument("--digits", type=int, default=24)
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()

    ctx.dps = args.dps
    ctx.cap = 4
    report = certificate(
        Fraction(args.r_lo),
        Fraction(args.r_hi),
        Fraction(args.W),
        args.subdivisions,
        Fraction(args.center),
        Fraction(args.beta),
        args.digits,
    )
    output = Path(args.out)
    if not args.no_write:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    print(
        json.dumps(
            {
                "all_checks_pass": report["all_checks_pass"],
                "config": report["config"],
                "scaled_t2_kappa3": report["scaled_t2_kappa3"],
                "margin": report["margin"],
                "errors": report["errors"],
                "out": None if args.no_write else str(output),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0 if report["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
