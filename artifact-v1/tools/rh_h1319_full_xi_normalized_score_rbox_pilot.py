#!/usr/bin/env python3
"""Single r-box H1319 pilot with the full-Xi score divided by ``a``.

Factoring ``a=B/A^(3/2)`` before integration avoids multiplying tiny score
moments by an exponentially large scale after interval dependency has already
been lost.  This is a rigorous certificate for one configured rational box,
not an adaptive cover.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import flint
from flint import arb

import rh_h1319_full_xi_fixed_r_arb_pilot as raw_base
import rh_h1319_full_xi_fixed_r_score_arb as score
from rh_h930_fixed_r_arb_tail_certificate import (
    arb_box,
    arb_from_fraction,
    interval_payload,
)


DEFAULT_OUT = Path(
    "research/riemann/h1319_full_xi_normalized_score_rbox_pilot.json"
)


def central_normalized_score_sums(
    model: score.ScoreFullXiModel,
    width: Fraction,
    subdivisions: int,
) -> tuple[list[arb], dict[str, arb]]:
    step = (2 * width) / subdivisions
    step_a = arb_from_fraction(step)
    sums = [arb(0) for _ in range(4)]
    min_kernel: arb | None = None
    max_kernel = arb(0)
    max_derivative_tail = arb(0)

    for index in range(subdivisions):
        left = -width + index * step
        w = arb_box(left, left + step)
        density = model.density(w)
        kernel, kernel_prime, derivative_tail = model.kernel_ratio_and_derivative(w)
        n_dom = model.dominant_score_excess(w)
        normalized_score_numerator = (
            density * (kernel * n_dom - kernel_prime) / model.a
        )
        w2 = w * w
        sums[0] += step_a * density * kernel
        sums[1] += step_a * w2 * density * kernel
        sums[2] += step_a * normalized_score_numerator
        sums[3] += step_a * (w2 + 2) * normalized_score_numerator

        if min_kernel is None or kernel.lower() < min_kernel:
            min_kernel = kernel.lower()
        if kernel.upper() > max_kernel:
            max_kernel = kernel.upper()
        if derivative_tail.upper() > max_derivative_tail:
            max_derivative_tail = derivative_tail.upper()

    assert min_kernel is not None
    return sums, {
        "kernel_min_lower": min_kernel,
        "kernel_max_upper": max_kernel,
        "max_kernel_derivative_tail_upper": max_derivative_tail,
    }


def attach_normalized_score_tails(
    model: score.ScoreFullXiModel,
    central: list[arb],
    width: Fraction,
) -> tuple[list[arb], dict[str, arb]]:
    if width != 24:
        raise ValueError("the imported H1272 bounds require W=24")
    width_a = arb_from_fraction(width)
    tail_0 = arb("2e-40")
    tail_2 = arb("2e-37")
    tail_4 = arb("7e-35")
    boundary = model.density(width_a).upper() + model.density(-width_a).upper()
    score_a_tail = boundary + tail_2 / width_a
    score_c_tail = (width_a * width_a + 2) * boundary + tail_4 / width_a
    a_lower = model.a.lower()
    if not a_lower > 0:
        raise ValueError("the r-box does not certify a>0")

    normalized_a_tail = score_a_tail / a_lower
    normalized_c_tail = score_c_tail / a_lower
    totals = [
        central[0] + raw_base.nonnegative_error(tail_0),
        central[1] + raw_base.nonnegative_error(tail_2),
        central[2] + raw_base.symmetric_error(normalized_a_tail),
        central[3] + raw_base.symmetric_error(normalized_c_tail),
    ]
    return totals, {
        "a_lower": a_lower,
        "normalized_A_tail_radius": normalized_a_tail,
        "normalized_C_tail_radius": normalized_c_tail,
    }


def certificate(
    r_lo: Fraction,
    r_hi: Fraction,
    width: Fraction,
    subdivisions: int,
    m_cutoff: int,
    series_degree: int,
    digits: int,
) -> dict[str, Any]:
    dependency = json.loads(raw_base.H1272_JSON.read_text(encoding="utf-8"))
    model = score.ScoreFullXiModel(
        arb_box(r_lo, r_hi), m_cutoff, series_degree
    )
    central, diagnostics = central_normalized_score_sums(
        model, width, subdivisions
    )
    totals, tail_data = attach_normalized_score_tails(model, central, width)
    mass, moment2_num, normalized_a_num, normalized_c_num = totals
    second = moment2_num / mass
    normalized_a = normalized_a_num / mass
    normalized_c = normalized_c_num / mass
    third_over_a = (
        -normalized_c
        + 3 * normalized_a * second
        - 2 * model.a**2 * normalized_a**3
    )
    t = model.q / 2
    main_scale = 8 * t**2 * model.B / model.A**3
    scaled = main_scale * third_over_a
    margin = scaled + arb(3) / 4

    theta_tail, theta_gate = score.global_kernel_unit_bound()
    checks = {
        "H1272_dependency": bool(dependency.get("all_checks_pass")),
        "global_kernel_bound": bool(theta_tail.upper() < theta_gate.lower()),
        "mass_positive": bool(mass.lower() > 0),
        "a_positive": bool(model.a.lower() > 0),
        "kernel_positive_on_central_boxes": bool(
            diagnostics["kernel_min_lower"] > 0
        ),
        "kernel_below_one_on_central_boxes": bool(
            diagnostics["kernel_max_upper"] < 1
        ),
        "kernel_derivative_tail_small": bool(
            diagnostics["max_kernel_derivative_tail_upper"] < arb("1e-100")
        ),
        "finite_scaled_interval": (
            str(scaled.lower()) != "nan" and str(scaled.upper()) != "nan"
        ),
        "target_margin_positive": bool(margin.lower() > 0),
    }
    return {
        "id": "h1319_full_xi_normalized_score_rbox_pilot",
        "classification": "h1319_single_normalized_score_rbox_certificate_not_cover",
        "config": {
            "r_lo": str(r_lo),
            "r_hi": str(r_hi),
            "r_width": str(r_hi - r_lo),
            "W": str(width),
            "subdivisions": subdivisions,
            "m_cutoff": m_cutoff,
            "series_degree": series_degree,
            "finite_differences_used": False,
        },
        "mass": interval_payload(mass, digits),
        "second_moment": interval_payload(second, digits),
        "normalized_score_A": interval_payload(normalized_a, digits),
        "normalized_score_C": interval_payload(normalized_c, digits),
        "third_central_over_a": interval_payload(third_over_a, digits),
        "main_scale_8t2B_over_A3": interval_payload(main_scale, digits),
        "scaled_t2_kappa3": interval_payload(scaled, digits),
        "margin": interval_payload(margin, digits),
        "tail_data": {
            key: interval_payload(value, digits) for key, value in tail_data.items()
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "scope": "one exact rational r-box only; adaptive cover remains open",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r-lo", default="10")
    parser.add_argument("--r-hi", default="101/10")
    parser.add_argument("--W", default="24")
    parser.add_argument("--subdivisions", type=int, default=2000)
    parser.add_argument("--m-cutoff", type=int, default=8)
    parser.add_argument("--series-degree", type=int, default=48)
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--digits", type=int, default=24)
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()

    flint.ctx.dps = args.dps
    report = certificate(
        Fraction(args.r_lo),
        Fraction(args.r_hi),
        Fraction(args.W),
        args.subdivisions,
        args.m_cutoff,
        args.series_degree,
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
                "third_central_over_a": report["third_central_over_a"],
                "main_scale": report["main_scale_8t2B_over_A3"],
                "scaled_t2_kappa3": report["scaled_t2_kappa3"],
                "margin": report["margin"],
                "out": None if args.no_write else str(output),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0 if report["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
