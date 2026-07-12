#!/usr/bin/env python3
"""Stable-scale variant of the normalized-score H1319 r-box pilot."""

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
import rh_h1319_full_xi_normalized_score_rbox_pilot as v1
from rh_h930_fixed_r_arb_tail_certificate import arb_box, interval_payload


DEFAULT_OUT = Path(
    "research/riemann/h1319_full_xi_normalized_score_rbox_v2.json"
)


def stable_main_scale(model: score.ScoreFullXiModel) -> arb:
    """Return the exact 8*t^2*B/A^3 after cancelling pi*r*exp(r)."""

    d = model.alpha / (model.pi * model.exp_r)
    q_ratio = 1 - d - 1 / (model.pi * model.r * model.exp_r)
    return (
        2
        * q_ratio**2
        * (model.r**2 + 3 * model.r + 1 - d)
        / (model.r + 1 - d) ** 3
    )


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
    central, diagnostics = v1.central_normalized_score_sums(
        model, width, subdivisions
    )
    totals, tail_data = v1.attach_normalized_score_tails(
        model, central, width
    )
    mass, moment2_num, normalized_a_num, normalized_c_num = totals
    second = moment2_num / mass
    normalized_a = normalized_a_num / mass
    normalized_c = normalized_c_num / mass
    third_over_a = (
        -normalized_c
        + 3 * normalized_a * second
        - 2 * model.a**2 * normalized_a**3
    )
    main_scale = stable_main_scale(model)
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
        "finite_scaled_interval": (
            str(scaled.lower()) != "nan" and str(scaled.upper()) != "nan"
        ),
        "target_margin_positive": bool(margin.lower() > 0),
    }
    return {
        "id": "h1319_full_xi_normalized_score_rbox_v2",
        "classification": "h1319_single_normalized_score_rbox_v2_not_cover",
        "config": {
            "r_lo": str(r_lo),
            "r_hi": str(r_hi),
            "r_width": str(r_hi - r_lo),
            "W": str(width),
            "subdivisions": subdivisions,
            "m_cutoff": m_cutoff,
            "series_degree": series_degree,
            "main_scale_form": "common-factor-cancelled",
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
