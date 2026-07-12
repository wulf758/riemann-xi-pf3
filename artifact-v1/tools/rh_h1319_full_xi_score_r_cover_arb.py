#!/usr/bin/env python3
"""Adaptive Arb cover for the compact H1319 full-Xi saddle interval.

The covered interval is [647/200, 71].  Its left endpoint lies strictly below
the saddle r solving t(r)=125, so a successful cover proves the desired
compact statement without a numerical root dependency.

Each r-box uses the cancellation-aware score integrals from
``rh_h1319_full_xi_fixed_r_score_arb``.  Both r and W are Arb intervals; no
point sampling or finite differencing is used.
"""

from __future__ import annotations

import argparse
import json
from collections import deque
from fractions import Fraction
from pathlib import Path
from typing import Any

import flint
from flint import arb

import rh_h1319_full_xi_fixed_r_arb_pilot as raw_base
import rh_h1319_full_xi_fixed_r_score_arb as score
from rh_h930_fixed_r_arb_tail_certificate import arb_box, interval_payload


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = (
    ROOT / "research" / "riemann" / "h1319_full_xi_score_r_cover_arb.json"
)


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def box_certificate(
    r_lo: Fraction,
    r_hi: Fraction,
    width: Fraction,
    subdivisions: int,
    m_cutoff: int,
    series_degree: int,
    digits: int,
) -> dict[str, Any]:
    """Certify one closed rational r-box."""

    model = score.ScoreFullXiModel(
        arb_box(r_lo, r_hi), m_cutoff, series_degree
    )
    central, diagnostics = score.central_score_sums(model, width, subdivisions)
    totals, tail_data = score.attach_score_tails(model, central, width)
    mass, moment2_num, score_a_num, score_c_num = totals
    second = moment2_num / mass
    score_a = score_a_num / mass
    score_c = score_c_num / mass
    third = -score_c + 3 * score_a * second - 2 * score_a**3
    t = model.q / 2
    scaled = t**2 * 8 * third / (model.A * model.sqrt_A)
    margin = scaled + arb(3) / 4

    checks = {
        "mass_positive": bool(mass.lower() > 0),
        "kernel_positive_on_central_boxes": bool(
            diagnostics["kernel_min_lower"] > 0
        ),
        "kernel_below_one_on_central_boxes": bool(
            diagnostics["kernel_max_upper"] < 1
        ),
        "kernel_m_tail_small": bool(model.kernel_tail.upper() < arb("1e-100")),
        "kernel_derivative_tail_small": bool(
            diagnostics["max_kernel_derivative_tail_upper"] < arb("1e-100")
        ),
        "finite_scaled_interval": (
            str(scaled.lower()) != "nan" and str(scaled.upper()) != "nan"
        ),
        "target_minus_three_over_four": bool(margin.lower() > 0),
    }
    return {
        "r_lo": fraction_text(r_lo),
        "r_hi": fraction_text(r_hi),
        "r_width": fraction_text(r_hi - r_lo),
        "subdivisions_W": subdivisions,
        "scaled_t2_kappa_Xi_third": interval_payload(scaled, digits),
        "margin_over_minus_3_over_4": interval_payload(margin, digits),
        "mass": interval_payload(mass, digits),
        "score_A": interval_payload(score_a, digits),
        "score_C": interval_payload(score_c, digits),
        "second_moment": interval_payload(second, digits),
        "third_central_moment_W": interval_payload(third, digits),
        "kernel_min_lower": interval_payload(
            diagnostics["kernel_min_lower"], digits
        ),
        "kernel_max_upper": interval_payload(
            diagnostics["kernel_max_upper"], digits
        ),
        "max_kernel_derivative_tail_upper": interval_payload(
            diagnostics["max_kernel_derivative_tail_upper"], digits
        ),
        "score_tail_radii": {
            "A": interval_payload(tail_data["score_A_tail_radius"], digits),
            "C": interval_payload(tail_data["score_C_tail_radius"], digits),
        },
        "checks": checks,
        "passes": all(checks.values()),
    }


def seed_boxes(
    left: Fraction, right: Fraction, seed_width: Fraction
) -> list[tuple[Fraction, Fraction]]:
    boxes = []
    cursor = left
    while cursor < right:
        endpoint = min(cursor + seed_width, right)
        boxes.append((cursor, endpoint))
        cursor = endpoint
    return boxes


def adaptive_cover(args: argparse.Namespace) -> dict[str, Any]:
    dependency = json.loads(raw_base.H1272_JSON.read_text(encoding="utf-8"))
    theta_tail, theta_gate = score.global_kernel_unit_bound()
    r_min = Fraction(args.r_min)
    r_max = Fraction(args.r_max)
    seed_width = Fraction(args.seed_width)
    min_r_width = Fraction(args.min_r_width)
    width = Fraction(args.W)
    if not r_min < r_max:
        raise ValueError("r_min must be strictly smaller than r_max")
    if seed_width <= 0 or min_r_width <= 0:
        raise ValueError("box widths must be positive")

    pending = deque(
        (lo, hi, 0, args.subdivisions)
        for lo, hi in seed_boxes(r_min, r_max, seed_width)
    )
    accepted: list[dict[str, Any]] = []
    unresolved: list[dict[str, Any]] = []
    evaluations = 0

    while pending:
        r_lo, r_hi, depth, subdivisions = pending.popleft()
        result = box_certificate(
            r_lo,
            r_hi,
            width,
            subdivisions,
            args.m_cutoff,
            args.series_degree,
            args.digits,
        )
        result["depth"] = depth
        evaluations += 1
        if result["passes"]:
            accepted.append(result)
            continue

        current_width = r_hi - r_lo
        if current_width > min_r_width and depth < args.max_r_depth:
            midpoint = (r_lo + r_hi) / 2
            pending.appendleft((midpoint, r_hi, depth + 1, subdivisions))
            pending.appendleft((r_lo, midpoint, depth + 1, subdivisions))
            continue

        if subdivisions < args.max_subdivisions:
            pending.appendleft(
                (
                    r_lo,
                    r_hi,
                    depth,
                    min(2 * subdivisions, args.max_subdivisions),
                )
            )
            continue

        unresolved.append(result)
        if len(unresolved) >= args.max_unresolved:
            break

    accepted.sort(key=lambda row: Fraction(row["r_lo"]))
    unresolved.sort(key=lambda row: Fraction(row["r_lo"]))
    exact_chain = bool(accepted) and accepted[0]["r_lo"] == fraction_text(r_min)
    if exact_chain:
        for left, right in zip(accepted, accepted[1:]):
            if left["r_hi"] != right["r_lo"]:
                exact_chain = False
                break
        exact_chain = exact_chain and accepted[-1]["r_hi"] == fraction_text(r_max)

    worst = None
    if accepted:
        worst = min(
            accepted,
            key=lambda row: arb(row["margin_over_minus_3_over_4"]["lower"]),
        )

    global_checks = {
        "H1272_dependency": bool(dependency.get("all_checks_pass")),
        "global_kernel_zero_lt_H_lt_one": bool(
            theta_tail.upper() < theta_gate.lower()
        ),
        "left_endpoint_below_r_of_t_125": bool(
            # t(647/200)=124.9637... <125, certified directly in Arb.
            score.ScoreFullXiModel(arb_box(r_min, r_min), 1, 8).q.upper()
            < arb(250)
        ),
        "accepted_boxes_nonempty": bool(accepted),
        "no_unresolved_boxes": not unresolved and not pending,
        "exact_gap_free_r_chain": exact_chain,
        "every_accepted_box_passes": all(row["passes"] for row in accepted),
    }
    return {
        "id": "h1319_full_xi_score_r_cover_arb",
        "classification": "h1319_full_xi_compact_r_interval_certificate",
        "statement": "t(r)^2*kappa_Xi'''(t(r)) >= -3/4",
        "covered_r_interval": [fraction_text(r_min), fraction_text(r_max)],
        "configuration": {
            "W": fraction_text(width),
            "initial_subdivisions_W": args.subdivisions,
            "max_subdivisions_W": args.max_subdivisions,
            "m_cutoff": args.m_cutoff,
            "series_degree": args.series_degree,
            "seed_width": fraction_text(seed_width),
            "min_r_width": fraction_text(min_r_width),
            "max_r_depth": args.max_r_depth,
            "flint_dps": flint.ctx.dps,
            "finite_differences_used": False,
            "point_sampling_used": False,
        },
        "summary": {
            "evaluations": evaluations,
            "accepted_box_count": len(accepted),
            "unresolved_box_count": len(unresolved),
            "smallest_accepted_r_width": (
                fraction_text(
                    min(Fraction(row["r_width"]) for row in accepted)
                )
                if accepted
                else None
            ),
            "largest_accepted_r_width": (
                fraction_text(
                    max(Fraction(row["r_width"]) for row in accepted)
                )
                if accepted
                else None
            ),
            "worst_margin_box": worst,
        },
        "global_kernel_bound": {
            "theta_tail_at_pi": interval_payload(theta_tail, args.digits),
            "three_over_two_pi": interval_payload(theta_gate, args.digits),
        },
        "global_checks": global_checks,
        "accepted_boxes": accepted,
        "unresolved_boxes": unresolved,
        "all_checks_pass": all(global_checks.values()),
        "scope": (
            "compact full-Xi curvature certificate only; H920 piecewise "
            "assembly and finite n<=126 are separate dependencies"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r-min", default="647/200")
    parser.add_argument("--r-max", default="71")
    parser.add_argument("--seed-width", default="1/2")
    parser.add_argument("--min-r-width", default="1/1024")
    parser.add_argument("--max-r-depth", type=int, default=12)
    parser.add_argument("--W", default="24")
    parser.add_argument("--subdivisions", type=int, default=2000)
    parser.add_argument("--max-subdivisions", type=int, default=16000)
    parser.add_argument("--m-cutoff", type=int, default=8)
    parser.add_argument("--series-degree", type=int, default=48)
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--digits", type=int, default=20)
    parser.add_argument("--max-unresolved", type=int, default=8)
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    flint.ctx.dps = args.dps
    report = adaptive_cover(args)
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "id": report["id"],
                "summary": report["summary"],
                "global_checks": report["global_checks"],
                "all_checks_pass": report["all_checks_pass"],
                "out": str(output),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0 if report["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
