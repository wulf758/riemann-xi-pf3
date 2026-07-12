"""Canonical correlated mean-value runner for the H1319 variance audit.

The first materialization deliberately exposed a failed componentwise quotient:
forming M2(R)/Z(R) from two independently widened ranges loses the r-correlation
and is far too broad.  H1319's own proof uses scalar mean-value enclosures, so
this runner applies the same construction directly to

    G= M2/Z,
    V= M2/Z-(SA/Z)^2.

Only this correlation repair and a cancellation-free A-2tr identity differ
from the base report generator.
"""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
from typing import Any

from flint import arb

import claude_h1319_box_variance_upper_certificate as engine


def abs_upper(value: arb) -> arb:
    return max(abs(value.lower()), abs(value.upper()))


def saddle_values(left: Fraction, right: Fraction) -> tuple[arb, arb, arb, arb]:
    r = engine.rball(left, right)
    alpha = arb(9) / 4
    p = arb.pi() * r.exp()
    q = p * r - alpha * r - 1
    t = q / 2
    big_a = r * (p * (r + 1) - alpha)
    # Exact cancellation-free form of A-2tr.
    gap = p * r + alpha * r * (r - 1) + r
    return r, t, big_a, gap


def scalar_mean_value_ranges(
    point: list[arb],
    derivative: list[arb],
    ranges: list[arb],
    errors: list[arb],
    radius: Fraction,
) -> tuple[arb, arb, dict[str, arb]]:
    z0, m20, sa0 = point[:3]
    z, m2, sa = ranges[:3]
    dz, dm2, dsa = derivative[:3]
    ez, em2, esa = [item.upper() for item in errors[:3]]
    offset = engine.symmetric(engine.qarb(radius))

    raw0 = m20 / z0
    raw_derivative = (dm2 * z - m2 * dz) / z**2
    finite_raw = raw0 + offset * raw_derivative

    var0 = m20 / z0 - (sa0 / z0) ** 2
    dvar_dz = -m2 / z**2 + 2 * sa**2 / z**3
    dvar_dm2 = 1 / z
    dvar_dsa = -2 * sa / z**2
    var_derivative = dvar_dz * dz + dvar_dm2 * dm2 + dvar_dsa * dsa
    finite_var = var0 + offset * var_derivative

    expanded_z = z + engine.symmetric(errors[0])
    expanded_m2 = m2 + engine.symmetric(errors[1])
    expanded_sa = sa + engine.symmetric(errors[2])
    if expanded_z.lower() <= 0:
        raise ValueError("expanded mass range contains zero")

    raw_tail_radius = (
        abs_upper(-expanded_m2 / expanded_z**2) * ez
        + abs_upper(1 / expanded_z) * em2
    )
    raw_full = finite_raw + engine.symmetric(raw_tail_radius)

    tail_dvar_dz = (
        -expanded_m2 / expanded_z**2
        + 2 * expanded_sa**2 / expanded_z**3
    )
    tail_dvar_dm2 = 1 / expanded_z
    tail_dvar_dsa = -2 * expanded_sa / expanded_z**2
    var_tail_radius = (
        abs_upper(tail_dvar_dz) * ez
        + abs_upper(tail_dvar_dm2) * em2
        + abs_upper(tail_dvar_dsa) * esa
    )
    var_full = finite_var + engine.symmetric(var_tail_radius)
    diagnostics = {
        "finite_raw_at_midpoint": raw0,
        "finite_raw_derivative": raw_derivative,
        "finite_raw_mean_value": finite_raw,
        "raw_tail_radius": raw_tail_radius,
        "finite_variance_at_midpoint": var0,
        "finite_variance_derivative": var_derivative,
        "finite_variance_mean_value": finite_var,
        "variance_tail_radius": var_tail_radius,
    }
    return raw_full, var_full, diagnostics


def middle_rows(path: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    report = json.loads(path.read_text(encoding="utf-8"))
    rows = []
    for stored in report.get("accepted_boxes", []):
        point = [
            engine.payload_ball(item)
            for item in stored["finite_integrals_at_midpoint"]
        ]
        derivatives = [
            engine.payload_ball(item)
            for item in stored["finite_integral_r_derivatives"]
        ]
        ranges = [
            engine.payload_ball(item)
            for item in stored["finite_integral_ranges_mean_value"]
        ]
        errors = [
            engine.payload_ball(item)
            for item in stored["uniform_integral_error_radii"]
        ]
        left = Fraction(stored["r_lo"])
        right = Fraction(stored["r_hi"])
        raw_second, variance, diagnostics = scalar_mean_value_ranges(
            point,
            derivatives,
            ranges,
            errors,
            (right - left) / 2,
        )
        expanded_mass = ranges[0] + engine.symmetric(errors[0])
        expanded_second = ranges[1] + engine.symmetric(errors[1])
        expanded_score = ranges[2] + engine.symmetric(errors[2])
        r, t, big_a, a_gap = saddle_values(left, right)
        kappa2_upper = 4 * raw_second / big_a
        tr_kappa2_upper = t * r * kappa2_upper
        row = {
            "source": str(path).replace("\\", "/"),
            "kind": "score_scalar_mean_value_variance",
            "r_lo": str(left),
            "r_hi": str(right),
            "center_description": (
                "G=M2/Z and V=M2/Z-(SA/Z)^2 use scalar mean-value forms; "
                "SA/Z=E[N_Xi]=-E[W]"
            ),
            "mass": expanded_mass,
            "first_or_score": expanded_score,
            "second": expanded_second,
            "raw_second": raw_second,
            "centered_variance": variance,
            "r": r,
            "t": t,
            "A": big_a,
            "A_minus_2tr": a_gap,
            "kappa2_upper": kappa2_upper,
            "tr_kappa2_upper": tr_kappa2_upper,
            "variance_margin": arb(26) / 25 - raw_second,
            "tr_margin": arb(52) / 25 - tr_kappa2_upper,
            "mean_value_diagnostics": diagnostics,
        }
        rows.append(row)

    configuration = report.get("configuration", {})
    global_checks = report.get("global_checks", {})
    checks = {
        "classification": report.get("classification")
        == "h1319_full_xi_compact_mean_value_r_interval_certificate",
        "all_checks_pass": report.get("all_checks_pass") is True,
        "unresolved_boxes_empty": report.get("unresolved_boxes") == [],
        "accepted_boxes_nonempty": bool(rows),
        "global_checks_all_true": bool(global_checks) and all(global_checks.values()),
        "mean_value_not_sampling": configuration.get("point_sampling_used") is False,
        "no_finite_differences": configuration.get("finite_differences_used") is False,
        "stored_midpoints_present": all(
            "finite_integrals_at_midpoint" in row
            for row in report.get("accepted_boxes", [])
        ),
        "stored_derivatives_present": all(
            "finite_integral_r_derivatives" in row
            for row in report.get("accepted_boxes", [])
        ),
    }
    return rows, {
        "path": str(path).replace("\\", "/"),
        "sha256": engine.sha256(path),
        "covered_r_interval": report.get("covered_r_interval"),
        "accepted_box_count": len(rows),
        "checks": checks,
        "passes": all(checks.values()),
    }


def mutation_canary(worst: dict[str, Any], factor: Fraction) -> dict[str, Any]:
    mutated_raw = worst["raw_second"] * engine.qarb(factor)
    margin = arb(26) / 25 - mutated_raw
    return {
        "target_box": [worst["r_lo"], worst["r_hi"]],
        "target_source": worst["source"],
        "mutated_field": "certified_raw_second_envelope",
        "factor": str(factor),
        "mutated_raw_second": engine.ball_record(mutated_raw),
        "mutated_margin": engine.ball_record(margin),
        "kills_26_over_25_certificate": bool(margin.upper() < 0),
    }


engine.saddle_values = saddle_values
engine.middle_rows = middle_rows
engine.mutation_canary = mutation_canary


if __name__ == "__main__":
    raise SystemExit(engine.main())
