#!/usr/bin/env python3
"""Independent audit checks for the H1319 score mean-value certificates.

This script does not generate a curvature certificate.  It checks the three
already materialized reports, re-evaluates their worst boxes, and compares the
forward-mode Arb derivatives against an independently written mpmath model.
The symbolic checks cover the stable score identity, the kernel w-derivative,
and the perturbation gradient.
"""

from __future__ import annotations

import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable

import flint
import mpmath as mp
import sympy as sp
from flint import arb


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import rh_h1319_full_xi_score_mean_value_r_cover_arb as engine


REPORTS = (
    ROOT
    / "research"
    / "riemann"
    / "h1319_full_xi_score_mean_value_r_cover_4_5.json",
    ROOT
    / "research"
    / "riemann"
    / "h1319_full_xi_score_mean_value_r_cover_5_10.json",
    ROOT
    / "research"
    / "riemann"
    / "h1319_full_xi_score_mean_value_r10_to_22_arb.json",
)
H1272 = (
    ROOT
    / "research"
    / "riemann"
    / "h1272_gamma_log_even_moment_tail_global.json"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def payload_ball(payload: dict[str, str]) -> arb:
    return arb(payload["interval"])


def endpoint_lower(payload: dict[str, str]) -> arb:
    return arb(payload["lower"]).lower()


def endpoint_upper(payload: dict[str, str]) -> arb:
    return arb(payload["upper"]).upper()


def symbolic_checks() -> dict[str, bool]:
    r, h, p, alpha, delta, a = sp.symbols(
        "r h p alpha delta a", positive=True
    )
    x, m, d, xp = sp.symbols("x m d xp", positive=True)
    z, m2, sa, sc = sp.symbols("z m2 sa sc", nonzero=True)

    phi = lambda value: sp.exp(value) - 1 - value
    delta_def = r * (sp.exp(h) - 1)
    big_a = r * (p * (r + 1) - alpha)
    direct_score = (
        p * r * sp.exp(h) * (sp.exp(delta_def) - 1)
        + (p * r - alpha * r) * (sp.exp(h) - 1)
        - big_a * h
    )
    stable_score = big_a * phi(h) + p * (
        r * phi(delta_def)
        + delta_def * (sp.exp(delta_def) - 1)
    )

    coefficient = m**4 - 3 * m**2 / (2 * x)
    kernel_term = coefficient * sp.exp(-d * x)
    kernel_expected = xp * sp.exp(-d * x) * (
        3 * m**2 / (2 * x**2) - d * coefficient
    )
    kernel_actual = sp.diff(kernel_term, x) * xp

    cumulant = -sc / z + 3 * sa * m2 / z**2 - 2 * sa**3 / z**3
    gradient_expected = (
        sc / z**2 - 6 * sa * m2 / z**3 + 6 * sa**3 / z**4,
        3 * sa / z**2,
        3 * m2 / z**2 - 6 * sa**2 / z**3,
        -1 / z,
    )
    gradient_actual = tuple(
        sp.diff(cumulant, variable) for variable in (z, m2, sa, sc)
    )

    return {
        "phi2_derivative": sp.simplify(
            sp.diff(phi(a), a) - (sp.exp(a) - 1)
        )
        == 0,
        "stable_dominant_score_identity": sp.simplify(
            direct_score - stable_score
        )
        == 0,
        "finite_kernel_w_derivative": sp.simplify(
            kernel_actual - kernel_expected
        )
        == 0,
        "cumulant_perturbation_gradient": all(
            sp.simplify(left - right) == 0
            for left, right in zip(gradient_actual, gradient_expected)
        ),
    }


def exact_fields(r: mp.mpf, w: mp.mpf, m_cutoff: int = 8) -> list[mp.mpf]:
    pi = mp.pi
    alpha = mp.mpf(9) / 4
    exp_r = mp.exp(r)
    p = pi * exp_r
    q = p * r - alpha * r - 1
    t = q / 2
    big_a = r * (p * (r + 1) - alpha)
    sqrt_a = mp.sqrt(big_a)
    h = w / sqrt_a
    delta = r * mp.expm1(h)
    phi2 = lambda value: mp.expm1(value) - value
    psi = p * phi2(delta) + (p * r - alpha * r) * phi2(h)
    density = mp.exp(-psi)
    n_dom = (
        big_a * phi2(h)
        + p * (r * phi2(delta) + delta * mp.expm1(delta))
    ) / sqrt_a
    rho = r * mp.exp(h)
    x = pi * mp.exp(rho)
    x_prime = x * rho / sqrt_a
    kernel = mp.mpf(0)
    kernel_prime = mp.mpf(0)
    for m_int in range(1, m_cutoff + 1):
        m = mp.mpf(m_int)
        d = m * m - 1
        exponential = mp.exp(-x * d)
        coefficient = m**4 - 3 * m**2 / (2 * x)
        kernel += coefficient * exponential
        kernel_prime += x_prime * exponential * (
            3 * m**2 / (2 * x**2) - d * coefficient
        )
    score_numerator = density * (kernel * n_dom - kernel_prime)
    mass_density = density * kernel
    w2 = w * w
    scale = 8 * t * t / (big_a * sqrt_a)
    return [
        psi,
        density,
        n_dom,
        kernel,
        kernel_prime,
        mass_density,
        w2 * mass_density,
        score_numerator,
        (w2 + 2) * score_numerator,
        scale,
    ]


def arb_fields(model: engine.DualScoreModel, w_value: arb) -> list[engine.Dual]:
    kernel, kernel_prime = model.finite_kernel_and_w_derivative(w_value)
    integrands = model.finite_integrands(w_value)
    return [
        model.psi(w_value),
        model.density(w_value),
        model.dominant_score_excess(w_value),
        kernel,
        kernel_prime,
        *integrands,
        model.scale(),
    ]


def mp_ball(value: mp.mpf) -> arb:
    return arb(mp.nstr(value, 80))


def dual_probe_checks() -> dict[str, Any]:
    mp.mp.dps = 100
    probes = (
        (mp.mpf("4.025"), mp.mpf("-24")),
        (mp.mpf("4.025"), mp.mpf("0")),
        (mp.mpf("4.025"), mp.mpf("24")),
        (mp.mpf("10.025"), mp.mpf("-24")),
        (mp.mpf("10.025"), mp.mpf("5")),
        (mp.mpf("10.025"), mp.mpf("24")),
        (mp.mpf("21.975"), mp.mpf("-24")),
        (mp.mpf("21.975"), mp.mpf("0")),
        (mp.mpf("21.975"), mp.mpf("24")),
    )
    names = (
        "psi",
        "density",
        "dominant_score_excess",
        "finite_kernel",
        "finite_kernel_w_derivative",
        "integrand_Z",
        "integrand_M2",
        "integrand_SA",
        "integrand_SC",
        "scale_8t2_over_A32",
    )
    failures: list[dict[str, str]] = []
    comparisons = 0
    for r_value, w_value in probes:
        model = engine.DualScoreModel(
            arb(mp.nstr(r_value, 30)), arb(1), 8, 48
        )
        balls = arb_fields(model, arb(mp.nstr(w_value, 30)))
        exact = exact_fields(r_value, w_value)
        exact_derivatives = [
            mp.diff(
                lambda variable, index=index: exact_fields(
                    variable, w_value
                )[index],
                r_value,
            )
            for index in range(len(names))
        ]
        for name, ball, value, derivative in zip(
            names, balls, exact, exact_derivatives
        ):
            comparisons += 2
            if not ball.val.overlaps(mp_ball(value)):
                failures.append(
                    {
                        "r": mp.nstr(r_value),
                        "w": mp.nstr(w_value),
                        "field": name,
                        "component": "value",
                    }
                )
            if not ball.der.overlaps(mp_ball(derivative)):
                failures.append(
                    {
                        "r": mp.nstr(r_value),
                        "w": mp.nstr(w_value),
                        "field": name,
                        "component": "r_derivative",
                    }
                )
    return {
        "probe_count": len(probes),
        "field_count_per_probe": len(names),
        "comparison_count": comparisons,
        "failures": failures,
        "all_pass": not failures,
    }


def report_checks(path: Path) -> dict[str, Any]:
    report = json.loads(path.read_text(encoding="utf-8"))
    left = Fraction(report["covered_r_interval"][0])
    right = Fraction(report["covered_r_interval"][1])
    boxes = report["accepted_boxes"]
    cursor = left
    chain = True
    rows_pass = True
    margin_positive = True
    full_mass_positive = True
    errors_nonnegative = True
    ordered = True
    previous_lo: Fraction | None = None
    for row in boxes:
        row_lo = Fraction(row["r_lo"])
        row_hi = Fraction(row["r_hi"])
        chain = chain and row_lo == cursor and row_hi > row_lo
        cursor = row_hi
        if previous_lo is not None:
            ordered = ordered and previous_lo < row_lo
        previous_lo = row_lo
        rows_pass = rows_pass and row["passes"] and all(
            row["checks"].values()
        )
        margin_positive = margin_positive and bool(
            endpoint_lower(row["margin_over_minus_3_over_4"]) > 0
        )
        finite_mass_lower = endpoint_lower(
            row["finite_integral_ranges_mean_value"][0]
        )
        mass_error_upper = endpoint_upper(
            row["uniform_integral_error_radii"][0]
        )
        full_mass_positive = full_mass_positive and bool(
            finite_mass_lower - mass_error_upper > 0
        )
        errors_nonnegative = errors_nonnegative and all(
            endpoint_lower(payload) >= 0
            for payload in row["uniform_integral_error_radii"]
        )
    chain = chain and cursor == right

    config = report["configuration"]
    worst = report["summary"]["worst_margin_box"]
    recomputed = engine.box_certificate(
        Fraction(worst["r_lo"]),
        Fraction(worst["r_hi"]),
        Fraction(config["W"]),
        int(config["subdivisions_W"]),
        int(config["m_cutoff"]),
        int(config["series_degree"]),
        20,
    )
    recomputed_overlap = payload_ball(
        worst["scaled_t2_kappa_Xi_third"]
    ).overlaps(payload_ball(recomputed["scaled_t2_kappa_Xi_third"]))

    checks = {
        "declared_all_checks_pass": bool(report["all_checks_pass"]),
        "global_checks_all_true": all(report["global_checks"].values()),
        "covered_range_in_H1272_domain": left >= Fraction(5, 2),
        "W_equals_24": Fraction(config["W"]) == 24,
        "no_finite_differences": not config["finite_differences_used"],
        "no_point_sampling": not config["point_sampling_used"],
        "accepted_count_matches": len(boxes)
        == report["summary"]["accepted_box_count"],
        "no_unresolved": not report["unresolved_boxes"],
        "rows_ordered": ordered,
        "independent_exact_fraction_chain": chain,
        "every_row_passes": rows_pass,
        "every_margin_strictly_positive": margin_positive,
        "every_full_mass_strictly_positive": full_mass_positive,
        "every_uniform_error_nonnegative": errors_nonnegative,
        "worst_box_recomputation_overlaps": bool(recomputed_overlap),
        "worst_box_recomputation_passes": bool(recomputed["passes"]),
    }
    return {
        "path": str(path.relative_to(ROOT)).replace("\\", "/"),
        "sha256": sha256(path),
        "covered_r_interval": [str(left), str(right)],
        "accepted_box_count": len(boxes),
        "worst_box": [worst["r_lo"], worst["r_hi"]],
        "worst_scaled_lower": worst["scaled_t2_kappa_Xi_third"]["lower"],
        "worst_margin_lower": worst["margin_over_minus_3_over_4"]["lower"],
        "checks": checks,
        "all_pass": all(checks.values()),
    }


def h1272_checks() -> dict[str, Any]:
    report = json.loads(H1272.read_text(encoding="utf-8"))
    bounds = report["global_bounds"]
    theta_tail, theta_gate = engine.score_base.global_kernel_unit_bound()
    checks = {
        "dependency_all_checks_pass": bool(report["all_checks_pass"]),
        "declared_range_r_ge_5_over_2": "r>=5/2" in report["statement"],
        "tail_0_below_imported_target": mp.mpf(
            bounds["0"]["global_upper"]
        )
        < mp.mpf("2e-40"),
        "tail_2_below_imported_target": mp.mpf(
            bounds["2"]["global_upper"]
        )
        < mp.mpf("2e-37"),
        "tail_4_below_imported_target": mp.mpf(
            bounds["4"]["global_upper"]
        )
        < mp.mpf("7e-35"),
        "global_kernel_theta_tail_gate": bool(
            theta_tail.upper() < theta_gate.lower()
        ),
    }
    return {
        "path": str(H1272.relative_to(ROOT)).replace("\\", "/"),
        "sha256": sha256(H1272),
        "checks": checks,
        "all_pass": all(checks.values()),
    }


def main() -> int:
    flint.ctx.dps = 90
    symbolic = symbolic_checks()
    reports = [report_checks(path) for path in REPORTS]
    dual = dual_probe_checks()
    dependency = h1272_checks()
    engine_path = TOOLS / "rh_h1319_full_xi_score_mean_value_r_cover_arb.py"
    result = {
        "id": "h1319_score_mean_value_independent_audit",
        "classification": "independent_audit_not_a_new_curvature_certificate",
        "engine": {
            "path": str(engine_path.relative_to(ROOT)).replace("\\", "/"),
            "sha256": sha256(engine_path),
        },
        "symbolic_identity_checks": symbolic,
        "dual_against_independent_mpmath": dual,
        "H1272_dependency": dependency,
        "reports": reports,
        "known_nonblocking_limitation": (
            "The generic CLI does not reject r_min<5/2 even though H1272 "
            "starts at r=5/2; all audited reports start at r>=4."
        ),
    }
    result["all_checks_pass"] = (
        all(symbolic.values())
        and dual["all_pass"]
        and dependency["all_pass"]
        and all(report["all_pass"] for report in reports)
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
