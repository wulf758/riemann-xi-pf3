#!/usr/bin/env python3
"""Mean-value Arb r-cover for the compact H1319 full-Xi curvature bound.

The direct interval-r score engine is stable at a fixed r but loses the
correlation between the saddle scale and the score integrals on an r-box.
This engine keeps that correlation by using the rigorous mean-value form

    F(r) in F(r0) + (r-r0) F'([r_lo,r_hi]),

where

    F(r) = 8 t(r)^2 A(r)^(-3/2) kappa_3(W_r).

The four finite-m score integrals and their r derivatives are enclosed by
Arb box quadrature.  Derivatives are obtained by forward automatic
differentiation, not finite differences.  The omitted m-tail and the H1272
outer tails are uniform perturbations; their effect on kappa_3 is bounded by
the explicit gradient in (Z,M2,SA,SC).

A successful gap-free chain proves F(r) >= -3/4 on its configured range.
Failed leaves are interval-method failures, never counterexamples.
"""

from __future__ import annotations

import argparse
import json
from collections import deque
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any

import flint
from flint import arb

import rh_h1319_full_xi_fixed_r_arb_pilot as raw_base
import rh_h1319_full_xi_fixed_r_score_arb as score_base
from rh_h1319_full_xi_fixed_r_series_arb import expm1_minus_x_series
from rh_h930_fixed_r_arb_tail_certificate import (
    arb_box,
    arb_from_fraction,
    interval_payload,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = (
    ROOT
    / "research"
    / "riemann"
    / "h1319_full_xi_score_mean_value_r_cover_arb.json"
)


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def abs_upper(value: arb) -> arb:
    return max(abs(value.lower()), abs(value.upper()))


def symmetric(radius: arb) -> arb:
    return arb(0, 1) * radius.upper()


def as_arb(value: Any) -> arb:
    if isinstance(value, arb):
        return value
    if isinstance(value, Fraction):
        return arb_from_fraction(value)
    return arb(value)


@dataclass(frozen=True)
class Dual:
    """First-order forward-mode dual number with Arb coefficients."""

    val: arb
    der: arb

    @staticmethod
    def coerce(value: Any) -> "Dual":
        if isinstance(value, Dual):
            return value
        return Dual(as_arb(value), arb(0))

    def __add__(self, other: Any) -> "Dual":
        rhs = Dual.coerce(other)
        return Dual(self.val + rhs.val, self.der + rhs.der)

    def __radd__(self, other: Any) -> "Dual":
        return self + other

    def __sub__(self, other: Any) -> "Dual":
        rhs = Dual.coerce(other)
        return Dual(self.val - rhs.val, self.der - rhs.der)

    def __rsub__(self, other: Any) -> "Dual":
        return Dual.coerce(other) - self

    def __neg__(self) -> "Dual":
        return Dual(-self.val, -self.der)

    def __mul__(self, other: Any) -> "Dual":
        rhs = Dual.coerce(other)
        return Dual(
            self.val * rhs.val,
            self.der * rhs.val + self.val * rhs.der,
        )

    def __rmul__(self, other: Any) -> "Dual":
        return self * other

    def reciprocal(self) -> "Dual":
        return Dual(1 / self.val, -self.der / (self.val * self.val))

    def __truediv__(self, other: Any) -> "Dual":
        return self * Dual.coerce(other).reciprocal()

    def __rtruediv__(self, other: Any) -> "Dual":
        return Dual.coerce(other) / self

    def __pow__(self, exponent: int) -> "Dual":
        if not isinstance(exponent, int):
            raise TypeError("Dual powers must be integral")
        if exponent == 0:
            return Dual(arb(1), arb(0))
        if exponent < 0:
            return (self.reciprocal()) ** (-exponent)
        return Dual(
            self.val**exponent,
            arb(exponent) * self.val ** (exponent - 1) * self.der,
        )

    def exp(self) -> "Dual":
        exponential = self.val.exp()
        return Dual(exponential, exponential * self.der)

    def expm1(self) -> "Dual":
        return Dual(self.val.expm1(), self.val.exp() * self.der)

    def sqrt(self) -> "Dual":
        root = self.val.sqrt()
        return Dual(root, self.der / (2 * root))


def phi2(value: Dual, degree: int) -> Dual:
    """Cancellation-free phi2(x)=exp(x)-1-x and its exact derivative."""

    return Dual(
        expm1_minus_x_series(value.val, degree),
        value.val.expm1() * value.der,
    )


class DualScoreModel:
    """Finite-m full-Xi score model, differentiated with respect to r."""

    def __init__(
        self,
        r_value: arb,
        r_derivative: arb,
        m_cutoff: int,
        series_degree: int,
    ) -> None:
        self.pi = arb.pi()
        self.alpha = arb(9) / 4
        self.r = Dual(r_value, r_derivative)
        self.exp_r = self.r.exp()
        self.q = self.pi * self.r * self.exp_r - self.alpha * self.r - 1
        self.t = self.q / 2
        self.A = self.r * (
            self.pi * self.exp_r * (self.r + 1) - self.alpha
        )
        self.sqrt_A = self.A.sqrt()
        self.m_cutoff = m_cutoff
        self.series_degree = series_degree

    def psi(self, w: arb) -> Dual:
        h = Dual(w, arb(0)) / self.sqrt_A
        delta = self.r * h.expm1()
        return (
            self.pi * self.exp_r * phi2(delta, self.series_degree)
            + (self.pi * self.r * self.exp_r - self.alpha * self.r)
            * phi2(h, self.series_degree)
        )

    def density(self, w: arb) -> Dual:
        return (-self.psi(w)).exp()

    def dominant_score_excess(self, w: arb) -> Dual:
        h = Dual(w, arb(0)) / self.sqrt_A
        delta = self.r * h.expm1()
        return (
            self.A * phi2(h, self.series_degree)
            + self.pi
            * self.exp_r
            * (
                self.r * phi2(delta, self.series_degree)
                + delta * delta.expm1()
            )
        ) / self.sqrt_A

    def finite_kernel_and_w_derivative(self, w: arb) -> tuple[Dual, Dual]:
        h = Dual(w, arb(0)) / self.sqrt_A
        rho = self.r * h.exp()
        x = self.pi * rho.exp()
        x_prime = x * rho / self.sqrt_A
        value = Dual(arb(0), arb(0))
        w_derivative = Dual(arb(0), arb(0))
        for m_int in range(1, self.m_cutoff + 1):
            m = arb(m_int)
            d = arb(m_int * m_int - 1)
            exponential = (-x * d).exp()
            coefficient = m**4 - arb(3) * m**2 / (2 * x)
            value += coefficient * exponential
            w_derivative += x_prime * exponential * (
                arb(3) * m**2 / (2 * x * x) - d * coefficient
            )
        return value, w_derivative

    def finite_integrands(self, w: arb) -> list[Dual]:
        density = self.density(w)
        kernel, kernel_prime = self.finite_kernel_and_w_derivative(w)
        n_dom = self.dominant_score_excess(w)
        score_numerator = density * (kernel * n_dom - kernel_prime)
        w2 = w * w
        mass_density = density * kernel
        return [
            mass_density,
            w2 * mass_density,
            score_numerator,
            (w2 + 2) * score_numerator,
        ]

    def scale(self) -> Dual:
        return 8 * self.t * self.t / (self.A * self.sqrt_A)


def kernel_derivative_tail_upper(
    model: DualScoreModel, w: arb
) -> arb:
    """Uniform |dH_tail/dw| bound copied from the fixed-r score proof."""

    rho = model.r.val * (w / model.sqrt_A.val).exp()
    x = arb.pi() * rho.exp()
    x_lo = x.lower()
    scale_hi = (x * rho / model.sqrt_A.val).upper()
    first_m = model.m_cutoff + 1
    m = arb(first_m)
    coefficient = 1 + arb(3) / (2 * x_lo) + arb(3) / (2 * x_lo * x_lo)
    first = coefficient * m**6 * (-x_lo * (first_m * first_m - 1)).exp()
    ratio = (
        arb(first_m + 1) / arb(first_m)
    ) ** 6 * (-x_lo * (2 * first_m + 1)).exp()
    return scale_hi * first / (1 - ratio)


def central_point_and_derivative_sums(
    midpoint: Fraction,
    r_lo: Fraction,
    r_hi: Fraction,
    width: Fraction,
    subdivisions: int,
    m_cutoff: int,
    series_degree: int,
) -> tuple[list[arb], list[arb], list[arb], dict[str, arb]]:
    """Enclose I(r0), I'([rlo,rhi]), and uniform omitted-m errors."""

    point_model = DualScoreModel(
        arb_from_fraction(midpoint), arb(0), m_cutoff, series_degree
    )
    box_model = DualScoreModel(
        arb_box(r_lo, r_hi), arb(1), m_cutoff, series_degree
    )
    step = 2 * width / subdivisions
    step_a = arb_from_fraction(step)
    point = [arb(0) for _ in range(4)]
    derivative = [arb(0) for _ in range(4)]
    omitted = [arb(0) for _ in range(4)]
    kernel_tail = raw_base.kernel_m_tail_upper(m_cutoff).upper()
    max_kernel_prime_tail = arb(0)

    for index in range(subdivisions):
        left = -width + index * step
        w = arb_box(left, left + step)
        point_values = point_model.finite_integrands(w)
        box_values = box_model.finite_integrands(w)
        for order in range(4):
            point[order] += step_a * point_values[order].val
            derivative[order] += step_a * box_values[order].der

        density_upper = box_model.density(w).val.upper()
        n_dom_upper = abs_upper(box_model.dominant_score_excess(w).val)
        kernel_prime_tail = kernel_derivative_tail_upper(box_model, w).upper()
        if kernel_prime_tail > max_kernel_prime_tail:
            max_kernel_prime_tail = kernel_prime_tail
        w2_upper = (w * w).upper()
        mass_error = step_a * density_upper * kernel_tail
        score_error = step_a * density_upper * (
            kernel_tail * n_dom_upper + kernel_prime_tail
        )
        omitted[0] += mass_error
        omitted[1] += w2_upper * mass_error
        omitted[2] += score_error
        omitted[3] += (w2_upper + 2) * score_error

    diagnostics = {
        "kernel_tail_uniform": kernel_tail,
        "max_kernel_w_derivative_tail": max_kernel_prime_tail,
    }
    return point, derivative, omitted, diagnostics


def outer_tail_errors(model: DualScoreModel, width: Fraction) -> list[arb]:
    """Uniform H1272 score-tail errors on an r-box."""

    if width != 24:
        raise ValueError("the imported H1272 bounds require W=24")
    width_a = arb_from_fraction(width)
    tail_0 = arb("2e-40")
    tail_2 = arb("2e-37")
    tail_4 = arb("7e-35")
    tail_1 = tail_2 / width_a
    tail_3 = tail_4 / width_a
    boundary = (
        model.density(width_a).val.upper()
        + model.density(-width_a).val.upper()
    )
    return [
        tail_0,
        tail_2,
        boundary + tail_1,
        (width_a * width_a + 2) * boundary + tail_3,
    ]


def third_from_integrals(values: list[arb]) -> arb:
    mass, moment2_num, score_a_num, score_c_num = values
    second = moment2_num / mass
    score_a = score_a_num / mass
    score_c = score_c_num / mass
    return -score_c + 3 * score_a * second - 2 * score_a**3


def third_dual(values: list[Dual]) -> Dual:
    mass, moment2_num, score_a_num, score_c_num = values
    second = moment2_num / mass
    score_a = score_a_num / mass
    score_c = score_c_num / mass
    return -score_c + 3 * score_a * second - 2 * score_a**3


def tail_third_error(integral_ranges: list[arb], errors: list[arb]) -> arb:
    """Gradient bound for perturbing (Z,M2,SA,SC) by uniform errors."""

    expanded = [
        value + symmetric(error)
        for value, error in zip(integral_ranges, errors)
    ]
    z, m2, sa, sc = expanded
    if z.lower() <= 0:
        return arb("+inf")
    derivatives = [
        sc / z**2 - 6 * sa * m2 / z**3 + 6 * sa**3 / z**4,
        3 * sa / z**2,
        3 * m2 / z**2 - 6 * sa**2 / z**3,
        -1 / z,
    ]
    total = arb(0)
    for derivative, error in zip(derivatives, errors):
        total += abs_upper(derivative) * error.upper()
    return total


def box_certificate(
    r_lo: Fraction,
    r_hi: Fraction,
    width: Fraction,
    subdivisions: int,
    m_cutoff: int,
    series_degree: int,
    digits: int,
) -> dict[str, Any]:
    midpoint = (r_lo + r_hi) / 2
    radius = (r_hi - r_lo) / 2
    point, derivatives, omitted, diagnostics = central_point_and_derivative_sums(
        midpoint,
        r_lo,
        r_hi,
        width,
        subdivisions,
        m_cutoff,
        series_degree,
    )
    box_model = DualScoreModel(
        arb_box(r_lo, r_hi), arb(1), m_cutoff, series_degree
    )
    point_model = DualScoreModel(
        arb_from_fraction(midpoint), arb(0), m_cutoff, series_degree
    )
    outer = outer_tail_errors(box_model, width)
    errors = [left + right for left, right in zip(omitted, outer)]
    offset = symmetric(arb_from_fraction(radius))
    integral_ranges = [
        value + offset * derivative
        for value, derivative in zip(point, derivatives)
    ]

    point_third = third_from_integrals(point)
    integral_duals = [
        Dual(value, derivative)
        for value, derivative in zip(integral_ranges, derivatives)
    ]
    finite_third_dual = third_dual(integral_duals)
    scale_box = box_model.scale()
    scale_point = point_model.scale().val
    finite_f_point = scale_point * point_third
    finite_f_derivative = (
        scale_box.der * finite_third_dual.val
        + scale_box.val * finite_third_dual.der
    )
    perturbation_third = tail_third_error(integral_ranges, errors)
    perturbation_f = abs_upper(scale_box.val) * perturbation_third
    scaled = (
        finite_f_point
        + offset * finite_f_derivative
        + symmetric(perturbation_f)
    )
    margin = scaled + arb(3) / 4
    full_mass_lower = (
        integral_ranges[0] - errors[0]
    ).lower()

    checks = {
        "finite_integral_mass_positive": bool(integral_ranges[0].lower() > 0),
        "full_mass_positive_under_uniform_errors": bool(full_mass_lower > 0),
        "finite_scaled_interval": (
            str(scaled.lower()) != "nan" and str(scaled.upper()) != "nan"
        ),
        "target_minus_three_over_four": bool(margin.lower() > 0),
    }
    return {
        "r_lo": fraction_text(r_lo),
        "r_hi": fraction_text(r_hi),
        "r_mid": fraction_text(midpoint),
        "r_width": fraction_text(r_hi - r_lo),
        "subdivisions_W": subdivisions,
        "finite_integrals_at_midpoint": [
            interval_payload(value, digits) for value in point
        ],
        "finite_integral_r_derivatives": [
            interval_payload(value, digits) for value in derivatives
        ],
        "finite_integral_ranges_mean_value": [
            interval_payload(value, digits) for value in integral_ranges
        ],
        "uniform_integral_error_radii": [
            interval_payload(value, digits) for value in errors
        ],
        "finite_F_at_midpoint": interval_payload(finite_f_point, digits),
        "finite_F_r_derivative": interval_payload(finite_f_derivative, digits),
        "tail_third_error_radius": interval_payload(
            perturbation_third, digits
        ),
        "tail_F_error_radius": interval_payload(perturbation_f, digits),
        "scaled_t2_kappa_Xi_third": interval_payload(scaled, digits),
        "margin_over_minus_3_over_4": interval_payload(margin, digits),
        "omitted_m_diagnostics": {
            key: interval_payload(value, digits)
            for key, value in diagnostics.items()
        },
        "checks": checks,
        "passes": all(checks.values()),
    }


def seed_boxes(
    left: Fraction, right: Fraction, seed_width: Fraction
) -> list[tuple[Fraction, Fraction]]:
    boxes: list[tuple[Fraction, Fraction]] = []
    cursor = left
    while cursor < right:
        endpoint = min(cursor + seed_width, right)
        boxes.append((cursor, endpoint))
        cursor = endpoint
    return boxes


def adaptive_cover(args: argparse.Namespace) -> dict[str, Any]:
    dependency = json.loads(raw_base.H1272_JSON.read_text(encoding="utf-8"))
    theta_tail, theta_gate = score_base.global_kernel_unit_bound()
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
        (lo, hi, 0)
        for lo, hi in seed_boxes(r_min, r_max, seed_width)
    )
    accepted: list[dict[str, Any]] = []
    unresolved: list[dict[str, Any]] = []
    evaluations = 0
    while pending:
        r_lo, r_hi, depth = pending.popleft()
        result = box_certificate(
            r_lo,
            r_hi,
            width,
            args.subdivisions,
            args.m_cutoff,
            args.series_degree,
            args.digits,
        )
        result["depth"] = depth
        evaluations += 1
        if result["passes"]:
            accepted.append(result)
            continue
        if r_hi - r_lo > min_r_width and depth < args.max_r_depth:
            midpoint = (r_lo + r_hi) / 2
            pending.appendleft((midpoint, r_hi, depth + 1))
            pending.appendleft((r_lo, midpoint, depth + 1))
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
        "accepted_boxes_nonempty": bool(accepted),
        "no_unresolved_boxes": not unresolved and not pending,
        "exact_gap_free_r_chain": exact_chain,
        "every_accepted_box_passes": all(row["passes"] for row in accepted),
    }
    return {
        "id": "h1319_full_xi_score_mean_value_r_cover_arb",
        "classification": "h1319_full_xi_compact_mean_value_r_interval_certificate",
        "statement": "t(r)^2*kappa_Xi'''(t(r)) >= -3/4",
        "covered_r_interval": [fraction_text(r_min), fraction_text(r_max)],
        "configuration": {
            "W": fraction_text(width),
            "subdivisions_W": args.subdivisions,
            "m_cutoff": args.m_cutoff,
            "series_degree": args.series_degree,
            "seed_width": fraction_text(seed_width),
            "min_r_width": fraction_text(min_r_width),
            "max_r_depth": args.max_r_depth,
            "flint_dps": flint.ctx.dps,
            "finite_differences_used": False,
            "point_sampling_used": False,
            "r_method": "rigorous midpoint plus interval derivative mean-value form",
        },
        "summary": {
            "evaluations": evaluations,
            "accepted_box_count": len(accepted),
            "unresolved_box_count": len(unresolved),
            "worst_margin_box": worst,
        },
        "global_checks": global_checks,
        "accepted_boxes": accepted,
        "unresolved_boxes": unresolved,
        "all_checks_pass": all(global_checks.values()),
        "scope": (
            "compact full-Xi curvature certificate only; H920 assembly and "
            "the finite n<=126 certificate are separate dependencies"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r-min", default="4")
    parser.add_argument("--r-max", default="71")
    parser.add_argument("--seed-width", default="1")
    parser.add_argument("--min-r-width", default="1/1024")
    parser.add_argument("--max-r-depth", type=int, default=12)
    parser.add_argument("--W", default="24")
    parser.add_argument("--subdivisions", type=int, default=4000)
    parser.add_argument("--m-cutoff", type=int, default=8)
    parser.add_argument("--series-degree", type=int, default=48)
    parser.add_argument("--dps", type=int, default=90)
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
