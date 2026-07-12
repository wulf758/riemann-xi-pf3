#!/usr/bin/env python3
"""Rigorous fixed-r H1319 pilot using a cancellation-aware Xi score.

For the full Xi saddle density

    g(w) = exp(-Psi(w)) H(w),

the score excess is N_Xi = Psi' - H'/H - w.  Stein integration by
parts gives

    E[W]   = -E[N_Xi],
    E[W^3] = -E[(W^2+2)N_Xi].

Consequently the third centred moment can be evaluated from integrals of
exp(-Psi) * (H*N_dom-H') without subtracting nearly equal raw moments.  The
dominant score excess N_dom is also evaluated without cancellation.

This script is a rigorous certificate at one rational r.  It is deliberately
not advertised as an r-interval cover.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import flint
from flint import arb

import rh_h1319_full_xi_fixed_r_arb_pilot as base
from rh_h1319_full_xi_fixed_r_series_arb import (
    SeriesFullXiModel,
    expm1_minus_x_series,
)
from rh_h930_fixed_r_arb_tail_certificate import (
    arb_box,
    arb_from_fraction,
    interval_payload,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = (
    ROOT / "research" / "riemann" / "h1319_full_xi_fixed_r_score_arb.json"
)


class ScoreFullXiModel(SeriesFullXiModel):
    """Full-Xi model with stable dominant score and kernel derivative."""

    def dominant_score_excess(self, w: arb) -> arb:
        """Return Psi'(w)-w with every small subtraction removed.

        Put h=w/sqrt(A), delta=r*(exp(h)-1), and phi2(x)=exp(x)-1-x.
        Direct algebra gives the exact identity

          sqrt(A) N_dom
            = A phi2(h)
              + pi exp(r) [r phi2(delta)+delta expm1(delta)].
        """

        h = w / self.sqrt_A
        delta = self.r * h.expm1()
        phi_h = expm1_minus_x_series(h, self.series_degree)
        phi_delta = expm1_minus_x_series(delta, self.series_degree)
        return (
            self.A * phi_h
            + self.pi
            * self.exp_r
            * (self.r * phi_delta + delta * delta.expm1())
        ) / self.sqrt_A

    def kernel_derivative_tail_upper(self, w: arb) -> arb:
        """Bound the absolute derivative tail after ``m_cutoff``.

        On a finite w-box, X=pi*exp(rho)>=pi.  For d=m^2-1 the absolute
        X-derivative of one omitted summand is bounded by

          exp(-d X) m^6 (1+3/(2X)+3/(2X^2)).

        The displayed m-series is bounded geometrically from its first term.
        """

        rho = self.r * (w / self.sqrt_A).exp()
        x = self.pi * rho.exp()
        x_lo = x.lower()
        scale_hi = (x * rho / self.sqrt_A).upper()
        first_m = self.m_cutoff + 1
        m = arb(first_m)
        coefficient = 1 + arb(3) / (2 * x_lo) + arb(3) / (2 * x_lo * x_lo)
        first = coefficient * m**6 * (-x_lo * (first_m * first_m - 1)).exp()
        ratio = (
            arb(first_m + 1) / arb(first_m)
        ) ** 6 * (-x_lo * (2 * first_m + 1)).exp()
        return scale_hi * first / (1 - ratio)

    def kernel_ratio_and_derivative(self, w: arb) -> tuple[arb, arb, arb]:
        """Enclose H(w), H'(w), and the H' omitted-m tail radius."""

        rho = self.r * (w / self.sqrt_A).exp()
        x = self.pi * rho.exp()
        x_prime = x * rho / self.sqrt_A
        value = arb(0)
        derivative = arb(0)
        for m_int in range(1, self.m_cutoff + 1):
            m = arb(m_int)
            d = arb(m_int * m_int - 1)
            exponential = (-x * d).exp()
            coefficient = m**4 - arb(3) * m**2 / (2 * x)
            value += coefficient * exponential
            derivative += x_prime * exponential * (
                arb(3) * m**2 / (2 * x * x) - d * coefficient
            )

        value += base.nonnegative_error(self.kernel_tail)
        derivative_tail = self.kernel_derivative_tail_upper(w)
        derivative += base.symmetric_error(derivative_tail)
        return value, derivative, derivative_tail


def global_kernel_unit_bound() -> tuple[arb, arb]:
    """Return a theta-tail bound proving 0<H(X)<1 for X>=pi."""

    pi = arb.pi()
    finite = sum(
        arb(m) ** 4 * (-pi * (m * m - 1)).exp()
        for m in range(2, 20)
    )
    omitted = arb(20) ** 4 * (-pi * (20**2 - 1)).exp() / (
        1 - (-4 * pi).exp()
    )
    return finite + omitted, arb(3) / (2 * pi)


def central_score_sums(
    model: ScoreFullXiModel, width: Fraction, subdivisions: int
) -> tuple[list[arb], dict[str, arb]]:
    """Integrate Z, M2, numerator(A), numerator(C) on [-W,W]."""

    step = (2 * width) / subdivisions
    step_a = arb_from_fraction(step)
    sums = [arb(0) for _ in range(4)]
    min_h_lower: arb | None = None
    max_h_upper = arb(0)
    max_derivative_tail = arb(0)

    for index in range(subdivisions):
        left = -width + index * step
        w = arb_box(left, left + step)
        density = model.density(w)
        kernel, kernel_prime, derivative_tail = model.kernel_ratio_and_derivative(w)
        n_dom = model.dominant_score_excess(w)
        score_numerator = density * (kernel * n_dom - kernel_prime)
        w2 = w * w

        sums[0] += step_a * density * kernel
        sums[1] += step_a * w2 * density * kernel
        sums[2] += step_a * score_numerator
        sums[3] += step_a * (w2 + 2) * score_numerator

        if min_h_lower is None or kernel.lower() < min_h_lower:
            min_h_lower = kernel.lower()
        if kernel.upper() > max_h_upper:
            max_h_upper = kernel.upper()
        if derivative_tail.upper() > max_derivative_tail:
            max_derivative_tail = derivative_tail.upper()

    assert min_h_lower is not None
    diagnostics = {
        "kernel_min_lower": min_h_lower,
        "kernel_max_upper": max_h_upper,
        "max_kernel_derivative_tail_upper": max_derivative_tail,
    }
    return sums, diagnostics


def attach_score_tails(
    model: ScoreFullXiModel, central: list[arb], width: Fraction
) -> tuple[list[arb], dict[str, arb]]:
    """Attach H1272 tails, using score integration by parts.

    If g=exp(-Psi)H, then g*N_Xi=-g'-w*g.  Over both tails this gives

      |int g N_Xi| <= g(W)+g(-W)+int_tail |w|g,
      |int (w^2+2)g N_Xi|
        <= (W^2+2)(g(W)+g(-W))+int_tail |w|^3 g.

    Since 0<H<1, g<=exp(-Psi), and H1272 supplies the even tails.
    """

    if width != 24:
        raise ValueError("the imported H1272 bounds require W=24")
    width_a = arb_from_fraction(width)
    tail_0 = arb("2e-40")
    tail_2 = arb("2e-37")
    tail_4 = arb("7e-35")
    tail_1 = tail_2 / width_a
    tail_3 = tail_4 / width_a
    boundary = model.density(width_a).upper() + model.density(-width_a).upper()
    tail_a = boundary + tail_1
    tail_c = (width_a * width_a + 2) * boundary + tail_3

    totals = [
        central[0] + base.nonnegative_error(tail_0),
        central[1] + base.nonnegative_error(tail_2),
        central[2] + base.symmetric_error(tail_a),
        central[3] + base.symmetric_error(tail_c),
    ]
    return totals, {
        "dominant_boundary_density_sum_upper": boundary,
        "mass_tail_upper": tail_0,
        "second_moment_tail_upper": tail_2,
        "score_A_tail_radius": tail_a,
        "score_C_tail_radius": tail_c,
    }


def certificate(
    r_value: Fraction,
    width: Fraction,
    subdivisions: int,
    m_cutoff: int,
    series_degree: int,
    digits: int,
) -> dict[str, Any]:
    dependency = json.loads(base.H1272_JSON.read_text(encoding="utf-8"))
    model = ScoreFullXiModel(
        arb_from_fraction(r_value), m_cutoff, series_degree
    )
    central, diagnostics = central_score_sums(model, width, subdivisions)
    totals, tail_data = attach_score_tails(model, central, width)

    mass, moment2_num, score_a_num, score_c_num = totals
    second = moment2_num / mass
    score_a = score_a_num / mass
    score_c = score_c_num / mass
    third = -score_c + 3 * score_a * second - 2 * score_a**3
    t = model.q / 2
    kappa_3_t = 8 * third / (model.A * model.sqrt_A)
    scaled = t**2 * kappa_3_t
    margin = scaled + arb(3) / 4

    theta_tail, theta_gate = global_kernel_unit_bound()
    checks = {
        "H1272_dependency": bool(dependency.get("all_checks_pass")),
        "mass_positive": bool(mass.lower() > 0),
        "kernel_positive_on_central_boxes": bool(
            diagnostics["kernel_min_lower"] > 0
        ),
        "kernel_below_one_on_central_boxes": bool(
            diagnostics["kernel_max_upper"] < 1
        ),
        "global_kernel_unit_bound": bool(theta_tail.upper() < theta_gate.lower()),
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
        "id": "h1319_full_xi_fixed_r_score_arb",
        "classification": "h1319_full_xi_fixed_r_score_point_certificate_not_r_cover",
        "config": {
            "r": str(r_value),
            "W": str(width),
            "subdivisions": subdivisions,
            "m_cutoff": m_cutoff,
            "series_degree": series_degree,
            "flint_dps": flint.ctx.dps,
            "finite_differences_used": False,
        },
        "identity": {
            "score": "N_Xi=Psi'-H'/H-w",
            "A": "E[N_Xi]=-E[W]",
            "C": "E[(W^2+2)N_Xi]=-E[W^3]",
            "central_third": "-C+3*A*E[W^2]-2*A^3",
        },
        "saddle": {
            "q": interval_payload(model.q, digits),
            "t": interval_payload(t, digits),
            "A": interval_payload(model.A, digits),
            "a": interval_payload(model.a, digits),
        },
        "central_integrals": {
            "mass": interval_payload(central[0], digits),
            "second_moment_numerator": interval_payload(central[1], digits),
            "score_A_numerator": interval_payload(central[2], digits),
            "score_C_numerator": interval_payload(central[3], digits),
        },
        "full_integrals_with_tails": {
            "mass": interval_payload(mass, digits),
            "second_moment_numerator": interval_payload(moment2_num, digits),
            "score_A_numerator": interval_payload(score_a_num, digits),
            "score_C_numerator": interval_payload(score_c_num, digits),
        },
        "normalized": {
            "second_moment": interval_payload(second, digits),
            "score_A": interval_payload(score_a, digits),
            "score_C": interval_payload(score_c, digits),
            "third_central_moment_W": interval_payload(third, digits),
        },
        "scaled_t2_kappa_Xi_third": interval_payload(scaled, digits),
        "margin_over_minus_3_over_4": interval_payload(margin, digits),
        "kernel_diagnostics": {
            key: interval_payload(value, digits)
            for key, value in diagnostics.items()
        },
        "tail_data": {
            key: interval_payload(value, digits) for key, value in tail_data.items()
        },
        "global_kernel_bound": {
            "theta_tail_at_pi": interval_payload(theta_tail, digits),
            "three_over_two_pi": interval_payload(theta_gate, digits),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "scope": "rigorous at one rational r only; no interval r-cover",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r", default="647/200")
    parser.add_argument("--W", default="24")
    parser.add_argument("--subdivisions", type=int, default=12000)
    parser.add_argument("--m-cutoff", type=int, default=8)
    parser.add_argument("--series-degree", type=int, default=48)
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--digits", type=int, default=24)
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    flint.ctx.dps = args.dps
    report = certificate(
        Fraction(args.r),
        Fraction(args.W),
        args.subdivisions,
        args.m_cutoff,
        args.series_degree,
        args.digits,
    )
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
