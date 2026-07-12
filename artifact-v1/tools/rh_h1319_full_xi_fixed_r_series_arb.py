#!/usr/bin/env python3
"""Cancellation-free fixed-r H1319 pilot using a certified exp remainder series."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from math import factorial
from pathlib import Path
from typing import Any

import flint
from flint import arb

import rh_h1319_full_xi_fixed_r_arb_pilot as base
from rh_h930_fixed_r_arb_tail_certificate import arb_from_fraction, interval_payload


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "research" / "riemann" / "h1319_full_xi_fixed_r_series_arb.json"


def expm1_minus_x_series(x: arb, degree: int) -> arb:
    """Enclose exp(x)-1-x without interval cancellation."""

    if degree < 2:
        raise ValueError("degree must be at least two")
    term = x * x / 2
    total = term
    for order in range(3, degree + 1):
        term = term * x / order
        total += term
    magnitude = abs(x).upper()
    remainder = (
        magnitude.exp()
        * magnitude ** (degree + 1)
        / arb(factorial(degree + 1))
    )
    return total + arb(0, 1) * remainder


class SeriesFullXiModel(base.FullXiModel):
    def __init__(self, r: arb, m_cutoff: int, series_degree: int):
        super().__init__(r, m_cutoff)
        self.series_degree = series_degree

    def psi(self, w: arb) -> arb:
        h = w / self.sqrt_A
        delta = self.r * h.expm1()
        return (
            self.pi
            * self.exp_r
            * expm1_minus_x_series(delta, self.series_degree)
            + (self.pi * self.r * self.exp_r - self.alpha * self.r)
            * expm1_minus_x_series(h, self.series_degree)
        )


def certificate(
    r_value: Fraction,
    width: Fraction,
    subdivisions: int,
    m_cutoff: int,
    series_degree: int,
    digits: int,
) -> dict[str, Any]:
    dependency = json.loads(base.H1272_JSON.read_text(encoding="utf-8"))
    model = SeriesFullXiModel(
        arb_from_fraction(r_value), m_cutoff, series_degree
    )
    central = base.central_moment_sums(model, width, subdivisions)
    moments = base.attach_global_tails(central, width)

    mass = moments[0]
    raw_1 = moments[1] / mass
    raw_2 = moments[2] / mass
    raw_3 = moments[3] / mass
    third = raw_3 - 3 * raw_1 * raw_2 + 2 * raw_1**3
    t = model.q / 2
    kappa_3_t = 8 * third / (model.A * model.sqrt_A)
    scaled = t**2 * kappa_3_t
    margin = scaled + arb(3) / 4

    checks = {
        "H1272_dependency": bool(dependency.get("all_checks_pass")),
        "mass_positive": bool(mass.lower() > 0),
        "kernel_tail_positive": bool(model.kernel_tail.lower() > 0),
        "kernel_tail_small": bool(model.kernel_tail.upper() < arb("1e-100")),
        "finite_scaled_interval": str(scaled.lower()) != "nan" and str(scaled.upper()) != "nan",
        "target_minus_three_over_four": bool(margin.lower() > 0),
    }
    return {
        "id": "h1319_full_xi_fixed_r_series_arb",
        "classification": "h1319_full_xi_fixed_r_series_point_certificate_not_r_cover",
        "config": {
            "r": str(r_value),
            "W": str(width),
            "subdivisions": subdivisions,
            "m_cutoff": m_cutoff,
            "series_degree": series_degree,
            "flint_dps": flint.ctx.dps,
            "psi_form": (
                "pi*exp(r)*phi2(delta)+(pi*r*exp(r)-alpha*r)*phi2(h), "
                "phi2(x)=sum_(k>=2)x^k/k!"
            ),
        },
        "saddle": {
            "q": interval_payload(model.q, digits),
            "t": interval_payload(t, digits),
            "A": interval_payload(model.A, digits),
        },
        "kernel_m_tail_upper": interval_payload(model.kernel_tail, digits),
        "central_integrals": {
            str(order): interval_payload(value, digits)
            for order, value in enumerate(central)
        },
        "full_integrals_with_tails": {
            str(order): interval_payload(value, digits)
            for order, value in enumerate(moments)
        },
        "third_central_moment_W": interval_payload(third, digits),
        "kappa_Xi_third_t": interval_payload(kappa_3_t, digits),
        "scaled_t2_kappa_Xi_third": interval_payload(scaled, digits),
        "margin_over_minus_3_over_4": interval_payload(margin, digits),
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "scope": "rigorous at one rational r only; no interval r-cover",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r", default="647/200")
    parser.add_argument("--W", default="24")
    parser.add_argument("--subdivisions", type=int, default=50000)
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
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
