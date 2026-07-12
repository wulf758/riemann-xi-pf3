#!/usr/bin/env python3
"""Rigorous fixed-r Arb pilot for the full-Xi third cumulant.

This is a point certificate, not yet a cover in r.  It integrates the exact
full Xi kernel as a positive multiplicative tilt of the dominant gamma-log
law, uses first-order Arb boxes on |W|<=24, and imports the proved H1272
dominant tail bounds.  No numerical differentiation is used.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import flint
from flint import arb

from rh_h930_fixed_r_arb_tail_certificate import (
    GammaLogFixedR,
    arb_box,
    arb_from_fraction,
    interval_payload,
)


ROOT = Path(__file__).resolve().parents[1]
H1272_JSON = ROOT / "research" / "riemann" / "h1272_gamma_log_even_moment_tail_global.json"
DEFAULT_OUT = ROOT / "research" / "riemann" / "h1319_full_xi_fixed_r_arb_pilot.json"


def symmetric_error(radius: arb) -> arb:
    return arb(0, 1) * radius


def nonnegative_error(upper: arb) -> arb:
    return upper / 2 + arb(0, 1) * (upper / 2)


def kernel_m_tail_upper(m_cutoff: int) -> arb:
    """Uniform upper bound for the omitted positive m-sum, using X>=pi."""

    first_m = m_cutoff + 1
    pi = arb.pi()
    first = arb(first_m) ** 4 * (-pi * (first_m * first_m - 1)).exp()
    ratio = arb(16) * (-pi * (2 * first_m + 1)).exp()
    return first / (1 - ratio)


class FullXiModel(GammaLogFixedR):
    def __init__(self, r: arb, m_cutoff: int):
        super().__init__(r)
        self.m_cutoff = m_cutoff
        self.kernel_tail = kernel_m_tail_upper(m_cutoff)

    def kernel_ratio(self, w: arb) -> arb:
        rho = self.r * (w / self.sqrt_A).exp()
        x = self.pi * rho.exp()
        total = arb(0)
        for m in range(1, self.m_cutoff + 1):
            m_a = arb(m)
            coefficient = m_a**4 - arb(3) * m_a**2 / (2 * x)
            total += coefficient * (-x * (m * m - 1)).exp()
        return total + nonnegative_error(self.kernel_tail)


def central_moment_sums(
    model: FullXiModel, width: Fraction, subdivisions: int
) -> list[arb]:
    step = (2 * width) / subdivisions
    step_a = arb_from_fraction(step)
    sums = [arb(0) for _ in range(4)]
    for index in range(subdivisions):
        left = -width + index * step
        w = arb_box(left, left + step)
        weight = model.density(w) * model.kernel_ratio(w)
        power = arb(1)
        for order in range(4):
            sums[order] += step_a * power * weight
            power *= w
    return sums


def attach_global_tails(central: list[arb], width: Fraction) -> list[arb]:
    # H1272 proves the coarser strict targets below for the dominant density.
    # The pointwise Xi ratio satisfies 0<H<1.  On |W|>=24,
    # |W|<=W^2/24 and |W|^3<=W^4/24.
    if width != 24:
        raise ValueError("the imported H1272 bounds require W=24")
    tail_0 = arb("2e-40")
    tail_2 = arb("2e-37")
    tail_4 = arb("7e-35")
    tail_1 = tail_2 / 24
    tail_3 = tail_4 / 24
    return [
        central[0] + nonnegative_error(tail_0),
        central[1] + symmetric_error(tail_1),
        central[2] + nonnegative_error(tail_2),
        central[3] + symmetric_error(tail_3),
    ]


def certificate(
    r_value: Fraction, width: Fraction, subdivisions: int, m_cutoff: int, digits: int
) -> dict[str, Any]:
    dependency = json.loads(H1272_JSON.read_text(encoding="utf-8"))
    model = FullXiModel(arb_from_fraction(r_value), m_cutoff)
    central = central_moment_sums(model, width, subdivisions)
    moments = attach_global_tails(central, width)

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
        "target_minus_three_over_four": bool(margin.lower() > 0),
    }
    return {
        "id": "h1319_full_xi_fixed_r_arb_pilot",
        "classification": "h1319_full_xi_fixed_r_point_certificate_not_r_cover",
        "config": {
            "r": str(r_value),
            "W": str(width),
            "subdivisions": subdivisions,
            "m_cutoff": m_cutoff,
            "flint_dps": flint.ctx.dps,
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
        "scope": (
            "rigorous at one rational r only; an interval r-cover and the "
            "piecewise H920 assembly remain separate"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r", default="647/200")
    parser.add_argument("--W", default="24")
    parser.add_argument("--subdivisions", type=int, default=4000)
    parser.add_argument("--m-cutoff", type=int, default=8)
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
        args.digits,
    )
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
