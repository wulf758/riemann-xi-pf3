#!/usr/bin/env python3
"""Stable high-r variant of the H1319 fixed-r full-Xi Arb pilot."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import flint
from flint import arb

import rh_h1319_full_xi_fixed_r_arb_pilot as base
from rh_h930_fixed_r_arb_tail_certificate import arb_from_fraction, interval_payload


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "research" / "riemann" / "h1319_full_xi_fixed_r_stable_arb.json"


class StableFullXiModel(base.FullXiModel):
    def psi(self, w: arb) -> arb:
        h = w / self.sqrt_A
        expm1_h = h.expm1()
        delta = self.r * expm1_h
        return (
            self.pi * self.exp_r * (delta.expm1() - delta)
            + (self.pi * self.r * self.exp_r - self.alpha * self.r)
            * (expm1_h - h)
        )


def certificate(
    r_value: Fraction,
    width: Fraction,
    subdivisions: int,
    m_cutoff: int,
    digits: int,
) -> dict[str, Any]:
    dependency = json.loads(base.H1272_JSON.read_text(encoding="utf-8"))
    model = StableFullXiModel(arb_from_fraction(r_value), m_cutoff)
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
        "id": "h1319_full_xi_fixed_r_stable_arb",
        "classification": "h1319_full_xi_fixed_r_stable_point_certificate_not_r_cover",
        "config": {
            "r": str(r_value),
            "W": str(width),
            "subdivisions": subdivisions,
            "m_cutoff": m_cutoff,
            "flint_dps": flint.ctx.dps,
            "psi_form": (
                "pi*exp(r)*(expm1(delta)-delta)+"
                "(pi*r*exp(r)-alpha*r)*(expm1(h)-h)"
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
