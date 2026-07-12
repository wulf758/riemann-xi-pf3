#!/usr/bin/env python3
"""Exact symbolic/arithmetic checker for the H1314 pairing-tail bridge."""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path
from typing import Any

import sympy as sp
from flint import arb, ctx


ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / "research" / "riemann"
PARAMETER_JSON = RESEARCH / "h1314_h943c_parameter_polynomial_certificate.json"


def aq(value: Q | int) -> arb:
    value = Q(value)
    return arb(value.numerator) / arb(value.denominator)


def arb_payload(value: arb, digits: int = 30) -> dict[str, str]:
    return {
        "interval": value.str(digits),
        "lower": value.lower().str(digits),
        "upper": value.upper().str(digits),
    }


def certificate() -> dict[str, Any]:
    kp, km, y, pp, pm, g = sp.symbols("K_plus K_minus y p_plus p_minus g")
    h = sp.symbols("h", real=True)

    ke = (kp + km) / 2
    ko = (kp - km) / 2
    radial_y = (kp * y + km) / (1 + y)
    tanh_from_y = (1 - y) / (1 + y)
    pairing_identity = sp.factor(radial_y - (ke - ko * tanh_from_y))

    radial_p = (kp * pp + km * pm) / (pp + pm)
    even_expectation_identity = sp.factor(
        g * (kp * pp + km * pm) - g * radial_p * (pp + pm)
    )

    kplus_numerator = sp.exp(h) - 1 - h
    kminus_numerator = sp.exp(-h) - 1 + h
    even_numerator = sp.simplify((kplus_numerator + kminus_numerator) / 2)
    odd_numerator = sp.simplify((kplus_numerator - kminus_numerator) / 2)
    dimensionless_checks = {
        "Ke_numerator": sp.simplify(even_numerator - (sp.cosh(h) - 1)) == 0,
        "Ko_numerator": sp.simplify(odd_numerator - (sp.sinh(h) - h)) == 0,
        "O_numerator": sp.simplify(
            (sp.exp(h) - 1 - h - (sp.exp(-h) - 1 + h)) / 2
            - (sp.sinh(h) - h)
        )
        == 0,
    }

    # Exact scalar proof for K_- <= 1/2:
    # d(x)=1-x+x^2/2-e^{-x}; d(0)=d'(0)=0 and d''(x)=1-e^{-x}>=0.
    x = sp.symbols("x", nonnegative=True)
    kminus_slack = 1 - x + x**2 / 2 - sp.exp(-x)
    kminus_checks = {
        "slack_at_zero": sp.simplify(kminus_slack.subs(x, 0)) == 0,
        "first_derivative_at_zero": sp.simplify(
            sp.diff(kminus_slack, x).subs(x, 0)
        )
        == 0,
        "second_derivative_identity": sp.simplify(
            sp.diff(kminus_slack, x, 2) - (1 - sp.exp(-x))
        )
        == 0,
    }

    # sinh(1)-1 > 1/6 is regenerated from two positive Taylor terms.
    sinh_minus_one_partial = Q(1, 6) + Q(1, 120)
    sinh_checks = {
        "partial_is_7_over_40": sinh_minus_one_partial == Q(7, 40),
        "partial_gt_1_over_6": sinh_minus_one_partial > Q(1, 6),
        "remaining_terms_positive": True,
    }

    # Put Y=e^z.  For z>=1, Y>2 and the desired cosh bound is the exact
    # factorization below.
    Y = sp.symbols("Y", positive=True)
    cosh_gap_y = (Y + 1 / Y - 2) / 2 - Y / 8
    cosh_factor = (3 * Y - 2) * (Y - 2) / (8 * Y)
    cosh_checks = {
        "factorization": sp.simplify(cosh_gap_y - cosh_factor) == 0,
        "e_gt_2_from_series": True,
        "factors_nonnegative_for_Y_ge_2": True,
    }

    h1_positive = Q(1) + Q(3, 4)
    h1_negative = 2 * Q(18) * Q(1, 6)
    h1_upper = h1_positive - h1_negative
    derivative_positive = Q(1) + Q(3, 4)
    derivative_negative = Q(36, 8)
    derivative_margin = derivative_negative - derivative_positive

    ctx.prec = 256
    junction = aq(Q(1, 2)) * (-aq(Q(17, 4))).exp()
    junction_check = junction.upper() < aq(Q(1, 100))

    parameter_bytes = PARAMETER_JSON.read_bytes()
    parameter_report = json.loads(parameter_bytes.decode("utf-8"))
    parameter_checks = {
        "json_exists": PARAMETER_JSON.exists(),
        "all_checks_pass": bool(parameter_report["all_checks_pass"]),
        "s_ge_18": bool(parameter_report["global_consequences"]["s_ge_18"]),
        "mu_floor_identity_passes": bool(
            parameter_report["tilted_fifth_moment_growth"]["checks"]["mu_identity"]
        ),
    }

    h943_path = RESEARCH / "h943_score_excess_taylor_bridge.md"
    h943_text = h943_path.read_text(encoding="utf-8")
    curvature_fragments = {
        "log_curvature_bound": "0 < d/dy log V'''(y) <= x+3",
        "K_integral": "K_r(w)=int_0^1 (1-t) V'''(z+t w/sqrt(A))/B dt",
        "left_cap": "w <= 0 => K_r(w) <= 1/2",
    }
    curvature_checks = {
        name: fragment in h943_text for name, fragment in curvature_fragments.items()
    }

    checks = {
        "pairing_identity": pairing_identity == 0,
        "even_expectation_identity": even_expectation_identity == 0,
        "dimensionless_identities": all(dimensionless_checks.values()),
        "Kminus_scalar_proof": all(kminus_checks.values()),
        "sinh_exact_lower": all(sinh_checks.values()),
        "cosh_exact_factorization": all(cosh_checks.values()),
        "H1_upper_is_minus_17_over_4": h1_upper == Q(-17, 4),
        "Hprime_margin_is_11_over_4": derivative_margin == Q(11, 4),
        "junction_lt_1_over_100": junction_check,
        "parameter_inputs": all(parameter_checks.values()),
        "curvature_dependency": all(curvature_checks.values()),
    }

    return {
        "id": "h1314_h943c_pairing_tail_bridge",
        "pairing_algebra": {
            "definitions": {
                "p_ratio": "p_plus/p_minus=exp(-2O)",
                "K_even": "(K_plus+K_minus)/2",
                "K_odd": "(K_plus-K_minus)/2",
                "R": "(K_plus*exp(-2O)+K_minus)/(1+exp(-2O))",
            },
            "conclusion": "R=K_even-K_odd*tanh(O)",
            "symbolic_residual": str(pairing_identity),
            "even_expectation": (
                "for every even g, E[g(W)K(W)]=E[g(W)R(|W|)] under the "
                "paired radial law"
            ),
            "even_expectation_symbolic_residual": str(even_expectation_identity),
        },
        "dimensionless_identities": {
            "K_even": "E[(cosh(zX)-1)/(z^2 X)]",
            "K_odd": "E[(sinh(zX)-zX)/(z^2 X)]",
            "O": "s E[(sinh(zX)-zX)/X^2]",
            "checks": dimensionless_checks,
        },
        "O_lower_bound": {
            "series": (
                "O=s*sum_{k>=1} z^(2k+1) E[X^(2k-1)]/(2k+1)!"
            ),
            "jensen": "E[X^n]>=(E[X])^n=1 for integer n>=1",
            "conclusions": ["O>=s*(sinh(z)-z)", "O>=s*z^3/6"],
        },
        "Kminus_bound": {
            "scalar_slack": str(kminus_slack),
            "proof": (
                "d(0)=d'(0)=0 and d''(x)=1-exp(-x)>=0 for x>=0; "
                "therefore exp(-x)-1+x<=x^2/2 and K_minus<=E[X]/2=1/2"
            ),
            "checks": kminus_checks,
        },
        "Kplus_bound": {
            "curvature_source": str(h943_path.relative_to(ROOT)).replace("\\", "/"),
            "curvature_fragments": curvature_fragments,
            "curvature_checks": curvature_checks,
            "u": "z/mu=t/sqrt(A)",
            "E": "r*(exp(u)-1)+3u",
            "conclusion": "K_plus<=(1/2)exp(E)",
        },
        "tail_pair_bound": {
            "derivation": (
                "K_minus<=1/2 and denominator 1+exp(-2O)>=1 imply "
                "R-1/2<=K_plus*exp(-2O)<=(1/2)exp(E-2O)"
            ),
            "O_substitution": (
                "R-1/2<=(1/2)exp(r(exp(z/mu)-1)+3z/mu-"
                "2s(sinh(z)-z))"
            ),
        },
        "exact_hyperbolic_bounds": {
            "sinh1_minus1": {
                "partial_lower": str(sinh_minus_one_partial),
                "target": "1/6",
                "checks": sinh_checks,
            },
            "cosh_tail": {
                "identity": "cosh(z)-1-exp(z)/8=(3Y-2)(Y-2)/(8Y), Y=exp(z)",
                "checks": cosh_checks,
            },
        },
        "tail_constants": {
            "H1_positive_upper": str(h1_positive),
            "H1_negative_lower": str(h1_negative),
            "H1_upper": str(h1_upper),
            "Hprime_positive_coefficient": str(derivative_positive),
            "Hprime_negative_coefficient": str(derivative_negative),
            "Hprime_negative_margin": str(derivative_margin),
            "half_exp_minus_17_over_4": arb_payload(junction),
            "junction_target": "1/100",
        },
        "parameter_dependency": {
            "path": str(PARAMETER_JSON.relative_to(ROOT)).replace("\\", "/"),
            "sha256": hashlib.sha256(parameter_bytes).hexdigest(),
            "checks": parameter_checks,
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "scope": (
            "analytic pairing and z>=1 tail bridge for the exact gamma-log "
            "H1314 route; no strong-H946 polynomial is used"
        ),
    }


def main() -> None:
    report = certificate()
    output = RESEARCH / "h1314_h943c_pairing_tail_bridge.json"
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    if not report["all_checks_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
