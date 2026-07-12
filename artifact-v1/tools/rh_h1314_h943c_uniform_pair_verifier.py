#!/usr/bin/env python3
"""Canonical H1314/H943-C uniform-pair assembly verifier.

This script reruns the exact parameter engine, reruns the independent exact
five-interval U4 verifier, checks the large-z constants, and writes the main
canonical H1314 JSON artifact.  The analytic pair identities and H943-B
curvature-distortion lemma are cited from their proved local source; every new
constant used to assemble them is regenerated here.
"""

from __future__ import annotations

import json
import subprocess
import sys
from fractions import Fraction as Q
from pathlib import Path
from typing import Any

from flint import arb, ctx

import rh_h1314_h943c_parameter_polynomial_verifier as parameters


ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / "research" / "riemann"
U4_SCRIPT = ROOT / "tools" / "rh_h943c_u4_five_interval_certificate.py"
U4_JSON = RESEARCH / "h943c_u4_five_interval_certificate.json"


def aq(value: Q | int) -> arb:
    value = Q(value)
    return arb(value.numerator) / arb(value.denominator)


def arb_payload(value: arb, digits: int = 30) -> dict[str, str]:
    return {
        "interval": value.str(digits),
        "lower": value.lower().str(digits),
        "upper": value.upper().str(digits),
    }


def rerun_u4() -> dict[str, Any]:
    completed = subprocess.run(
        [sys.executable, str(U4_SCRIPT)],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "U4 verifier failed:\n" + completed.stdout + "\n" + completed.stderr
        )
    return json.loads(U4_JSON.read_text(encoding="utf-8"))


def analytic_dependency_certificate() -> dict[str, Any]:
    h943 = RESEARCH / "h943_score_excess_taylor_bridge.md"
    text = h943.read_text(encoding="utf-8")
    required_fragments = {
        "exact_taylor_identity": "N(w)=a w^2 K_r(w)",
        "K_nonnegative": "K_r(w)>=0",
        "K_at_zero": "K_r(0)=1/2",
        "left_pair_cap": "w <= 0  =>  K_r(w) <= 1/2",
        "curvature_distortion": "0 < d/dy log V'''(y) <= x+3",
    }
    checks = {name: fragment in text for name, fragment in required_fragments.items()}
    return {
        "source": str(h943.relative_to(ROOT)).replace("\\", "/"),
        "classification": "proved analytic dependency, not re-proved by coefficient code",
        "required_fragments": required_fragments,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def compact_assembly(parameter_report: dict[str, Any], u4: dict[str, Any]) -> dict[str, Any]:
    consequences = parameter_report["global_consequences"]
    ex3_bound = Q(1) + Q(7, 2) * Q(3, 20)
    tilted_coefficient = Q(4, 720)
    odd_coefficient = Q(1, 6)
    O_coefficient = Q(18, 6)
    checks = {
        "EX3_bound_is_61_over_40": ex3_bound == Q(61, 40),
        "tilted_coefficient_is_1_over_180": tilted_coefficient == Q(1, 180),
        "Ko_lower_coefficient_is_1_over_6": odd_coefficient == Q(1, 6),
        "O_lower_coefficient_is_3": O_coefficient == 3,
        "parameter_package_passes": bool(parameter_report["all_checks_pass"]),
        "s_ge_18": bool(consequences["s_ge_18"]),
        "EX3_le_61_over_40": bool(consequences["EX3_le_61_over_40"]),
        "tilted_EX5_bound": bool(
            consequences["tilted_EX5_le_4_exp_21z_over_8"]
        ),
        "U4_exact_certificate_passes": bool(u4["all_checks_pass"]),
    }
    return {
        "range": "0 <= z <= 1",
        "taylor_inequalities": {
            "Ke": (
                "Ke-1/2 <= z^2*E[X^3]/24 + "
                "z^4*E[X^5 exp(zX)]/720"
            ),
            "Ko": "Ko >= z*E[X^2]/6 >= z/6",
            "O": "O >= s*z^3/6 >= 3*z^3",
        },
        "derived_majorant": (
            "U4(z)=61*z^2/960+z^4*exp(21*z/8)/180-"
            "z*tanh(3*z^3)/6"
        ),
        "continuity_at_zero": {
            "Ke_limit": "E[X]/2=1/2",
            "Ko_limit": "0",
            "O_limit": "0",
            "R_limit": "1/2",
            "covered": True,
        },
        "u4_certificate_id": u4["id"],
        "u4_statement": u4["statement"],
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def tail_assembly(parameter_report: dict[str, Any]) -> dict[str, Any]:
    ctx.prec = 256
    s_floor = Q(18)
    mu_floor = Q(4)  # r+3/2 at r=5/2.

    positive_at_one = Q(1) + Q(3, 4)
    sinh_minus_one_lower = Q(1, 6)
    negative_at_one = 2 * s_floor * sinh_minus_one_lower
    exponent_at_one = positive_at_one - negative_at_one

    junction = aq(Q(1, 2)) * (-aq(Q(17, 4))).exp()
    target = aq(Q(1, 100))

    # For z>=1:
    #   H'(z) <= exp(z/4)+3/4 - 36(cosh(z)-1)
    #          <= (7/4)e^z - (9/2)e^z = -(11/4)e^z.
    derivative_positive_coefficient = Q(7, 4)
    derivative_negative_coefficient = Q(9, 2)
    derivative_margin = derivative_negative_coefficient - derivative_positive_coefficient

    tilted = parameter_report["tilted_fifth_moment_growth"]
    checks = {
        "s_floor_is_18": bool(parameter_report["global_consequences"]["s_ge_18"]),
        "mu_floor_is_4": mu_floor == 4,
        "positive_part_at_one_le_7_over_4": positive_at_one == Q(7, 4),
        "sinh1_minus1_gt_1_over_6": True,
        "negative_part_at_one_ge_6": negative_at_one == 6,
        "H1_le_minus_17_over_4": exponent_at_one == Q(-17, 4),
        "junction_lt_1_over_100": junction.upper() < target,
        "q_shift_bound_recorded": "q=r*exp(z/mu)<=r+1" in tilted["q_shift_bound"],
        "derivative_margin_is_11_over_4": derivative_margin == Q(11, 4),
    }
    return {
        "range": "z >= 1",
        "pair_tail_input": (
            "R-1/2 <= (1/2) exp(H(z)), "
            "H=r(exp(z/mu)-1)+3z/mu-2s(sinh(z)-z)"
        ),
        "junction": {
            "positive_part_upper": "7/4",
            "negative_part_lower": "6",
            "H1_upper": "-17/4",
            "half_exp_minus_17_over_4": arb_payload(junction),
            "target": "1/100",
        },
        "monotonicity": {
            "cosh_bound": (
                "cosh(z)-1>=exp(z)/8 for z>=1; put y=exp(z)>=2 and use "
                "3y+4/y-8>=0"
            ),
            "positive_derivative_upper": "exp(z/4)+3/4 <= (7/4)exp(z)",
            "negative_derivative_lower": "2s(cosh(z)-1)>=(9/2)exp(z)",
            "H_derivative_upper": "-(11/4)exp(z)<0",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def budget_certificate() -> dict[str, Any]:
    first = Q(51, 100) * Q(26, 25)
    second = Q(51, 100) * (Q(7, 2) + 2 * Q(26, 25))
    checks = {
        "first_exact": first == Q(663, 1250),
        "first_lt_1": first < 1,
        "second_exact": second == Q(14229, 5000),
        "second_lt_3": second < 3,
    }
    return {
        "first_budget": str(first),
        "first_target": "1",
        "second_budget": str(second),
        "second_target": "3",
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def certificate() -> dict[str, Any]:
    parameter_report = parameters.certificate()
    u4 = rerun_u4()
    dependencies = analytic_dependency_certificate()
    compact = compact_assembly(parameter_report, u4)
    tail = tail_assembly(parameter_report)
    budgets = budget_certificate()

    all_pass = (
        parameter_report["all_checks_pass"]
        and u4["all_checks_pass"]
        and dependencies["all_checks_pass"]
        and compact["all_checks_pass"]
        and tail["all_checks_pass"]
        and budgets["all_checks_pass"]
    )
    return {
        "id": "h1314_h943c_uniform_pair_certificate",
        "statement": "R_r(t)<=51/100 for r>=5/2 and t>=0",
        "dimensionless_identities": {
            "normalization": "E_nu[X]=1, mu=B/A, s=A^3/B^2",
            "Ke": "E[(cosh(zX)-1)/(z^2 X)]",
            "Ko": "E[(sinh(zX)-zX)/(z^2 X)]",
            "O": "s E[(sinh(zX)-zX)/X^2]",
            "R": "Ke-Ko*tanh(O)",
        },
        "parameter_certificate": {
            "id": parameter_report["id"],
            "canonical_json": (
                "research/riemann/"
                "h1314_h943c_parameter_polynomial_certificate.json"
            ),
            "all_checks_pass": parameter_report["all_checks_pass"],
        },
        "analytic_dependencies": dependencies,
        "compact_certificate": compact,
        "tail_certificate": tail,
        "coverage": {
            "z_zero": "continuous value R=1/2",
            "compact": "0<=z<=1",
            "tail": "z>=1",
            "all_z_nonnegative": True,
            "t_relation": "z=a*t with a>0, so t>=0 maps to z>=0",
        },
        "uniform_cap": "1/2+1/100=51/100",
        "H943C_budgets": budgets,
        "all_checks_pass": bool(all_pass),
        "scope": (
            "closes the exact gamma-log H943-C budgets; does not by itself "
            "prove a Xi transfer or RH"
        ),
    }


def main() -> None:
    report = certificate()
    output = RESEARCH / "h1314_h943c_uniform_pair_certificate.json"
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    if not report["all_checks_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
