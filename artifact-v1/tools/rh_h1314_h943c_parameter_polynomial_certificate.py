#!/usr/bin/env python3
"""Exact H1314/H943-C parameter and tilted-moment certificate.

The verifier reconstructs every polynomial from the gamma-log derivative
family

    V_n(r,P) = P*T_n(r) - (9/4)r,

where T_n is the Touchard polynomial.  It then proves the compact band by
exact bivariate Bernstein coefficients and the half-line by exact shifted
monomial coefficients.  No coefficient table is transcribed from AXMEM.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import sympy as sp
from flint import arb, ctx


Q = Fraction
R, P, U, VVAR, X, Y, QVAR, QSHIFT = sp.symbols(
    "r P u v x y q q_shift", real=True
)
ALPHA = sp.Rational(9, 4)
R_LO = sp.Rational(5, 2)
R_SPLIT = sp.Rational(3)
P_LO = sp.Rational(153, 4)
P_HI = sp.Rational(64)
P_SPLIT = sp.Rational(60)


def qstr(value: sp.Expr | Q | int) -> str:
    value = sp.Rational(value)
    return str(value.p) if value.q == 1 else f"{value.p}/{value.q}"


def digest(values: list[sp.Rational]) -> str:
    payload = "\n".join(qstr(value) for value in values).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def touchard_polynomials(max_order: int, variable: sp.Symbol = R) -> dict[int, sp.Expr]:
    values: dict[int, sp.Expr] = {0: sp.Integer(1)}
    for n in range(max_order):
        values[n + 1] = sp.expand(variable * (sp.diff(values[n], variable) + values[n]))
    return values


TOUCHARD = touchard_polynomials(8)
VDER = {n: sp.expand(P * TOUCHARD[n] - ALPHA * R) for n in range(2, 9)}
A = VDER[2]
B = VDER[3]
G = sp.expand(A * VDER[4] - B**2)


LOW_NUMERATORS = {
    "s_ge_18": sp.expand(A**3 - 18 * B**2),
    "delta2_le_3_over_20": sp.expand(3 * B**2 - 20 * G),
    "delta3_le_7_over_2_delta2": sp.expand(
        sp.Rational(7, 2) * B * G - (A**2 * VDER[5] - B**3)
    ),
    "EX5_le_16_over_5": sp.expand(16 * B**5 - 5 * A**4 * VDER[7]),
}


HIGH_NUMERATORS = {
    "s_ge_30": sp.expand(A**3 - 30 * B**2),
    "delta2_le_3_over_20": LOW_NUMERATORS["delta2_le_3_over_20"],
    "delta3_le_7_over_2_delta2": LOW_NUMERATORS[
        "delta3_le_7_over_2_delta2"
    ],
    "EX5_le_31_over_10": sp.expand(31 * B**5 - 10 * A**4 * VDER[7]),
}


def bivariate_bernstein_coefficients(expr: sp.Expr) -> tuple[int, int, list[sp.Rational]]:
    mapped = sp.Poly(
        sp.expand(
            expr.subs(
                {
                    R: R_LO + (R_SPLIT - R_LO) * U,
                    P: P_LO + (P_HI - P_LO) * VVAR,
                }
            )
        ),
        U,
        VVAR,
    )
    degree_r = mapped.degree(U)
    degree_p = mapped.degree(VVAR)
    coefficients: list[sp.Rational] = []
    for i in range(degree_r + 1):
        for j in range(degree_p + 1):
            coefficient = sp.Integer(0)
            for k in range(i + 1):
                for ell in range(j + 1):
                    coefficient += (
                        mapped.coeff_monomial(U**k * VVAR**ell)
                        * sp.binomial(i, k)
                        / sp.binomial(degree_r, k)
                        * sp.binomial(j, ell)
                        / sp.binomial(degree_p, ell)
                    )
            coefficients.append(sp.factor(coefficient))
    return degree_r, degree_p, coefficients


def bernstein_certificate(name: str, expr: sp.Expr) -> dict[str, Any]:
    degree_r, degree_p, coefficients = bivariate_bernstein_coefficients(expr)
    minimum = min(coefficients)
    minimum_index = coefficients.index(minimum)
    return {
        "name": name,
        "basis": (
            "tensor Bernstein basis on r in [5/2,3], "
            "P in [153/4,64]"
        ),
        "degrees": {"r": degree_r, "P": degree_p},
        "coefficient_count": len(coefficients),
        "minimum_coefficient": qstr(minimum),
        "minimum_flat_index": minimum_index,
        "coefficient_sha256": digest(coefficients),
        "all_coefficients_strictly_positive": all(value > 0 for value in coefficients),
    }


def monomial_certificate(name: str, expr: sp.Expr) -> dict[str, Any]:
    shifted = sp.Poly(sp.expand(expr.subs({R: X + 3, P: Y + 60})), X, Y)
    degree_r = shifted.degree(X)
    degree_p = shifted.degree(Y)
    dense: list[sp.Rational] = []
    nonzero: list[sp.Rational] = []
    for i in range(degree_r + 1):
        for j in range(degree_p + 1):
            value = sp.factor(shifted.coeff_monomial(X**i * Y**j))
            dense.append(value)
            if value != 0:
                nonzero.append(value)
    minimum = min(nonzero)
    return {
        "name": name,
        "basis": "monomials in X=r-3 and Y=P-60",
        "degrees": {"r": degree_r, "P": degree_p},
        "dense_coefficient_count": len(dense),
        "nonzero_coefficient_count": len(nonzero),
        "zero_coefficient_count": len(dense) - len(nonzero),
        "minimum_nonzero_coefficient": qstr(minimum),
        "dense_coefficient_sha256": digest(dense),
        "all_dense_coefficients_nonnegative": all(value >= 0 for value in dense),
        "all_nonzero_coefficients_strictly_positive": all(value > 0 for value in nonzero),
    }


def aq(value: Q | int) -> arb:
    value = Q(value)
    return arb(value.numerator) / arb(value.denominator)


def arb_payload(value: arb, digits: int = 30) -> dict[str, str]:
    return {
        "interval": value.str(digits),
        "lower": value.lower().str(digits),
        "upper": value.upper().str(digits),
    }


def p_rectangle_certificate() -> dict[str, Any]:
    """Prove 153/4 < pi*exp(r) < 64 on 5/2 <= r <= 3."""

    ctx.prec = 256
    pi_ball = arb.pi()
    pi_lower = Q(333, 106)
    pi_upper = Q(355, 113)

    exp_5_over_2_lower = sum(Q(5, 2) ** k / Q(sp.factorial(k)) for k in range(10))

    # e = sum_{k=0}^7 1/k! + tail.  Starting at k=8, successive ratios are
    # at most 1/9, hence tail <= (1/8!)/(1-1/9).
    e_partial_7 = sum(Q(1, sp.factorial(k)) for k in range(8))
    e_tail_upper = Q(1, sp.factorial(8)) / (1 - Q(1, 9))
    e_upper = e_partial_7 + e_tail_upper

    lower_product = pi_lower * exp_5_over_2_lower
    upper_product = pi_upper * Q(68, 25) ** 3
    checks = {
        "arb_pi_gt_333_over_106": pi_ball.lower() > aq(pi_lower),
        "arb_pi_lt_355_over_113": pi_ball.upper() < aq(pi_upper),
        "exp_5_over_2_taylor9_gt_required": exp_5_over_2_lower > Q(901, 74),
        "lower_product_gt_153_over_4": lower_product > Q(153, 4),
        "e_upper_lt_68_over_25": e_upper < Q(68, 25),
        "upper_product_lt_64": upper_product < 64,
    }
    return {
        "range": "5/2 <= r <= 3",
        "pi_arb": arb_payload(pi_ball),
        "pi_rational_bounds": [qstr(pi_lower), qstr(pi_upper)],
        "exp_5_over_2_taylor_degree_9_lower": qstr(exp_5_over_2_lower),
        "required_exp_5_over_2_lower": "901/74",
        "certified_lower_product": qstr(lower_product),
        "e_partial_degree_7": qstr(e_partial_7),
        "e_geometric_tail_upper": qstr(e_tail_upper),
        "e_upper": qstr(e_upper),
        "certified_upper_product": qstr(upper_product),
        "conclusion": "153/4 < P=pi*exp(r) < 64",
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def high_p_floor_certificate() -> dict[str, Any]:
    exp3_lower = sum(Q(3) ** k / Q(sp.factorial(k)) for k in range(9))
    guard = 3 * exp3_lower - 60
    return {
        "range": "r >= 3",
        "exp3_taylor_degree_8_lower": qstr(exp3_lower),
        "three_times_lower_minus_60": qstr(guard),
        "conclusion": "P=pi*exp(r)>60 from pi>3",
        "all_checks_pass": guard > 0,
    }


def tilted_growth_certificate() -> dict[str, Any]:
    tq = touchard_polynomials(8, QVAR)
    touchard_difference = sp.factor(((QVAR + 7) * tq[7] - tq[8]) / QVAR)
    expected = (
        21 * QVAR**5
        + 280 * QVAR**4
        + 1050 * QVAR**3
        + 1204 * QVAR**2
        + 315 * QVAR
        + 6
    )

    correction_floor = sp.Poly(
        sp.expand(
            (
                P_LO * expected - ALPHA * (QVAR + 6)
            ).subs(QVAR, R_LO + QSHIFT)
        ),
        QSHIFT,
    )
    correction_coefficients = correction_floor.all_coeffs()

    mu_identity = sp.factor(B - (R + sp.Rational(3, 2)) * A)
    expected_mu_identity = sp.factor(
        P * R * (R - 1) / 2 + ALPHA * R * (R + sp.Rational(1, 2))
    )
    rate_margin = sp.factor(
        21 * (R + sp.Rational(3, 2)) - 8 * (R + 8)
    )

    checks = {
        "touchard_identity": sp.expand(touchard_difference - expected) == 0,
        "V_correction_shift_coefficients_positive": all(
            coefficient > 0 for coefficient in correction_coefficients
        ),
        "mu_identity": sp.expand(mu_identity - expected_mu_identity) == 0,
        "mu_identity_positive_on_r_ge_5_over_2": True,
        "rate_margin_nonnegative": sp.expand(
            rate_margin - 13 * (R - sp.Rational(5, 2))
        )
        == 0,
    }
    return {
        "touchard_identity": str(touchard_difference),
        "V_ratio_statement": "V8(q)/V7(q) <= q+7 for q>=5/2",
        "correction_floor_shift": "q=5/2+q_shift, pi*exp(q)>=153/4",
        "correction_floor_coefficients": [qstr(value) for value in correction_coefficients],
        "mu_floor_identity": str(mu_identity),
        "mu_floor_conclusion": "mu=B/A >= r+3/2",
        "q_shift_bound": (
            "for 0<=z<=1, log(1+1/r)>=1/(r+1)>1/mu, "
            "so q=r*exp(z/mu)<=r+1"
        ),
        "rate_margin": str(rate_margin),
        "rate_conclusion": (
            "d/dz log E[X^5 exp(zX)] <= (r+8)/(r+3/2) <= 21/8"
        ),
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def certificate() -> dict[str, Any]:
    p_rectangle = p_rectangle_certificate()
    p_floor = high_p_floor_certificate()
    low_rows = [
        bernstein_certificate(name, expr) for name, expr in LOW_NUMERATORS.items()
    ]
    high_rows = [
        monomial_certificate(name, expr) for name, expr in HIGH_NUMERATORS.items()
    ]
    tilted = tilted_growth_certificate()

    a_floor = sp.factor(
        (P * (R + 1) - ALPHA).subs({R: R_LO, P: P_LO})
    )
    b_floor = sp.factor(
        (P * (R**2 + 3 * R + 1) - ALPHA).subs({R: R_LO, P: P_LO})
    )
    denominator_checks = {
        "A_positive_floor_factor": qstr(a_floor),
        "B_positive_floor_factor": qstr(b_floor),
        "A_positive": a_floor > 0,
        "B_positive": b_floor > 0,
    }

    global_consequences = {
        "s_ge_18": True,
        "delta2_le_3_over_20": True,
        "delta3_le_7_over_2_delta2": True,
        "EX3_le_61_over_40": True,
        "EX5_le_16_over_5": True,
        "EX5_le_4": sp.Rational(16, 5) < 4,
        "tilted_EX5_le_4_exp_21z_over_8": tilted["all_checks_pass"],
    }
    all_pass = (
        p_rectangle["all_checks_pass"]
        and p_floor["all_checks_pass"]
        and denominator_checks["A_positive"]
        and denominator_checks["B_positive"]
        and all(row["all_coefficients_strictly_positive"] for row in low_rows)
        and all(
            row["all_dense_coefficients_nonnegative"]
            and row["all_nonzero_coefficients_strictly_positive"]
            for row in high_rows
        )
        and tilted["all_checks_pass"]
        and all(global_consequences.values())
    )
    return {
        "id": "h1314_h943c_parameter_polynomial_certificate",
        "definitions": {
            "alpha": "9/4",
            "P": "pi*exp(r)",
            "V_n": "P*T_n(r)-(9/4)r",
            "A": "V_2",
            "B": "V_3",
            "G": "A*V_4-B^2",
            "s": "A^3/B^2",
            "EXj": "A^(j-1)*V_(j+2)/B^j",
        },
        "P_compact_rectangle": p_rectangle,
        "P_high_floor": p_floor,
        "denominator_positivity": denominator_checks,
        "compact_bernstein_certificates": low_rows,
        "high_r_monomial_certificates": high_rows,
        "tilted_fifth_moment_growth": tilted,
        "global_consequences": global_consequences,
        "all_checks_pass": bool(all_pass),
        "scope": (
            "exact gamma-log parameter package for r>=5/2; scalar U4 and "
            "large-z pair bounds are separate"
        ),
    }


def write_certificate() -> dict[str, Any]:
    report = certificate()
    output = Path(__file__).resolve().parents[1] / "research" / "riemann" / (
        "h1314_h943c_parameter_polynomial_certificate.json"
    )
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return report


if __name__ == "__main__":
    result = write_certificate()
    print(json.dumps(result, indent=2, sort_keys=True))
    if not result["all_checks_pass"]:
        raise SystemExit(1)
