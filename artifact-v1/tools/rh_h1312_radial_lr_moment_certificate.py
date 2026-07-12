"""H1312 exact moment and one-dimensional radial-LR certificate.

This verifier proves the coarse parameter inequalities used by the H1312
radial likelihood-ratio comparison on

    3 <= r <= 49/5.

The algebraic part is exact over QQ.  Three numerators become polynomials
with positive monomial coefficients after

    X = r-3,  Y = pi*exp(r)-60.

The two inequalities that need the upper endpoint r<=49/5 are checked by
positive exact Bernstein coefficients in r, coefficient by coefficient in
Y.  The final scalar envelope in z is enclosed by Arb on rational boxes.
"""

from __future__ import annotations

import json
from fractions import Fraction
from typing import Any

import sympy as sp
from flint import arb, ctx


R, P, X, Y, U, QVAR = sp.symbols("r P X Y u q", real=True)
ALPHA = sp.Rational(9, 4)
R_LO = sp.Rational(3)
R_HI = sp.Rational(49, 5)
P_FLOOR = sp.Rational(60)


def touchard_polynomials(max_order: int) -> dict[int, sp.Expr]:
    values: dict[int, sp.Expr] = {0: sp.Integer(1)}
    for n in range(max_order):
        values[n + 1] = sp.expand(R * (sp.diff(values[n], R) + values[n]))
    return values


TOUCHARD = touchard_polynomials(8)
V = {n: sp.expand(P * TOUCHARD[n] - ALPHA * R) for n in range(2, 8)}
A = V[2]
B = V[3]
G = sp.expand(A * V[4] - B**2)


def bernstein_coefficients(poly: sp.Expr) -> list[sp.Rational]:
    """Exact Bernstein coefficients of poly(r) on [3,49/5]."""

    mapped = sp.Poly(
        sp.expand(poly.subs(R, R_LO + (R_HI - R_LO) * U)), U
    )
    degree = mapped.degree()
    power = [mapped.nth(k) for k in range(degree + 1)]
    return [
        sp.factor(
            sum(
                power[k]
                * sp.binomial(i, k)
                / sp.binomial(degree, k)
                for k in range(i + 1)
            )
        )
        for i in range(degree + 1)
    ]


def positive_monomial_certificate(name: str, expr: sp.Expr) -> dict[str, Any]:
    shifted = sp.Poly(
        sp.expand(expr.subs({R: X + R_LO, P: Y + P_FLOOR})), X, Y
    )
    coefficients = shifted.coeffs()
    return {
        "name": name,
        "basis": "monomial in X=r-3 and Y=P-60",
        "coefficient_count": len(coefficients),
        "minimum_coefficient": str(min(coefficients)),
        "all_coefficients_positive": all(c > 0 for c in coefficients),
        "total_degree": shifted.total_degree(),
    }


def positive_bernstein_in_y_certificate(
    name: str, expr: sp.Expr
) -> dict[str, Any]:
    shifted = sp.Poly(sp.expand(expr.subs(P, Y + P_FLOOR)), Y)
    rows: list[dict[str, Any]] = []
    for power in range(shifted.degree() + 1):
        coefficients = bernstein_coefficients(shifted.nth(power))
        rows.append(
            {
                "Y_power": power,
                "bernstein_degree": len(coefficients) - 1,
                "minimum_coefficient": str(min(coefficients)),
                "all_coefficients_positive": all(c > 0 for c in coefficients),
            }
        )
    return {
        "name": name,
        "basis": "powers of Y=P-60 with Bernstein coefficients in r",
        "Y_degree": shifted.degree(),
        "rows": rows,
        "all_coefficients_positive": all(
            row["all_coefficients_positive"] for row in rows
        ),
    }


def aq(value: Fraction | int) -> arb:
    value = Fraction(value)
    return arb(value.numerator) / arb(value.denominator)


def abox(left: Fraction, right: Fraction) -> arb:
    lo = aq(left)
    hi = aq(right)
    return (lo + hi) / 2 + arb(0, 1) * (hi - lo) / 2


def arb_payload(value: arb, digits: int = 24) -> dict[str, str]:
    return {
        "interval": value.str(digits),
        "lower": value.lower().str(digits),
        "upper": value.upper().str(digits),
    }


def scalar_envelope_certificate(subdivisions: int = 64) -> dict[str, Any]:
    """Certify the explicit z-envelope below 199/250 on [0,1]."""

    ctx.prec = 192
    split = Fraction(13, 25)
    target = aq(Fraction(199, 250))
    log4 = arb(4).log()
    constant = 4 * (arb(1).cosh() - 1)
    worst: tuple[arb, Fraction, Fraction, str, arb] | None = None
    failures: list[dict[str, str]] = []

    branches = [
        (Fraction(0), split, "constant"),
        (split, Fraction(1), "tail"),
    ]
    for start, stop, mode in branches:
        for index in range(subdivisions):
            left = start + (stop - start) * index / subdivisions
            right = start + (stop - start) * (index + 1) / subdivisions
            z = abox(left, right)
            q_poly = (
                arb(1) / 6
                - arb(7) * z / 48
                + z**2 / 24
                - z**3 / 24
            )
            k_majorant = arb(7) / 240 + z**2 * (arb(5) * z / 2).exp() / 168
            if mode == "constant":
                m_majorant = arb(3) * log4 / 4
            else:
                m_majorant = 30 * z**3 * (-10 * z**3).exp()
            ratio = constant * z * k_majorant * m_majorant / q_poly
            if q_poly.lower() <= arb(0) or ratio.upper() >= target:
                failures.append(
                    {
                        "left": str(left),
                        "right": str(right),
                        "mode": mode,
                        "Q": q_poly.str(20),
                        "ratio": ratio.str(20),
                    }
                )
            if worst is None or ratio.upper() > worst[0]:
                worst = (ratio.upper(), left, right, mode, ratio)

    assert worst is not None
    split_guard = 30 * aq(split) ** 3 - 3 * log4
    tail_floor = aq(Fraction(1, 960))
    tail_rhs_at_one = (
        2
        * (arb(1).cosh() - 1)
        * (-60 * (arb(1).sinh() - 1)).exp()
    )
    tail_log_derivative_at_one = (
        arb(1).sinh() / (arb(1).cosh() - 1)
        - 60 * (arb(1).cosh() - 1)
    )
    tail_pass = (
        split_guard.lower() > arb(0)
        and tail_rhs_at_one.upper() < tail_floor
        and tail_log_derivative_at_one.upper() < arb(0)
    )
    return {
        "subdivisions_per_branch": subdivisions,
        "target": "199/250",
        "split": "13/25",
        "split_guard_30z3_minus_3log4": arb_payload(split_guard),
        "worst_box": {
            "left": str(worst[1]),
            "right": str(worst[2]),
            "mode": worst[3],
            "ratio": arb_payload(worst[4]),
            "upper": worst[0].str(24),
        },
        "failure_count": len(failures),
        "failures": failures[:10],
        "compact_pass": not failures,
        "tail_l_floor": "1/960",
        "tail_rhs_at_z1": arb_payload(tail_rhs_at_one),
        "tail_log_derivative_at_z1": arb_payload(tail_log_derivative_at_one),
        "tail_pass": bool(tail_pass),
    }


def certificate() -> dict[str, Any]:
    numerators = {
        "s_ge_30": sp.expand(A**3 - 30 * B**2),
        "delta3_over_delta2_le_7_over_2": sp.expand(
            sp.Rational(7, 2) * B * G - (A**2 * V[5] - B**3)
        ),
        "delta4_over_delta2_ge_5": sp.expand(
            A**3 * V[6] - B**4 - 5 * B**2 * G
        ),
    }
    monomial_rows = [
        positive_monomial_certificate(name, expr)
        for name, expr in numerators.items()
    ]

    endpoint_numerators = {
        "delta2_ge_1_over_20": sp.expand(20 * G - B**2),
        "EX5_over_delta2_le_30": sp.expand(30 * B**3 * G - A**4 * V[7]),
    }
    bernstein_rows = [
        positive_bernstein_in_y_certificate(name, expr)
        for name, expr in endpoint_numerators.items()
    ]

    exp3_floor = sum(sp.Rational(3) ** n / sp.factorial(n) for n in range(9))
    p_floor_guard = sp.expand(3 * exp3_floor - 60)

    mu_guard_identity = sp.factor(B - (R + sp.Rational(3, 2)) * A)
    expected_mu_guard = sp.factor(
        P * R * (R - 1) / 2 + ALPHA * R * (R + sp.Rational(1, 2))
    )

    tq = {n: sp.expand(TOUCHARD[n].subs(R, QVAR)) for n in range(9)}
    touchard_ratio_identity = sp.factor(
        ((QVAR + 6) * tq[7] - tq[8]) / QVAR
    )
    expected_touchard_ratio = sp.expand(
        (140 - QVAR**2) * QVAR**4
        + 700 * QVAR**3
        + 903 * QVAR**2
        + 252 * QVAR
        + 5
    )
    q_square_guard = sp.Rational(140) - sp.Rational(54, 5) ** 2
    tilted_mean_guard = sp.factor(
        sp.Rational(5, 2) * (R + sp.Rational(3, 2)) - (R + 7)
    )

    scalar = scalar_envelope_certificate()
    all_pass = (
        p_floor_guard > 0
        and all(row["all_coefficients_positive"] for row in monomial_rows)
        and all(row["all_coefficients_positive"] for row in bernstein_rows)
        and sp.expand(mu_guard_identity - expected_mu_guard) == 0
        and sp.expand(touchard_ratio_identity - expected_touchard_ratio) == 0
        and q_square_guard > 0
        and tilted_mean_guard.subs(R, R_LO) > 0
        and scalar["compact_pass"]
        and scalar["tail_pass"]
    )
    return {
        "id": "h1312_radial_lr_moment_certificate",
        "range": "3 <= r <= 49/5",
        "P_definition": "P=pi*exp(r)",
        "P_floor": {
            "exp3_taylor_degree_8": str(exp3_floor),
            "three_times_floor_minus_60": str(p_floor_guard),
            "conclusion": "P>60 from pi>3 and exp(r)>=exp(3)",
        },
        "monomial_certificates": monomial_rows,
        "bernstein_certificates": bernstein_rows,
        "tilted_moment_growth": {
            "mu_guard_identity": str(mu_guard_identity),
            "identity_matches": bool(
                sp.expand(mu_guard_identity - expected_mu_guard) == 0
            ),
            "touchard_ratio_identity": str(touchard_ratio_identity),
            "touchard_identity_matches": bool(
                sp.expand(touchard_ratio_identity - expected_touchard_ratio) == 0
            ),
            "140_minus_54_over_5_squared": str(q_square_guard),
            "five_halves_mu_guard_from_mu_floor": str(tilted_mean_guard),
            "conclusion": "d/dz log E[X^5 exp(zX)] <= 5/2 on 0<=z<=1",
        },
        "scalar_envelope": scalar,
        "all_checks_pass": bool(all_pass),
        "limitation": "proves radial LR only for the exact gamma-log family on 3<=r<=49/5",
    }


if __name__ == "__main__":
    print(json.dumps(certificate(), indent=2))
