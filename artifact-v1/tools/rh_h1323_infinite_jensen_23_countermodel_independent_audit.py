#!/usr/bin/env python3
"""Independent stdlib/Fraction audit of the H1323 infinite countermodel."""

from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import TypeAlias


ROOT = Path(__file__).resolve().parents[1]
PRIMARY_PATH = ROOT / "research" / "riemann" / "h1323_infinite_jensen_23_countermodel.json"
AUDIT_PATH = (
    ROOT
    / "research"
    / "riemann"
    / "h1323_infinite_jensen_23_countermodel_independent_audit.json"
)
MAX_N = 256

Polynomial: TypeAlias = tuple[Fraction, ...]
RationalFunction: TypeAlias = tuple[Polynomial, Polynomial]


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def polynomial(*coefficients: int | Fraction) -> Polynomial:
    values = [Fraction(value) for value in coefficients]
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values)


def p_add(left: Polynomial, right: Polynomial) -> Polynomial:
    size = max(len(left), len(right))
    return polynomial(
        *(sum((values[index] if index < len(values) else Fraction(0)) for values in (left, right)) for index in range(size))
    )


def p_neg(value: Polynomial) -> Polynomial:
    return polynomial(*(-coefficient for coefficient in value))


def p_sub(left: Polynomial, right: Polynomial) -> Polynomial:
    return p_add(left, p_neg(right))


def p_mul(left: Polynomial, right: Polynomial) -> Polynomial:
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, left_coefficient in enumerate(left):
        for j, right_coefficient in enumerate(right):
            result[i + j] += left_coefficient * right_coefficient
    return polynomial(*result)


def p_scale(value: Polynomial, scalar: int | Fraction) -> Polynomial:
    return polynomial(*(Fraction(scalar) * coefficient for coefficient in value))


def p_pow(value: Polynomial, exponent: int) -> Polynomial:
    result = polynomial(1)
    for _ in range(exponent):
        result = p_mul(result, value)
    return result


def rf(numerator: Polynomial, denominator: Polynomial = polynomial(1)) -> RationalFunction:
    return numerator, denominator


def r_const(value: int | Fraction) -> RationalFunction:
    return rf(polynomial(value))


def r_add(left: RationalFunction, right: RationalFunction) -> RationalFunction:
    return rf(
        p_add(p_mul(left[0], right[1]), p_mul(right[0], left[1])),
        p_mul(left[1], right[1]),
    )


def r_neg(value: RationalFunction) -> RationalFunction:
    return rf(p_neg(value[0]), value[1])


def r_sub(left: RationalFunction, right: RationalFunction) -> RationalFunction:
    return r_add(left, r_neg(right))


def r_mul(left: RationalFunction, right: RationalFunction) -> RationalFunction:
    return rf(p_mul(left[0], right[0]), p_mul(left[1], right[1]))


def r_div(left: RationalFunction, right: RationalFunction) -> RationalFunction:
    return rf(p_mul(left[0], right[1]), p_mul(left[1], right[0]))


def r_pow(value: RationalFunction, exponent: int) -> RationalFunction:
    return rf(p_pow(value[0], exponent), p_pow(value[1], exponent))


def r_equal(left: RationalFunction, right: RationalFunction) -> bool:
    return p_mul(left[0], right[1]) == p_mul(right[0], left[1])


def tail_q(shift: int) -> RationalFunction:
    # (201/208)*(2*(n+shift)-3)/(2*(n+shift))
    return rf(
        polynomial(201 * (2 * shift - 3), 402),
        polynomial(416 * shift, 416),
    )


def shifted_n(shift: int) -> Polynomial:
    return polynomial(shift, 1)


def capital_C(shift: int) -> RationalFunction:
    z = shifted_n(shift)
    return rf(p_mul(p_scale(z, 2), p_sub(p_scale(z, 2), polynomial(1))))


PREFIX = {
    2: Fraction(13, 32),
    3: Fraction(141, 256),
    4: Fraction(85, 128),
    5: Fraction(89, 128),
    6: Fraction(193, 256),
    7: Fraction(201, 256),
    8: Fraction(201, 256),
}
C_LIMIT = Fraction(201, 208)


def q(index: int) -> Fraction:
    if index in PREFIX:
        return PREFIX[index]
    if index >= 8:
        return C_LIMIT * Fraction(2 * index - 3, 2 * index)
    raise ValueError(index)


def cubic_discriminant(leading: Fraction, quadratic: Fraction, linear: Fraction, constant: Fraction) -> Fraction:
    return (
        quadratic**2 * linear**2
        - 4 * leading * linear**3
        - 4 * quadratic**3 * constant
        - 27 * leading**2 * constant**2
        + 18 * leading * quadratic * linear * constant
    )


def quartic_discriminant(a: Fraction, b: Fraction, c: Fraction, d: Fraction, e: Fraction) -> Fraction:
    return (
        256 * a**3 * e**3
        - 192 * a**2 * b * d * e**2
        - 128 * a**2 * c**2 * e**2
        + 144 * a**2 * c * d**2 * e
        - 27 * a**2 * d**4
        + 144 * a * b**2 * c * e**2
        - 6 * a * b**2 * d**2 * e
        - 80 * a * b * c**2 * d * e
        + 18 * a * b * c * d**3
        + 16 * a * c**4 * e
        - 4 * a * c**3 * d**2
        - 27 * b**4 * e**2
        + 18 * b**3 * c * d * e
        - 4 * b**3 * d**3
        - 4 * b**2 * c**3 * e
        + b**2 * c**2 * d**2
    )


def build_coefficients(max_index: int) -> tuple[list[Fraction], list[Fraction], list[Fraction]]:
    ratios = [Fraction(0)] * (max_index + 1)
    coefficients = [Fraction(0)] * (max_index + 1)
    ratios[1] = coefficients[0] = coefficients[1] = Fraction(1)
    for index in range(2, max_index + 1):
        ratios[index] = ratios[index - 1] * q(index)
        coefficients[index] = coefficients[index - 1] * ratios[index]
    gamma = [Fraction(math.factorial(index)) * coefficients[index] for index in range(max_index + 1)]
    return ratios, coefficients, gamma


def build_audit() -> dict[str, object]:
    primary = json.loads(PRIMARY_PATH.read_text(encoding="utf-8"))
    ratios, coefficients, gamma = build_coefficients(20)

    q_delta = r_sub(tail_q(1), tail_q(0))
    q_delta_expected = rf(polynomial(603), polynomial(0, 416, 416))

    tau_ratio = r_mul(
        r_div(tail_q(1), tail_q(0)),
        r_div(r_mul(capital_C(1), capital_C(-1)), r_pow(capital_C(0), 2)),
    )
    tau_expected = rf(
        p_mul(polynomial(-1, 1), polynomial(1, 2)),
        p_mul(polynomial(0, 1), polynomial(-1, 2)),
    )
    unit_barrier = r_sub(r_const(1), rf(polynomial(1), polynomial(0, 0, 1)))
    unit_margin = r_sub(tau_ratio, unit_barrier)
    margin_expected = rf(
        polynomial(-1, 1),
        p_mul(polynomial(0, 0, 1), polynomial(-1, 2)),
    )

    one = r_const(1)
    D_derived = r_sub(
        r_pow(r_sub(one, tail_q(1)), 2),
        r_mul(
            r_pow(tail_q(1), 2),
            r_mul(r_sub(one, tail_q(0)), r_sub(one, tail_q(2))),
        ),
    )
    D_expected_num = polynomial(
        -15372297693,
        192552254264,
        13154707832,
        297357088,
        2244592,
    )
    D_expected_den = p_scale(
        p_mul(
            polynomial(0, 1),
            p_mul(p_pow(polynomial(1, 1), 2), polynomial(2, 1)),
        ),
        29948379136,
    )
    D_expected = rf(D_expected_num, D_expected_den)

    u = rf(polynomial(201, 402), polynomial(416, 416))
    v = rf(polynomial(603, 402), polynomial(832, 416))
    F_derived = r_add(
        r_sub(r_pow(r_mul(u, v), 2), r_mul(r_const(6), r_mul(u, v))),
        r_sub(r_add(r_mul(r_const(4), u), r_mul(r_const(4), v)), r_const(3)),
    )
    F_positive_poly = polynomial(5859746333, 7568113136, 1236652088, 73874752, 1509200)
    F_expected_num = p_scale(F_positive_poly, -3)
    F_expected_den = p_scale(
        p_mul(p_pow(polynomial(1, 1), 2), p_pow(polynomial(2, 1), 2)),
        29948379136,
    )
    F_expected = rf(F_expected_num, F_expected_den)

    prefix_D = []
    expected_prefix_D = [
        Fraction(37926823, 268435456),
        Fraction(28183907, 536870912),
        Fraction(28382139, 536870912),
        Fraction(50157087, 2147483648),
        Fraction(58256935, 4294967296),
        Fraction(142163945, 6979321856),
    ]
    for index in range(2, 8):
        prefix_D.append(
            (1 - q(index + 1)) ** 2
            - q(index + 1) ** 2 * (1 - q(index)) * (1 - q(index + 2))
        )

    degree2_prefix = []
    degree3_prefix = []
    for shift in range(6):
        u0 = Fraction(shift + 2, shift + 1) * q(shift + 2)
        v0 = Fraction(shift + 3, shift + 2) * q(shift + 3)
        degree2_prefix.append(4 * (1 - u0))
        degree3_prefix.append(
            cubic_discriminant(u0**2 * v0, 3 * u0, Fraction(3), Fraction(1))
        )

    j4_ascending = [Fraction(math.comb(4, index)) * gamma[index] for index in range(5)]
    j4_discriminant = quartic_discriminant(
        j4_ascending[4],
        j4_ascending[3],
        j4_ascending[2],
        j4_ascending[1],
        j4_ascending[0],
    )

    # Independent finite sweeps supplement, but do not replace, the exact
    # rational-function identities above.
    h804_sweep = []
    D_sweep = []
    jensen_sweep = []
    for index in range(8, MAX_N + 1):
        C_int = lambda z: 2 * z * (2 * z - 1)
        observed_tau = (
            q(index + 1)
            / q(index)
            * Fraction(C_int(index + 1) * C_int(index - 1), C_int(index) ** 2)
        )
        h804_sweep.append(observed_tau > 1 - Fraction(1, index**2))
        D_sweep.append(
            (1 - q(index + 1)) ** 2
            - q(index + 1) ** 2 * (1 - q(index)) * (1 - q(index + 2))
            > 0
        )
    for shift in range(MAX_N + 1):
        u0 = Fraction(shift + 2, shift + 1) * q(shift + 2)
        v0 = Fraction(shift + 3, shift + 2) * q(shift + 3)
        disc2 = 4 * (1 - u0)
        disc3 = cubic_discriminant(u0**2 * v0, 3 * u0, Fraction(3), Fraction(1))
        jensen_sweep.append(disc2 > 0 and disc3 > 0)

    toeplitz_pf2 = []
    indices = range(21)
    for rows in itertools.combinations(indices, 2):
        for columns in itertools.combinations(indices, 2):
            matrix = [
                [
                    Fraction(0) if column - row < 0 else coefficients[column - row]
                    for column in columns
                ]
                for row in rows
            ]
            toeplitz_pf2.append(matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0] >= 0)

    primary_checks = {
        "primary_all_checks_pass": primary["all_checks_pass"] is True,
        "primary_q_prefix_matches": primary["construction"]["q_2_to_q_8"]
        == [fraction_text(q(index)) for index in range(2, 9)],
        "primary_D_tail_matches": primary["h851_translated_D3"]["tail_formula"]
        == "(2244592*n**4 + 297357088*n**3 + 13154707832*n**2 + 192552254264*n - 15372297693)/(29948379136*n*(n + 1)**2*(n + 2))",
        "primary_F_tail_matches": primary["jensen_degree_2_3_all_shifts"]["F_tail_factorization"]
        == "-3*(1509200*n**4 + 73874752*n**3 + 1236652088*n**2 + 7568113136*n + 5859746333)/(29948379136*(n + 1)**2*(n + 2)**2)",
        "primary_J4_discriminant_matches": primary["degree4_falsifier"]["discriminant"]
        == fraction_text(j4_discriminant),
    }

    checks = {
        "prefix_q_nondecreasing_below_c": all(
            PREFIX[index] <= PREFIX[index + 1] for index in range(2, 8)
        )
        and all(value < C_LIMIT for value in PREFIX.values()),
        "tail_q_delta_symbolic_identity": r_equal(q_delta, q_delta_expected),
        "tail_q_numerical_sweep_n8_to_n256": all(
            q(index) <= q(index + 1) < C_LIMIT for index in range(8, MAX_N + 1)
        ),
        "PF2_all_minors_on_indices_0_to_20": all(toeplitz_pf2),
        "h804_tau_symbolic_identity": r_equal(tau_ratio, tau_expected),
        "h804_margin_symbolic_identity": r_equal(unit_margin, margin_expected),
        "h804_sweep_n8_to_n256": all(h804_sweep),
        "D_prefix_exact_and_positive": prefix_D == expected_prefix_D
        and all(value > 0 for value in prefix_D),
        "D_tail_symbolic_identity": r_equal(D_derived, D_expected),
        "D_tail_positive_coefficient_certificate": all(value > 0 for value in D_expected_num[1:])
        and D_expected_num[0] < 0
        and D_expected_num[0] + D_expected_num[1] > 0,
        "D_tail_sweep_n8_to_n256": all(D_sweep),
        "degree2_prefix_positive": all(value > 0 for value in degree2_prefix),
        "degree3_prefix_positive": all(value > 0 for value in degree3_prefix),
        "degree3_F_tail_symbolic_identity": r_equal(F_derived, F_expected),
        "degree3_F_tail_negative_coefficient_certificate": all(
            value > 0 for value in F_positive_poly
        ),
        "jensen_degree2_degree3_sweep_shifts_0_to_256": all(jensen_sweep),
        "J4_coefficients_exact": [fraction_text(value) for value in j4_ascending]
        == ["1", "4", "39/8", "71487/32768", "11138032035/34359738368"],
        "J4_discriminant_exact_negative": j4_discriminant
        == Fraction(
            -691976798016209978769928197,
            158456325028528675187087900672,
        ),
        "entire_majorant_inputs": C_LIMIT < 1 and all(value < C_LIMIT for value in PREFIX.values()),
        "pinned_primary_fields": all(primary_checks.values()),
    }

    return {
        "classification": "independent_fraction_polynomial_audit_pass",
        "implementation": "python_stdlib_Fraction_no_SymPy_no_primary_code_import",
        "exact_symbolic_engine": "univariate polynomial numerator/denominator cross-products",
        "finite_sweep": {"n_min": 8, "n_max": MAX_N, "jensen_shift_min": 0, "jensen_shift_max": MAX_N},
        "reconstructed": {
            "q_2_to_q_8": [fraction_text(q(index)) for index in range(2, 9)],
            "D_2_to_D_7": [fraction_text(value) for value in prefix_D],
            "degree2_prefix_normalized_discriminants": [fraction_text(value) for value in degree2_prefix],
            "degree3_prefix_normalized_discriminants": [fraction_text(value) for value in degree3_prefix],
            "J4_coefficients_ascending": [fraction_text(value) for value in j4_ascending],
            "J4_discriminant": fraction_text(j4_discriminant),
            "PF2_minor_count_indices_0_to_20": len(toeplitz_pf2),
        },
        "primary_field_checks": primary_checks,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-report", action="store_true")
    args = parser.parse_args()
    audit = build_audit()
    if args.check_report:
        pinned = json.loads(AUDIT_PATH.read_text(encoding="utf-8"))
        if pinned != audit:
            raise SystemExit(f"audit mismatch: {AUDIT_PATH}")
    print(json.dumps(audit, indent=2, sort_keys=True))
    return 0 if audit["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
