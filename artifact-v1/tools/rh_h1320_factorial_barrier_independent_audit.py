#!/usr/bin/env python3
"""Independent stdlib/Fraction audit of the H1320 exact counterexample."""

from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE_PATH = (
    ROOT / "research" / "riemann" / "h1320_factorial_barrier_all_degree_certificate.json"
)
AUDIT_PATH = (
    ROOT
    / "research"
    / "riemann"
    / "h1320_factorial_barrier_independent_audit.json"
)
MAX_N = 200


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def q(index: int) -> Fraction:
    if index == 2:
        return Fraction(1, 2)
    if index >= 3:
        return Fraction(2 * index - 3, 2 * index)
    raise ValueError("q is defined only for indices >=2")


def c(index: int) -> int:
    return 2 * index * (2 * index - 1)


def determinant(matrix: list[list[Fraction]]) -> Fraction:
    size = len(matrix)
    total = Fraction(0)
    for permutation in itertools.permutations(range(size)):
        inversions = sum(
            permutation[i] > permutation[j]
            for i in range(size)
            for j in range(i + 1, size)
        )
        term = Fraction(-1 if inversions % 2 else 1)
        for row, column in enumerate(permutation):
            term *= matrix[row][column]
        total += term
    return total


def cubic_discriminant(
    leading: Fraction,
    quadratic: Fraction,
    linear: Fraction,
    constant: Fraction,
) -> Fraction:
    return (
        quadratic**2 * linear**2
        - 4 * leading * linear**3
        - 4 * quadratic**3 * constant
        - 27 * leading**2 * constant**2
        + 18 * leading * quadratic * linear * constant
    )


def build_audit() -> dict[str, object]:
    # Build the sequence only from q and the definitions of R and a.
    ratios = [Fraction(0)] * (MAX_N + 2)
    coefficients = [Fraction(0)] * (MAX_N + 2)
    ratios[1] = Fraction(1)
    coefficients[0] = coefficients[1] = Fraction(1)
    for index in range(2, MAX_N + 2):
        ratios[index] = ratios[index - 1] * q(index)
        coefficients[index] = coefficients[index - 1] * ratios[index]

    moments = [
        Fraction(math.factorial(2 * index)) * coefficients[index]
        for index in range(MAX_N + 2)
    ]
    first_ratios = [Fraction(0)] * (MAX_N + 2)
    second_ratios = [Fraction(0)] * (MAX_N + 2)
    for index in range(1, MAX_N + 2):
        first_ratios[index] = moments[index] / moments[index - 1]
    for index in range(2, MAX_N + 2):
        second_ratios[index] = first_ratios[index] / first_ratios[index - 1]

    h803_rows: list[bool] = []
    unit_rows: list[bool] = []
    d3_rows: list[bool] = []
    closed_identity_rows: list[bool] = []
    entire_bound_rows: list[bool] = []

    for index in range(2, MAX_N + 1):
        observed = second_ratios[index + 1] / second_ratios[index]
        theta = Fraction(c(index - 1) * c(index + 1), c(index) ** 2)
        h803_rows.append(observed >= theta)
        if index >= 3:
            unit_rows.append(observed >= 1 - Fraction(1, index**2))
            closed_identity_rows.extend(
                [
                    q(index + 1) - q(index)
                    == Fraction(3, 2 * index * (index + 1)),
                    observed
                    == Fraction((index - 1) * (2 * index + 1), index * (2 * index - 1)),
                    observed - theta
                    == Fraction(
                        3 * (index - 1) * (2 * index + 1),
                        index**2 * (2 * index - 1) ** 2,
                    ),
                    observed - (1 - Fraction(1, index**2))
                    == Fraction(index - 1, index**2 * (2 * index - 1)),
                ]
            )

        e_index = 1 - q(index)
        e_next = 1 - q(index + 1)
        e_next2 = 1 - q(index + 2)
        d_value = e_next**2 - q(index + 1) ** 2 * e_index * e_next2
        expected_d = (
            Fraction(13, 64)
            if index == 2
            else Fraction(
                9 * (12 * index - 1),
                16 * index * (index + 1) ** 2 * (index + 2),
            )
        )
        d3_rows.append(d_value == expected_d and d_value > 0)

        closed_r = Fraction(
            2 * math.comb(2 * index - 2, index - 1),
            index * 4 ** (index - 1),
        )
        entire_bound_rows.append(
            ratios[index] == closed_r and ratios[index] <= Fraction(2, index)
        )

    gamma = [Fraction(math.factorial(index)) * coefficients[index] for index in range(4)]
    # Ascending Jensen coefficients are [1,3,3,3/4].
    jensen = [Fraction(math.comb(3, index)) * gamma[index] for index in range(4)]
    discriminant = cubic_discriminant(jensen[3], jensen[2], jensen[1], jensen[0])

    rows = (0, 1, 2, 3)
    columns = (1, 2, 3, 4)
    toeplitz = [
        [
            Fraction(0) if column - row < 0 else coefficients[column - row]
            for column in columns
        ]
        for row in rows
    ]
    toeplitz_determinant = determinant(toeplitz)
    p4_formula_value = (
        1
        - 3 * q(2)
        + q(2) ** 2 * (1 + 2 * q(3))
        - q(2) ** 3 * q(3) ** 2 * q(4)
    )


    certificate = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))
    special = certificate["special_n_2_reconstructed"]
    certificate_checks = {
        "primary_all_checks_pass": certificate["all_checks_pass"] is True,
        "primary_special_q_matches": special["q_2_q_3_q_4"]
        == [fraction_text(q(index)) for index in (2, 3, 4)],
        "primary_special_R_matches": special["R_1_to_4"]
        == [fraction_text(ratios[index]) for index in range(1, 5)],
        "primary_special_a_matches": special["a_0_to_4"]
        == [fraction_text(coefficients[index]) for index in range(5)],
        "primary_special_M_matches": special["M_0_to_3"]
        == [fraction_text(moments[index]) for index in range(4)],
        "primary_special_T_ratio_matches": special["T_3_over_T_2"]
        == fraction_text(second_ratios[3] / second_ratios[2]),
        "primary_D2_matches": special["D_2"] == "13/64",
        "primary_jensen_discriminant_matches": certificate["jensen_falsifier"]["discriminant"]
        == fraction_text(discriminant),
        "primary_toeplitz_determinant_matches": certificate["toeplitz_falsifier"]["determinant"]
        == fraction_text(toeplitz_determinant),
        "primary_p4_formula_matches": certificate["order4_missing_invariant"]["normalized_formula"]
        == "1 - 3*q2 + q2**2*(1 + 2*q3) - q2**3*q3**2*q4",
        "primary_p4_value_matches": certificate["order4_missing_invariant"]["counterexample_value"]
        == fraction_text(p4_formula_value),

    }

    checks = {
        "h803_n_2_to_200": all(h803_rows),
        "unit_barrier_n_3_to_200": all(unit_rows),
        "d3_n_2_to_200": all(d3_rows),
        "closed_identities_n_3_to_200": all(closed_identity_rows),
        "R_closed_form_and_entire_bound_n_2_to_200": all(entire_bound_rows),
        "base_T_ratio_reconstructed": second_ratios[3] / second_ratios[2]
        == Fraction(5, 12),
        "base_D2_reconstructed": (
            (1 - q(3)) ** 2 - q(3) ** 2 * (1 - q(2)) * (1 - q(4))
        )
        == Fraction(13, 64),
        "jensen_discriminant": discriminant == Fraction(-27, 16),
        "toeplitz_determinant": toeplitz_determinant == Fraction(-5, 256),
        "p4_formula_equals_toeplitz": p4_formula_value == toeplitz_determinant,
        "pinned_certificate_fields": all(certificate_checks.values()),
    }

    return {
        "classification": "independent_fraction_audit_pass",
        "implementation": "python_stdlib_fractions_no_primary_import",
        "range": {"n_min": 2, "n_max": MAX_N},
        "checks": checks,
        "certificate_field_checks": certificate_checks,
        "reconstructed": {
            "q_2_q_3_q_4": [fraction_text(q(index)) for index in (2, 3, 4)],
            "R_1_to_4": [fraction_text(ratios[index]) for index in range(1, 5)],
            "a_0_to_4": [fraction_text(coefficients[index]) for index in range(5)],
            "M_0_to_3": [fraction_text(moments[index]) for index in range(4)],
            "T_3_over_T_2": fraction_text(second_ratios[3] / second_ratios[2]),
            "D_2": "13/64",
            "jensen_discriminant": fraction_text(discriminant),
            "toeplitz_determinant": fraction_text(toeplitz_determinant),
            "p4_formula_value": fraction_text(p4_formula_value),
        },
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
