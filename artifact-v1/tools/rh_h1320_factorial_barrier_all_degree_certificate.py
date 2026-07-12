#!/usr/bin/env python3
"""Canonical exact H1320 certificate with non-tautological base reconstruction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
REPORT_PATH = (
    ROOT / "research" / "riemann" / "h1320_factorial_barrier_all_degree_certificate.json"
)


def factored(value: sp.Expr) -> str:
    return str(sp.factor(value))


def q_exact(index: int) -> sp.Rational:
    if index == 2:
        return sp.Rational(1, 2)
    if index >= 3:
        return sp.Rational(2 * index - 3, 2 * index)
    raise ValueError("q is defined only for indices >=2")


def build_certificate() -> dict[str, object]:
    n, x = sp.symbols("n x", integer=True, positive=True)

    qn = (2 * n - 3) / (2 * n)
    qnext = (2 * n - 1) / (2 * (n + 1))
    qnext2 = (2 * n + 1) / (2 * (n + 2))
    q_increment = sp.factor(qnext - qn)

    c_n = 2 * n * (2 * n - 1)
    c_prev = 2 * (n - 1) * (2 * n - 3)
    c_next = 2 * (n + 1) * (2 * n + 1)
    theta = sp.factor(c_prev * c_next / c_n**2)
    t_ratio = sp.factor(theta * qnext / qn)
    factorial_slack = sp.factor(t_ratio - theta)
    unit_slack = sp.factor(t_ratio - (1 - 1 / n**2))
    d3 = sp.factor(
        (1 - qnext) ** 2 - qnext**2 * (1 - qn) * (1 - qnext2)
    )

    # Reconstruct n=2 entirely from definitions, rather than comparing two
    # hard-coded copies of 5/12.
    q2, q3, q4 = (q_exact(index) for index in (2, 3, 4))
    r1 = sp.Integer(1)
    r2 = sp.factor(r1 * q2)
    r3 = sp.factor(r2 * q3)
    r4 = sp.factor(r3 * q4)
    a0 = sp.Integer(1)
    a1 = sp.Integer(1)
    a2 = sp.factor(a1 * r2)
    a3 = sp.factor(a2 * r3)
    a4 = sp.factor(a3 * r4)
    a = [a0, a1, a2, a3, a4]

    moments = [sp.factorial(2 * index) * a[index] for index in range(4)]
    s1 = sp.factor(moments[1] / moments[0])
    s2 = sp.factor(moments[2] / moments[1])
    s3 = sp.factor(moments[3] / moments[2])
    t2 = sp.factor(s2 / s1)
    t3 = sp.factor(s3 / s2)
    t3_over_t2 = sp.factor(t3 / t2)
    theta2 = sp.factor(
        sp.Rational((2 * 1 * 1) * (2 * 3 * 5), (2 * 2 * 3) ** 2)
    )
    d2 = sp.factor((1 - q3) ** 2 - q3**2 * (1 - q2) * (1 - q4))

    gamma = [sp.factorial(j) * a[j] for j in range(4)]
    jensen3 = sp.expand(
        sum(sp.binomial(3, j) * gamma[j] * x**j for j in range(4))
    )
    jensen3_discriminant = sp.factor(sp.discriminant(jensen3, x))

    rows = (0, 1, 2, 3)
    columns = (1, 2, 3, 4)

    def coefficient(index: int) -> sp.Rational:
        return sp.Integer(0) if index < 0 else a[index]

    toeplitz = sp.Matrix(
        [[coefficient(column - row) for column in columns] for row in rows]
    )
    toeplitz_determinant = sp.factor(toeplitz.det())
    q2_var, q3_var, q4_var = sp.symbols("q2 q3 q4", positive=True)
    symbolic_a = [
        sp.Integer(1),
        sp.Integer(1),
        q2_var,
        q2_var**2 * q3_var,
        q2_var**3 * q3_var**2 * q4_var,
    ]
    symbolic_toeplitz = sp.Matrix(
        [[sp.Integer(0) if column - row < 0 else symbolic_a[column - row]
          for column in columns] for row in rows]
    )
    p4_from_determinant = sp.factor(symbolic_toeplitz.det())
    p4_expected = 1 - 3*q2_var + q2_var**2*(1 + 2*q3_var) - q2_var**3*q3_var**2*q4_var
    p4_at_counterexample = sp.factor(
        p4_from_determinant.subs({q2_var: q2, q3_var: q3, q4_var: q4})
    )

    r_closed = sp.factor(
        2 * sp.binomial(2 * n - 2, n - 1) / (n * 4 ** (n - 1))
    )

    raw_checks = {
        "q_increment_symbolic": sp.simplify(
            q_increment - 3 / (2 * n * (n + 1))
        )
        == 0,
        "t_ratio_symbolic": sp.simplify(
            t_ratio - (n - 1) * (2 * n + 1) / (n * (2 * n - 1))
        )
        == 0,
        "factorial_slack_symbolic": sp.simplify(
            factorial_slack
            - 3 * (n - 1) * (2 * n + 1) / (n**2 * (2 * n - 1) ** 2)
        )
        == 0,
        "unit_slack_symbolic": sp.simplify(
            unit_slack - (n - 1) / (n**2 * (2 * n - 1))
        )
        == 0,
        "d3_symbolic": sp.simplify(
            d3 - 9 * (12 * n - 1) / (16 * n * (n + 1) ** 2 * (n + 2))
        )
        == 0,
        "n2_q_reconstructed": q2 == q3 == sp.Rational(1, 2),
        "n2_factorial_ratio_reconstructed": t3_over_t2 == theta2,
        "d2_reconstructed_positive": d2 == sp.Rational(13, 64) and d2 > 0,
        "jensen_degree3_nonhyperbolic": jensen3_discriminant < 0,
        "toeplitz_order4_negative": toeplitz_determinant < 0,
        "p4_symbolic_identity": sp.simplify(
            p4_from_determinant - p4_expected
        ) == 0,
        "p4_evaluation_matches_toeplitz": p4_at_counterexample == toeplitz_determinant,
    }
    checks = {name: bool(value) for name, value in raw_checks.items()}

    return {
        "classification": (
            "factorial_and_eventual_unit_barriers_do_not_imply_"
            "all_degree_jensen_hyperbolicity"
        ),
        "arithmetic": "exact_sympy_rationals",
        "sequence": {
            "a_0": "1",
            "a_1": "1",
            "q_2": "1/2",
            "q_n_for_n_ge_3": "(2*n - 3)/(2*n)",
            "R_n_for_n_ge_2": factored(r_closed),
            "entire_majorant": "a_n <= 2^(n-1)/n!",
        },
        "symbolic_identities_for_n_ge_3": {
            "q_next_minus_q": factored(q_increment),
            "theta_h803": factored(theta),
            "T_next_over_T": factored(t_ratio),
            "h803_factorial_slack": factored(factorial_slack),
            "h804_unit_slack": factored(unit_slack),
            "translated_D3": factored(d3),
        },
        "special_n_2_reconstructed": {
            "q_2_q_3_q_4": [factored(value) for value in (q2, q3, q4)],
            "R_1_to_4": [factored(value) for value in (r1, r2, r3, r4)],
            "a_0_to_4": [factored(value) for value in a],
            "M_0_to_3": [factored(value) for value in moments],
            "S_1_to_3": [factored(value) for value in (s1, s2, s3)],
            "T_2_T_3": [factored(value) for value in (t2, t3)],
            "T_3_over_T_2": factored(t3_over_t2),
            "theta_2": factored(theta2),
            "D_2": factored(d2),
        },
        "jensen_falsifier": {
            "convention": "gamma_j=j!*a_j; J=sum(binomial(3,j)*gamma_j*x^j)",
            "gamma_0_to_3": [factored(value) for value in gamma],
            "J_gamma_3_0": str(jensen3),
            "discriminant": factored(jensen3_discriminant),
        },
        "toeplitz_falsifier": {
            "rows": list(rows),
            "columns": list(columns),
            "matrix": [
                [factored(value) for value in row] for row in toeplitz.tolist()
            ],
            "determinant": factored(toeplitz_determinant),
        },
        "order4_missing_invariant": {
            "normalization": "a_0=a_1=1",
            "normalized_formula": "1 - 3*q2 + q2**2*(1 + 2*q3) - q2**3*q3**2*q4",
            "derived_determinant": str(sp.expand(p4_from_determinant)),
            "counterexample_q2_q3_q4": [factored(value) for value in (q2, q3, q4)],
            "counterexample_value": factored(p4_at_counterexample),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "scope": (
            "Exact abstract counterexample to a generic implication from the "
            "listed ratio/barrier conditions; not an Xi or RH counterexample."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-report", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate()
    if args.check_report:
        pinned = json.loads(REPORT_PATH.read_text(encoding="utf-8"))
        if pinned != certificate:
            raise SystemExit(f"certificate mismatch: {REPORT_PATH}")
    print(json.dumps(certificate, indent=2, sort_keys=True))
    return 0 if certificate["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
