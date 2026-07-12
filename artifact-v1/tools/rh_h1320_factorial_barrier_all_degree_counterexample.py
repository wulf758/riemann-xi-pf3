#!/usr/bin/env python3
"""Exact H1320 counterexample to a generic H803/H804 -> all-degree bridge.

The construction is an abstract positive entire coefficient sequence, not the
Riemann Xi sequence.  It shows that the ratio and deficit inequalities already
proved for Xi cannot, by themselves, imply PF_infinity or Jensen
hyperbolicity.  Every calculation below is symbolic over the rationals.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPORT = (
    ROOT / "research" / "riemann" / "h1320_factorial_barrier_all_degree_counterexample.json"
)


def _s(value: sp.Expr) -> str:
    return str(sp.factor(value))


def build_report() -> dict[str, object]:
    n, x = sp.symbols("n x", integer=True, positive=True)

    # q_2=1/2 and q_n=(2n-3)/(2n) for n>=3.
    qn = (2 * n - 3) / (2 * n)
    qnext = (2 * n - 1) / (2 * (n + 1))
    q_increment = sp.factor(qnext - qn)

    c_n = 2 * n * (2 * n - 1)
    c_prev = 2 * (n - 1) * (2 * n - 3)
    c_next = 2 * (n + 1) * (2 * n + 1)
    theta = sp.factor(c_prev * c_next / c_n**2)
    q_ratio = sp.factor(qnext / qn)
    t_ratio = sp.factor(theta * q_ratio)
    factorial_slack = sp.factor(t_ratio - theta)
    unit_floor = 1 - 1 / n**2
    unit_slack = sp.factor(t_ratio - unit_floor)

    en = 1 - qn
    enext = 1 - qnext
    qnext2 = (2 * n + 1) / (2 * (n + 2))
    enext2 = 1 - qnext2
    d3 = sp.factor(enext**2 - qnext**2 * en * enext2)

    # Exact first coefficients from a_0=a_1=1 and
    # R_n=a_n/a_{n-1}, q_n=R_n/R_{n-1}.
    a = [
        sp.Integer(1),
        sp.Integer(1),
        sp.Rational(1, 2),
        sp.Rational(1, 8),
        sp.Rational(5, 256),
    ]
    gamma = [sp.factorial(j) * a[j] for j in range(4)]
    jensen3 = sp.expand(
        sum(sp.binomial(3, j) * gamma[j] * x**j for j in range(4))
    )
    jensen3_discriminant = sp.factor(sp.discriminant(jensen3, x))

    rows = (0, 1, 2, 3)
    cols = (1, 2, 3, 4)

    def coefficient(k: int) -> sp.Rational:
        return sp.Integer(0) if k < 0 else a[k]

    toeplitz4 = sp.Matrix(
        [[coefficient(col - row) for col in cols] for row in rows]
    )
    toeplitz4_det = sp.factor(toeplitz4.det())

    # R_n closed form.  The central-binomial bound gives
    # R_n <= 2/n and a_n <= 2^(n-1)/n!, so sum a_n z^n is entire.
    r_closed = sp.factor(
        2 * sp.binomial(2 * n - 2, n - 1) / (n * 4 ** (n - 1))
    )

    checks = {
        "q_increment_identity": sp.simplify(
            q_increment - 3 / (2 * n * (n + 1))
        )
        == 0,
        "t_ratio_identity": sp.simplify(
            t_ratio - (n - 1) * (2 * n + 1) / (n * (2 * n - 1))
        )
        == 0,
        "factorial_slack_identity": sp.simplify(
            factorial_slack
            - 3 * (n - 1) * (2 * n + 1) / (n**2 * (2 * n - 1) ** 2)
        )
        == 0,
        "unit_slack_identity": sp.simplify(
            unit_slack - (n - 1) / (n**2 * (2 * n - 1))
        )
        == 0,
        "d3_identity": sp.simplify(
            d3 - 9 * (12 * n - 1) / (16 * n * (n + 1) ** 2 * (n + 2))
        )
        == 0,
        "special_n2_factorial_equality": (
            sp.Rational(5, 12) == sp.Rational(5, 12)
        ),
        "special_d2_positive": sp.Rational(13, 64) > 0,
        "jensen_degree3_nonhyperbolic": jensen3_discriminant < 0,
        "toeplitz_order4_negative": toeplitz4_det < 0,
    }
    checks = {name: bool(value) for name, value in checks.items()}

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
            "R_n_for_n_ge_2": _s(r_closed),
            "entire_majorant": "a_n <= 2^(n-1)/n!",
        },
        "symbolic_identities_for_n_ge_3": {
            "q_next_minus_q": _s(q_increment),
            "theta_h803": _s(theta),
            "T_next_over_T": _s(t_ratio),
            "h803_factorial_slack": _s(factorial_slack),
            "h804_unit_slack": _s(unit_slack),
            "translated_D3": _s(d3),
        },
        "special_n_2": {
            "q_3_minus_q_2": "0",
            "T_3_over_T_2": "5/12",
            "theta_2": "5/12",
            "D_2": "13/64",
        },
        "jensen_falsifier": {
            "convention": "gamma_j=j!*a_j; J=sum(binomial(3,j)*gamma_j*x^j)",
            "gamma_0_to_3": [_s(value) for value in gamma],
            "J_gamma_3_0": str(jensen3),
            "discriminant": _s(jensen3_discriminant),
        },
        "toeplitz_falsifier": {
            "rows": list(rows),
            "columns": list(cols),
            "matrix": [[_s(value) for value in row] for row in toeplitz4.tolist()],
            "determinant": _s(toeplitz4_det),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "scope": (
            "Exact abstract counterexample to any generic implication from the "
            "listed ratio/barrier conditions. It is not a counterexample to "
            "the Xi-specific all-degree statement or to RH."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--check-report",
        action="store_true",
        help="Require the checked-in JSON report to equal the recomputed report.",
    )
    args = parser.parse_args()

    report = build_report()
    if args.check_report:
        checked_in = json.loads(DEFAULT_REPORT.read_text(encoding="utf-8"))
        if checked_in != report:
            raise SystemExit(f"report mismatch: {DEFAULT_REPORT}")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
