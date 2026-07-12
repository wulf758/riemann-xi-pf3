#!/usr/bin/env python3
"""Canonical SymPy certificate for the H1323 infinite Jensen 2--3 countermodel."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPORT = ROOT / "research" / "riemann" / "h1323_infinite_jensen_23_countermodel.json"


def exact_text(value: sp.Expr) -> str:
    return str(sp.factor(value))


def build_report() -> dict[str, object]:
    n = sp.symbols("n", integer=True, positive=True)
    x = sp.symbols("x")
    c = sp.Rational(201, 208)
    prefix = {
        2: sp.Rational(13, 32),
        3: sp.Rational(141, 256),
        4: sp.Rational(85, 128),
        5: sp.Rational(89, 128),
        6: sp.Rational(193, 256),
        7: sp.Rational(201, 256),
        8: sp.Rational(201, 256),
    }

    def q(index: int) -> sp.Expr:
        if index in prefix:
            return prefix[index]
        if index >= 8:
            return sp.factor(c * sp.Rational(2 * index - 3, 2 * index))
        raise ValueError(index)

    ratios = {1: sp.Integer(1)}
    coefficients = {0: sp.Integer(1), 1: sp.Integer(1)}
    for index in range(2, 9):
        ratios[index] = sp.factor(ratios[index - 1] * q(index))
        coefficients[index] = sp.factor(coefficients[index - 1] * ratios[index])
    gamma = {
        index: sp.factor(sp.factorial(index) * coefficients[index]) for index in range(9)
    }

    qn = sp.factor(c * (2 * n - 3) / (2 * n))
    qn1 = sp.factor(c * (2 * n - 1) / (2 * (n + 1)))
    q_delta = sp.factor(qn1 - qn)
    prefix_differences = [sp.factor(q(k + 1) - q(k)) for k in range(2, 8)]

    C = lambda z: 2 * z * (2 * z - 1)
    tau_ratio = sp.factor((qn1 / qn) * C(n + 1) * C(n - 1) / C(n) ** 2)
    unit_margin = sp.factor(tau_ratio - (1 - 1 / n**2))

    prefix_D: list[dict[str, object]] = []
    for index in range(2, 8):
        value = sp.factor(
            (1 - q(index + 1)) ** 2
            - q(index + 1) ** 2 * (1 - q(index)) * (1 - q(index + 2))
        )
        prefix_D.append({"n": index, "D_n": exact_text(value), "positive": bool(value > 0)})

    def q_tail(shift: int) -> sp.Expr:
        z = n + shift
        return sp.factor(c * (2 * z - 3) / (2 * z))

    D_tail = sp.factor(
        (1 - q_tail(1)) ** 2
        - q_tail(1) ** 2 * (1 - q_tail(0)) * (1 - q_tail(2))
    )
    D_num, D_den = map(sp.factor, sp.fraction(D_tail))
    D_coefficients = sp.Poly(D_num, n).all_coeffs()

    def Q(index: int) -> sp.Expr:
        return sp.factor(sp.Rational(index, index - 1) * q(index))

    degree2_prefix: list[dict[str, object]] = []
    degree3_prefix: list[dict[str, object]] = []
    for shift in range(6):
        u0 = Q(shift + 2)
        v0 = Q(shift + 3)
        disc2 = sp.factor(4 * (1 - u0))
        F0 = sp.factor(u0**2 * v0**2 - 6 * u0 * v0 + 4 * u0 + 4 * v0 - 3)
        disc3 = sp.factor(-27 * u0**2 * F0)
        degree2_prefix.append(
            {"shift": shift, "u": exact_text(u0), "normalized_discriminant": exact_text(disc2)}
        )
        degree3_prefix.append(
            {
                "shift": shift,
                "u": exact_text(u0),
                "v": exact_text(v0),
                "normalized_discriminant": exact_text(disc3),
            }
        )

    u = sp.factor(c * (2 * n + 1) / (2 * (n + 1)))
    v = sp.factor(c * (2 * n + 3) / (2 * (n + 2)))
    F = sp.factor(u**2 * v**2 - 6 * u * v + 4 * u + 4 * v - 3)
    F_num, F_den = map(sp.factor, sp.fraction(F))
    disc2_tail = sp.factor(4 * (1 - u))
    disc3_tail = sp.factor(-27 * u**2 * F)

    J4 = sp.expand(
        sum(sp.binomial(4, index) * gamma[index] * x**index for index in range(5))
    )
    J4_coefficients = [sp.factor(J4.coeff(x, index)) for index in range(5)]
    J4_discriminant = sp.factor(sp.discriminant(J4, x))

    checks = {
        "tail_connects_at_q8": bool(sp.simplify(qn.subs(n, 8) - prefix[8]) == 0),
        "prefix_q_nondecreasing": all(value >= 0 for value in prefix_differences),
        "prefix_q_below_c": all(value < c for value in prefix.values()),
        "tail_q_delta_identity": bool(
            sp.simplify(q_delta - 3 * c / (2 * n * (n + 1))) == 0
        ),
        "tail_q_below_c_below_one": bool(c < 1),
        "h804_tau_ratio_identity": bool(
            sp.simplify(tau_ratio - (n - 1) * (2 * n + 1) / (n * (2 * n - 1))) == 0
        ),
        "h804_unit_margin_identity": bool(
            sp.simplify(unit_margin - (n - 1) / (n**2 * (2 * n - 1))) == 0
        ),
        "D_prefix_positive_n2_to_n7": all(row["positive"] for row in prefix_D),
        "D_tail_factorization_exact": exact_text(D_tail)
        == "(2244592*n**4 + 297357088*n**3 + 13154707832*n**2 + 192552254264*n - 15372297693)/(29948379136*n*(n + 1)**2*(n + 2))",
        "D_tail_positive_for_n_ge_8": bool(
            all(coefficient > 0 for coefficient in D_coefficients[:-1])
            and D_coefficients[-1] < 0
            and D_coefficients[-2] + D_coefficients[-1] > 0
            and D_den.subs(n, 8) > 0
        ),
        "degree2_prefix_positive": all(
            sp.Rational(row["normalized_discriminant"]) > 0 for row in degree2_prefix
        ),
        "degree3_prefix_positive": all(
            sp.Rational(row["normalized_discriminant"]) > 0 for row in degree3_prefix
        ),
        "degree2_tail_positive": bool(c < 1),
        "degree3_tail_F_factorization_exact": exact_text(F)
        == "-3*(1509200*n**4 + 73874752*n**3 + 1236652088*n**2 + 7568113136*n + 5859746333)/(29948379136*(n + 1)**2*(n + 2)**2)",
        "degree3_tail_F_negative": bool(
            all(coefficient > 0 for coefficient in sp.Poly(-F_num, n).all_coeffs())
            and F_den.subs(n, 6) > 0
        ),
        "jensen_degree4_discriminant_negative": bool(J4_discriminant < 0),
        "entire_root_bound_tends_to_zero": bool(
            sp.limit(c ** ((n - 1) / 2), n, sp.oo) == 0
        ),
    }

    return {
        "classification": "h1323_infinite_exact_countermodel_to_jensen_degree_2_3_patch",
        "construction": {
            "c": exact_text(c),
            "q_2_to_q_8": [exact_text(q(index)) for index in range(2, 9)],
            "q_n_for_n_ge_8": "(201/208)*(2*n-3)/(2*n)",
            "definitions": [
                "R_1=a_0=a_1=1",
                "R_n/R_(n-1)=q_n",
                "a_n/a_(n-1)=R_n",
                "gamma_n=n!*a_n",
            ],
            "a_0_to_a_8": [exact_text(coefficients[index]) for index in range(9)],
            "gamma_0_to_gamma_8": [exact_text(gamma[index]) for index in range(9)],
        },
        "global_q_pf2_h803": {
            "prefix_differences": [exact_text(value) for value in prefix_differences],
            "tail_difference": exact_text(q_delta),
            "tail_bound": "0<q_n<c=201/208<1",
            "PF2_reason": "R_(n+1)/R_n=q_(n+1)<1; hence a is positive log-concave, which is equivalent to one-sided PF2",
            "H803_reason": "q_(n+1)>=q_n for every integer n>=2",
        },
        "h804_unit_barrier_tail": {
            "range": "n>=8",
            "definitions": "C_n=(2*n)*(2*n-1), tau_n=q_n*C_n/C_(n-1)",
            "tau_(n+1)/tau_n": exact_text(tau_ratio),
            "barrier": "1-1/n**2",
            "strict_margin": exact_text(unit_margin),
        },
        "entire_function": {
            "Phi": "sum_(n>=0) a_n*z**n",
            "majorant": [
                "q_n<c",
                "R_n<=c**(n-1)",
                "a_n<=c**(n*(n-1)/2)",
                "a_n**(1/n)<=c**((n-1)/2)->0",
            ],
            "conclusion": "Phi is entire by the root test",
        },
        "h851_translated_D3": {
            "definition": "D_n=(1-q_(n+1))**2-q_(n+1)**2*(1-q_n)*(1-q_(n+2))",
            "prefix_n_2_to_7": prefix_D,
            "tail_range": "n>=8",
            "tail_formula": exact_text(D_tail),
            "tail_numerator": exact_text(D_num),
            "tail_positivity_reason": "all nonconstant coefficients are positive and already the linear coefficient exceeds the absolute constant for n>=1",
        },
        "jensen_degree_2_3_all_shifts": {
            "convention": "J_gamma^(d,n)(X)=sum_(j=0)^d binom(d,j)*gamma_(n+j)*X**j",
            "normalization": "with y=P_(n+1)X and Q_m=m*q_m/(m-1): d=2 gives 1+2y+u*y**2; d=3 gives 1+3y+3u*y**2+u**2*v*y**3",
            "prefix_shifts": "n=0,...,5",
            "degree2_prefix": degree2_prefix,
            "degree3_prefix": degree3_prefix,
            "tail_shifts": "n>=6",
            "tail_u": exact_text(u),
            "tail_v": exact_text(v),
            "degree2_tail_normalized_discriminant": exact_text(disc2_tail),
            "F_definition": "u**2*v**2-6*u*v+4*u+4*v-3",
            "F_tail_factorization": exact_text(F),
            "degree3_tail_normalized_discriminant": exact_text(disc3_tail),
            "conclusion": "every degree-2 and degree-3 Jensen polynomial is strictly hyperbolic",
        },
        "degree4_falsifier": {
            "shift": 0,
            "coefficients_ascending": [exact_text(value) for value in J4_coefficients],
            "polynomial": str(J4),
            "discriminant": exact_text(J4_discriminant),
            "conclusion": "the negative discriminant forces one conjugate nonreal pair",
        },
        "scope": {
            "invalidated_patch": "PF2 + H803 + H804 tail + global H851 D_n + entire positivity + all shifted Jensen degrees 2 and 3 does not imply Jensen degree 4",
            "not_claimed": "the constructed sequence is not Xi and does not refute any Xi-specific property",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-report", action="store_true")
    args = parser.parse_args()
    report = build_report()
    if args.check_report:
        pinned = json.loads(DEFAULT_REPORT.read_text(encoding="utf-8"))
        if pinned != report:
            raise SystemExit(f"report mismatch: {DEFAULT_REPORT}")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
