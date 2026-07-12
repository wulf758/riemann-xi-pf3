#!/usr/bin/env python3
"""Independent addendum: strict PF2 partner with arbitrary gaps."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "research" / "riemann" / "h1326_pf2_arbitrary_gap_independent_addendum.json"


def toeplitz(sequence: list[Fraction], row: int, col: int) -> Fraction:
    index = col - row
    return sequence[index] if 0 <= index < len(sequence) else Fraction(0)


def minor2(
    sequence: list[Fraction], rows: tuple[int, int], cols: tuple[int, int]
) -> Fraction:
    i, j = rows
    p, q = cols
    return toeplitz(sequence, i, p) * toeplitz(sequence, j, q) - toeplitz(
        sequence, i, q
    ) * toeplitz(sequence, j, p)


def run(max_gap: int) -> dict:
    width = 3 * max_gap + 3
    ratios = [sp.Integer(1), *sp.symbols(f"R1:{width}")]
    coefficients = [sp.Integer(1)]
    for index in range(1, width):
        coefficients.append(sp.expand(coefficients[-1] * ratios[index]))

    exact_ratios = [Fraction(1)] + [Fraction(1, index + 1) for index in range(1, width)]
    exact_coefficients = [Fraction(1)]
    for index in range(1, width):
        exact_coefficients.append(exact_coefficients[-1] * exact_ratios[index])

    failures: list[dict] = []
    symbolic_count = 0
    shifted_product_count = 0
    exact_strict_count = 0
    for x in range(max_gap + 1):
        for alpha in range(1, max_gap + 1):
            for beta in range(1, max_gap + 1):
                determinant = (
                    coefficients[x + alpha] * coefficients[x + beta]
                    - coefficients[x + alpha + beta] * coefficients[x]
                )
                early = sp.prod(ratios[x + offset] for offset in range(1, alpha + 1))
                late = sp.prod(
                    ratios[x + beta + offset] for offset in range(1, alpha + 1)
                )
                factorization = coefficients[x] * coefficients[x + beta] * (early - late)
                if sp.expand(determinant - factorization) == 0:
                    symbolic_count += 1
                else:
                    failures.append(
                        {"kind": "factorization", "x": x, "alpha": alpha, "beta": beta}
                    )

                if all(x + offset < x + beta + offset for offset in range(1, alpha + 1)):
                    shifted_product_count += 1
                else:
                    failures.append(
                        {"kind": "shift_order", "x": x, "alpha": alpha, "beta": beta}
                    )

                exact_determinant = (
                    exact_coefficients[x + alpha] * exact_coefficients[x + beta]
                    - exact_coefficients[x + alpha + beta] * exact_coefficients[x]
                )
                if exact_determinant > 0:
                    exact_strict_count += 1
                else:
                    failures.append(
                        {"kind": "strict_canary", "x": x, "alpha": alpha, "beta": beta}
                    )

    triangular_count = 0
    for i in range(max_gap + 1):
        for j in range(i + 1, i + max_gap + 1):
            for p in range(i, j):
                for q in range(j, j + max_gap):
                    determinant = minor2(exact_coefficients, (i, j), (p, q))
                    diagonal_product = exact_coefficients[p - i] * exact_coefficients[q - j]
                    if determinant == diagonal_product and determinant > 0:
                        triangular_count += 1
                    else:
                        failures.append(
                            {"kind": "triangular", "rows": [i, j], "cols": [p, q]}
                        )

    product_cases = (max_gap + 1) * max_gap * max_gap
    all_checks = (
        not failures
        and symbolic_count == product_cases
        and shifted_product_count == product_cases
        and exact_strict_count == product_cases
        and triangular_count > 0
    )
    return {
        "schema": "rh_h1326_pf2_arbitrary_gap_independent_addendum.v0",
        "classification": (
            "h1326_arbitrary_gap_strict_pf2_partner_pass"
            if all_checks
            else "h1326_arbitrary_gap_strict_pf2_partner_fail"
        ),
        "all_checks_pass": all_checks,
        "identity": "D=a_(x+alpha)*a_(x+beta)-a_(x+alpha+beta)*a_x",
        "ratio_comparison": (
            "D/(a_x*a_(x+beta))=prod_{k=1..alpha}R_(x+k)"
            "-prod_{k=1..alpha}R_(x+beta+k)>0"
        ),
        "max_gap": max_gap,
        "product_cases": product_cases,
        "symbolic_factorizations_passed": symbolic_count,
        "shifted_ratio_product_orderings_passed": shifted_product_count,
        "exact_strict_rational_canaries_passed": exact_strict_count,
        "triangular_i_le_p_lt_j_cases_passed": triangular_count,
        "failed_checks": failures[:16],
    }


def write_markdown(report: dict, path: Path) -> None:
    lines = [
        "# H1326 Arbitrary-Gap PF2 Partner — Independent Addendum",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "For arbitrary positive gaps `alpha` and `beta`, the partner is",
        "",
        "`D=a_(x+alpha)a_(x+beta)-a_(x+alpha+beta)a_x`.",
        "",
        "After division by the positive factor `a_x a_(x+beta)`, this compares two "
        "products of `alpha` adjacent ratios, the second shifted `beta` places right. "
        "Strict decrease of the ratios makes every factor in the first product larger.",
        "",
        "The boundary `i<=p<j<=q` is separate: the lower-left Toeplitz entry is zero, "
        "so the minor is the strictly positive diagonal product.",
        "",
        f"Symbolic arbitrary-gap cases: {report['symbolic_factorizations_passed']}/"
        f"{report['product_cases']}; triangular cases: "
        f"{report['triangular_i_le_p_lt_j_cases_passed']}; failures: "
        f"{len(report['failed_checks'])}.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-gap", type=int, default=5)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()
    report = run(args.max_gap)
    print(json.dumps(report, indent=2))
    if not args.no_write:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(report, indent=2), encoding="utf-8")
        write_markdown(report, args.out.with_suffix(".md"))
    return 0 if report["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
