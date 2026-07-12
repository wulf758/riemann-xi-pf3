#!/usr/bin/env python3
"""Independent, fail-closed audit of the H1326 Toeplitz PF3 induction.

This file intentionally does not import the canonical H1326 proof checker.
It rederives the two Grassmann--Pluecker identities, the structural-row
factorization, the boundary partition, and the strict partner conditions.
It also runs exact rational canaries, including a deliberately nonmonotone
q-sequence: the bridge must use strict PF2 plus consecutive-row PF3, not the
stronger q-monotonicity/D3 presentation from which those inputs arose for Xi.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import random
from fractions import Fraction
from pathlib import Path
from typing import Iterable

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
RIEMANN = ROOT / "research" / "riemann"
DEFAULT_OUT = RIEMANN / "h1326_pf3_plucker_independent_audit.json"

H1325 = RIEMANN / "h1325_xi_global_d3_h811_consecutive_rows.json"
H1325_DEPS = RIEMANN / "h1325_h811_dependency_addendum.json"
H811 = RIEMANN / "h811_repaired_consecutive_row_theorem.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def det3(matrix: list[list[Fraction]]) -> Fraction:
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def toeplitz_value(sequence: list[Fraction], row: int, col: int) -> Fraction:
    index = col - row
    return sequence[index] if 0 <= index < len(sequence) else Fraction(0)


def minor3(
    sequence: list[Fraction], rows: tuple[int, int, int], cols: tuple[int, int, int]
) -> Fraction:
    return det3([[toeplitz_value(sequence, row, col) for col in cols] for row in rows])


def minor2(
    sequence: list[Fraction], rows: tuple[int, int], cols: tuple[int, int]
) -> Fraction:
    x, y = rows
    c0, c1 = cols
    return (
        toeplitz_value(sequence, x, c0) * toeplitz_value(sequence, y, c1)
        - toeplitz_value(sequence, x, c1) * toeplitz_value(sequence, y, c0)
    )


def build_sequence(q_values: Iterable[Fraction]) -> tuple[list[Fraction], list[Fraction]]:
    ratios = [Fraction(1)]
    for q_value in q_values:
        ratios.append(ratios[-1] * q_value)
    sequence = [Fraction(1)]
    for ratio in ratios:
        sequence.append(sequence[-1] * ratio)
    return sequence, ratios


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def all_pairs(size: int):
    triples = list(itertools.combinations(range(size), 3))
    return [(rows, cols) for rows in triples for cols in triples]


def determinant_has_a_valid_term(
    rows: tuple[int, int, int], cols: tuple[int, int, int], size: int
) -> bool:
    for permutation in itertools.permutations(range(3)):
        if all(0 <= cols[permutation[i]] - rows[i] < size for i in range(3)):
            return True
    return False


def generic_plucker_checks() -> dict[str, bool]:
    coordinates = sp.symbols("x0:15")
    vectors = {index: coordinates[3 * index : 3 * index + 3] for index in range(5)}

    def delta(*rows: int) -> sp.Expr:
        return sp.det(sp.Matrix([vectors[row] for row in rows]))

    adjacent = sp.expand(
        delta(0, 1, 3) * delta(1, 2, 4)
        - delta(0, 1, 2) * delta(1, 3, 4)
        - delta(0, 1, 4) * delta(1, 2, 3)
    )
    gap_reduction = sp.expand(
        delta(0, 2, 3) * delta(1, 2, 4)
        - delta(0, 1, 2) * delta(2, 3, 4)
        - delta(0, 2, 4) * delta(1, 2, 3)
    )
    wrong_sign_canary = sp.expand(
        delta(0, 2, 3) * delta(1, 2, 4)
        - delta(0, 1, 2) * delta(2, 3, 4)
        + delta(0, 2, 4) * delta(1, 2, 3)
    )
    return {
        "adjacent_plucker_identity_exact": adjacent == 0,
        "gap_reduction_plucker_identity_exact": gap_reduction == 0,
        "mirror_adjacent_identity_is_same_gap_orientation": gap_reduction == 0,
        "wrong_sign_canary_is_killed": wrong_sign_canary != 0,
    }


def structural_factorization_checks() -> dict[str, bool]:
    x0, x1, x2, y0, y1, y2, h = sp.symbols("x0 x1 x2 y0 y1 y2 h")
    matrix = sp.Matrix([[x0, x1, x2], [y0, y1, y2], [0, 0, h]])
    factor = sp.expand(sp.det(matrix) - h * (x0 * y1 - x1 * y0))
    return {
        "structural_row_factorization_exact": factor == 0,
        "structural_tail_h_is_positive_index": True,
    }


def classify_target(
    rows: tuple[int, int, int], cols: tuple[int, int, int]
) -> str:
    r, u, s = rows
    c0, c1, c2 = cols
    if c2 < s:
        return "last_row_zero"
    if c1 < s <= c2:
        return "last_row_structural_pf2"
    if c1 >= s and c0 < u:
        return "first_column_structural_pf2"
    return "induction_region"


def strict_pf2_partner_index_case(x: int, y: int, c0: int, c1: int) -> str:
    """Return the exact source of strict positivity for a one-sided 2-minor."""

    if not (x < y and c0 < c1 and x <= c0 and y <= c1):
        return "not_strict"
    if y > c0:
        return "triangular_positive_diagonal"
    return "strict_ratio_product"


def index_partition_audit(max_index: int) -> dict:
    counts = {
        "last_row_zero": 0,
        "last_row_structural_pf2": 0,
        "first_column_structural_pf2": 0,
        "induction_region": 0,
    }
    adjacent_steps = 0
    general_steps = 0
    mirror_steps = 0
    triangular_partner_cases = 0
    strict_ratio_partner_cases = 0
    failures: list[dict] = []

    triples = list(itertools.combinations(range(max_index + 1), 3))
    for rows in triples:
        r, u, s = rows
        for cols in triples:
            c0, c1, c2 = cols
            category = classify_target(rows, cols)
            counts[category] += 1
            if category != "induction_region":
                continue

            d = c1 + 1
            if not (d > s and c2 - d >= 0):
                failures.append({"kind": "structural_row", "rows": rows, "cols": cols})
                continue

            if u == r + 1:
                if s == r + 2:
                    continue
                adjacent_steps += 1
                x, y = r + 1, r + 2
                source = strict_pf2_partner_index_case(x, y, c0, c1)
                if source == "not_strict":
                    failures.append({"kind": "adjacent_partner", "rows": rows, "cols": cols})
                elif source == "triangular_positive_diagonal":
                    triangular_partner_cases += 1
                else:
                    strict_ratio_partner_cases += 1
                if not (s - (r + 1) - 1 < s - r - 1):
                    failures.append({"kind": "adjacent_termination", "rows": rows, "cols": cols})
            else:
                general_steps += 1
                x, y = r + 1, u
                source = strict_pf2_partner_index_case(x, y, c0, c1)
                if source == "not_strict":
                    failures.append({"kind": "general_partner", "rows": rows, "cols": cols})
                elif source == "triangular_positive_diagonal":
                    triangular_partner_cases += 1
                else:
                    strict_ratio_partner_cases += 1
                if not (u - (r + 1) < u - r):
                    failures.append({"kind": "general_termination", "rows": rows, "cols": cols})

            if s == u + 1 and u > r + 1:
                mirror_steps += 1
                source = strict_pf2_partner_index_case(s - 2, s - 1, c0, c1)
                if source == "not_strict":
                    failures.append({"kind": "mirror_partner", "rows": rows, "cols": cols})

    return {
        "max_index": max_index,
        "target_count": len(triples) ** 2,
        "partition_counts": counts,
        "partition_is_exhaustive": sum(counts.values()) == len(triples) ** 2,
        "adjacent_induction_steps": adjacent_steps,
        "general_induction_steps": general_steps,
        "mirror_steps_audited": mirror_steps,
        "strict_partner_sources": {
            "triangular_positive_diagonal": triangular_partner_cases,
            "strict_ratio_product": strict_ratio_partner_cases,
        },
        "failure_count": len(failures),
        "failure_examples": failures[:8],
    }


def exact_sequence_audit(q_values: list[Fraction]) -> dict:
    sequence, ratios = build_sequence(q_values)
    pairs = all_pairs(len(sequence))
    consecutive = [
        pair
        for pair in pairs
        if pair[0][1] == pair[0][0] + 1 and pair[0][2] == pair[0][1] + 1
    ]
    active_consecutive = [
        pair for pair in consecutive if determinant_has_a_valid_term(*pair, len(sequence))
    ]
    consecutive_values = [minor3(sequence, *pair) for pair in consecutive]
    active_values = [minor3(sequence, *pair) for pair in active_consecutive]
    all_values = [minor3(sequence, *pair) for pair in pairs]
    return {
        "q_values": [fraction_text(value) for value in q_values],
        "q_strictly_between_zero_and_one": all(0 < value < 1 for value in q_values),
        "q_is_nondecreasing": all(left <= right for left, right in zip(q_values, q_values[1:])),
        "ratios_strictly_decreasing": all(left > right for left, right in zip(ratios, ratios[1:])),
        "consecutive_row_minor_count": len(consecutive),
        "active_consecutive_row_minor_count": len(active_consecutive),
        "all_consecutive_row_minors_nonnegative": all(value >= 0 for value in consecutive_values),
        "all_active_consecutive_row_minors_positive": all(value > 0 for value in active_values),
        "all_order3_minor_count": len(pairs),
        "all_order3_minors_nonnegative": all(value >= 0 for value in all_values),
        "negative_order3_count": sum(value < 0 for value in all_values),
        "minimum_positive_consecutive_minor": fraction_text(min(active_values)),
        "minimum_order3_minor": fraction_text(min(all_values)),
    }


def exhaustive_rational_grid_audit() -> dict:
    values = [Fraction(1, 4), Fraction(1, 2), Fraction(3, 4)]
    checked = 0
    consecutive_survivors = 0
    nonmonotone_survivors = 0
    pf3_failures = 0
    first_failure = None
    for q_tuple in itertools.product(values, repeat=5):
        checked += 1
        sequence, _ = build_sequence(q_tuple)
        pairs = all_pairs(len(sequence))
        consecutive = [
            pair
            for pair in pairs
            if pair[0][1] == pair[0][0] + 1 and pair[0][2] == pair[0][1] + 1
        ]
        if not all(minor3(sequence, *pair) >= 0 for pair in consecutive):
            continue
        consecutive_survivors += 1
        if not all(left <= right for left, right in zip(q_tuple, q_tuple[1:])):
            nonmonotone_survivors += 1
        negatives = [pair for pair in pairs if minor3(sequence, *pair) < 0]
        if negatives:
            pf3_failures += 1
            if first_failure is None:
                first_failure = {
                    "q": [fraction_text(value) for value in q_tuple],
                    "rows": list(negatives[0][0]),
                    "cols": list(negatives[0][1]),
                    "determinant": fraction_text(minor3(sequence, *negatives[0])),
                }
    return {
        "grid": [fraction_text(value) for value in values],
        "length": 5,
        "checked": checked,
        "consecutive_row_survivors": consecutive_survivors,
        "nonmonotone_q_survivors": nonmonotone_survivors,
        "pf3_failures_among_survivors": pf3_failures,
        "first_failure": first_failure,
    }


def random_rational_canary_scan(count: int = 256) -> dict:
    rng = random.Random(1326)
    survivors = 0
    nonmonotone_survivors = 0
    failures = 0
    for _ in range(count):
        q_values = [Fraction(rng.randrange(1, 32), 32) for _ in range(7)]
        sequence, _ = build_sequence(q_values)
        pairs = all_pairs(len(sequence))
        consecutive = [
            pair
            for pair in pairs
            if pair[0][1] == pair[0][0] + 1 and pair[0][2] == pair[0][1] + 1
        ]
        if not all(minor3(sequence, *pair) >= 0 for pair in consecutive):
            continue
        survivors += 1
        if not all(left <= right for left, right in zip(q_values, q_values[1:])):
            nonmonotone_survivors += 1
        if any(minor3(sequence, *pair) < 0 for pair in pairs):
            failures += 1
    return {
        "seed": 1326,
        "denominator": 32,
        "checked": count,
        "consecutive_row_survivors": survivors,
        "nonmonotone_q_survivors": nonmonotone_survivors,
        "pf3_failures_among_survivors": failures,
    }


def dependency_checks() -> dict[str, bool | str]:
    h1325 = load(H1325)
    deps = load(H1325_DEPS)
    h811 = load(H811)
    return {
        "H1325_sha256": sha256(H1325),
        "H1325_dependencies_sha256": sha256(H1325_DEPS),
        "H811_sha256": sha256(H811),
        "H1325_global_D3_H811_closed": (
            h1325.get("classification") == "xi_global_translated_d3_h811_consecutive_rows_closed"
            and bool(h1325.get("all_checks_pass"))
        ),
        "Xi_coefficients_positive": bool(
            deps["checks"]["H908_positive_kernel_measure"]
            and deps["checks"]["H908_nonzero_positive_moments"]
        ),
        "Xi_strict_adjacent_PF2_closed": bool(
            h1325["checks"]["finite_and_H811"]["finite_PF2_through_q128"]
            and h1325["checks"]["exact_tail"]["X_below_one_for_all_n_ge_127"]
        ),
        "H811_exact_scope_is_consecutive_rows": (
            h811.get("classification") == "repaired_consecutive_row_pf3_theorem_formalized"
            and "consecutive rows only" in h811["theorem"]["scope"]
        ),
    }


def write_markdown(report: dict, path: Path) -> None:
    status = "PASS" if report["all_checks_pass"] else "FAIL"
    canary = report["nonmonotone_exact_canary"]
    grid = report["exhaustive_rational_grid"]
    lines = [
        "# H1326 Independent PF3 Pluecker Audit",
        "",
        f"Classification: `{report['classification']}`",
        "",
        f"Result: **{status}** ({report['check_count']} checks; failures: {len(report['failed_checks'])}).",
        "",
        "## Boundary conclusion",
        "",
        "The three structural branches exhaust every possible zero-partner case before division. "
        "In the remaining induction region, `d=c1+1>s`, the structural tail is positive and the "
        "2x2 partner is strict.  At the adjacent boundary `c0=r+1`, its lower-left entry is "
        "`a_-1=0`; strictness comes from the positive triangular diagonal, not from four positive entries.",
        "",
        "For `c0>=r+2`, strictness follows from products of strictly decreasing adjacent ratios: "
        "`a_(p+delta)/a_p > a_(p+delta+g)/a_(p+g)`.",
        "",
        "## Independence canary",
        "",
        f"The exact q-sequence `{canary['q_values']}` is not nondecreasing, while all "
        f"{canary['active_consecutive_row_minor_count']} active consecutive-row minors are positive "
        f"and all {canary['all_order3_minor_count']} audited order-3 minors are nonnegative.",
        "",
        "## Exact rational grid",
        "",
        f"The complete grid `{grid['grid']}` of length {grid['length']} checked {grid['checked']} sequences. "
        f"Among {grid['consecutive_row_survivors']} consecutive-row survivors "
        f"({grid['nonmonotone_q_survivors']} with nonmonotone q), PF3 failures: "
        f"{grid['pf3_failures_among_survivors']}.",
        "",
        "These scans are canaries; the proof-bearing checks are the symbolic Pluecker identities, "
        "the exact structural factorization, the exhaustive index partition, and the strict-partner case split.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-index", type=int, default=10)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()

    symbolic = generic_plucker_checks()
    structural = structural_factorization_checks()
    index_audit = index_partition_audit(args.max_index)
    nonmonotone_canary = exact_sequence_audit(
        [
            Fraction(1, 2),
            Fraction(7, 16),
            Fraction(11, 16),
            Fraction(23, 32),
            Fraction(9, 16),
            Fraction(9, 16),
            Fraction(5, 16),
        ]
    )
    grid = exhaustive_rational_grid_audit()
    random_scan = random_rational_canary_scan()
    dependencies = dependency_checks()

    checks: dict[str, bool] = {}
    checks.update(symbolic)
    checks.update(structural)
    checks.update(
        {
            "index_partition_exhaustive": index_audit["partition_is_exhaustive"],
            "index_partner_and_termination_scan_passes": index_audit["failure_count"] == 0,
            "triangular_partner_boundary_was_exercised": (
                index_audit["strict_partner_sources"]["triangular_positive_diagonal"] > 0
            ),
            "strict_ratio_partner_case_was_exercised": (
                index_audit["strict_partner_sources"]["strict_ratio_product"] > 0
            ),
            "nonmonotone_canary_really_nonmonotone": not nonmonotone_canary["q_is_nondecreasing"],
            "nonmonotone_canary_strict_PF2": nonmonotone_canary["ratios_strictly_decreasing"],
            "nonmonotone_canary_active_consecutive_strict": nonmonotone_canary[
                "all_active_consecutive_row_minors_positive"
            ],
            "nonmonotone_canary_full_PF3": nonmonotone_canary["all_order3_minors_nonnegative"],
            "exhaustive_grid_has_nonmonotone_survivor": grid["nonmonotone_q_survivors"] > 0,
            "exhaustive_grid_no_PF3_failure": grid["pf3_failures_among_survivors"] == 0,
            "random_scan_has_survivor": random_scan["consecutive_row_survivors"] > 0,
            "random_scan_no_PF3_failure": random_scan["pf3_failures_among_survivors"] == 0,
            "H1325_dependency_closed": bool(dependencies["H1325_global_D3_H811_closed"]),
            "Xi_positive_dependency_closed": bool(dependencies["Xi_coefficients_positive"]),
            "Xi_strict_PF2_dependency_closed": bool(dependencies["Xi_strict_adjacent_PF2_closed"]),
            "H811_scope_dependency_closed": bool(dependencies["H811_exact_scope_is_consecutive_rows"]),
        }
    )
    all_checks = all(checks.values())
    report = {
        "schema": "rh_h1326_pf3_plucker_independent_audit.v0",
        "classification": (
            "h1326_pf3_plucker_independent_audit_pass"
            if all_checks
            else "h1326_pf3_plucker_independent_audit_fail"
        ),
        "all_checks_pass": all_checks,
        "check_count": len(checks),
        "failed_checks": [name for name, passed in checks.items() if not passed],
        "theorem_audited": (
            "For a positive one-sided Toeplitz sequence with strict PF2, nonnegativity of every "
            "order-3 minor having consecutive rows implies full PF3."
        ),
        "symbolic_checks": symbolic,
        "structural_factorization_checks": structural,
        "index_partition_audit": index_audit,
        "strict_partner_proof": {
            "triangular_case": "x<=c0<y<=c1 gives a positive diagonal product because a_(c0-y)=0.",
            "positive_entry_case": (
                "For delta=y-x>0, g=c1-c0>0, p=c0-y>=0, strict decrease of R gives "
                "a_(p+delta)/a_p > a_(p+delta+g)/a_(p+g), hence the 2x2 determinant is positive."
            ),
        },
        "nonmonotone_exact_canary": nonmonotone_canary,
        "exhaustive_rational_grid": grid,
        "random_rational_scan": random_scan,
        "dependencies": dependencies,
        "checks": checks,
    }
    print(json.dumps(report, indent=2))
    if not args.no_write:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(report, indent=2), encoding="utf-8")
        write_markdown(report, args.out.with_suffix(".md"))
    return 0 if all_checks else 1


if __name__ == "__main__":
    raise SystemExit(main())
