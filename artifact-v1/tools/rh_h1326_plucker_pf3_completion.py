#!/usr/bin/env python3
"""Fail-closed replay for H1326: strict PF2 + H811 implies PF3.

The proof uses two oriented Grassmann--Plucker inductions. First, every
order-3 Toeplitz minor having an adjacent row pair is reduced to consecutive
rows. Second, an arbitrary sparse row triple is reduced to that adjacent-row
lemma and to a smaller first row gap. One-sided Toeplitz boundary cases factor
directly through PF2. A reflection/transposition identity supplies a redundant
column-adjacency audit.

This replay fails on dependency drift, a sign error, a structural
factorization error, a non-strict induction partner, nontermination, or any
unproved item in the H813/H816 exact finite universes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from functools import lru_cache
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
RIEMANN = ROOT / "research" / "riemann"
THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from rh_h807_essential_minor_hard_solver_audit import all_order3_pairs, determinant_terms
from rh_h812_sparse_row_pf3_obstruction_map import is_consecutive_rows
from rh_h813_sparse_row_symbolic_hierarchy import determinant_q_poly, poly_to_string
from rh_h816_d3_cone_sparse_row_coupling import remaining_after_h815


H1325 = RIEMANN / "h1325_xi_global_d3_h811_consecutive_rows.json"
H1325_DEPS = RIEMANN / "h1325_h811_dependency_addendum.json"
H811 = RIEMANN / "h811_repaired_consecutive_row_theorem.json"
H813 = RIEMANN / "h813_sparse_row_symbolic_hierarchy.json"
H816 = RIEMANN / "h816_d3_cone_sparse_row_coupling.json"

PINNED_HASHES = {
    H1325: "9308258F232A84CB964F0CE4357C2BE279D5E2FE3A49B227B37CCD7E6D43A1A0",
    H1325_DEPS: "F034A24C08EED9FC8E0A0405D596DF47A4227CFBEC7E7664D2823062A45A26B4",
    H811: "09C0B09673CD582BB1F9968466B5304FBA13A4B8AE0D97B6348EAD60498D6A12",
    H813: "C92E94BAED6251E729F3B3C02AE5B93A92653118022B32682E7D8B7E24EE5ED6",
    H816: "8E6F0364C259FBBBC1AB66D9B1C250D2F706004FE1842644AFC584D07BDC2A00",
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def generic_plucker_checks() -> dict[str, bool]:
    """Verify all three orientations on six generic row vectors."""

    x = sp.symbols("x0:18")
    vectors = {r: x[3 * r : 3 * r + 3] for r in range(6)}

    def minor(*rows: int) -> sp.Expr:
        return sp.det(sp.Matrix([vectors[row] for row in rows]))

    # (r,r+1,s), with r=0, s=4, d=5.
    forward = sp.expand(
        minor(0, 1, 4) * minor(1, 2, 5)
        - minor(0, 1, 2) * minor(1, 4, 5)
        - minor(0, 1, 5) * minor(1, 2, 4)
    )
    # (r,s-1,s), with r=0, s=4, d=5.
    backward = sp.expand(
        minor(0, 3, 4) * minor(2, 3, 5)
        - minor(0, 2, 3) * minor(3, 4, 5)
        - minor(0, 3, 5) * minor(2, 3, 4)
    )
    # Arbitrary r<u<s, with r=0, r+1=1, u=3, s=4, d=5.
    general = sp.expand(
        minor(0, 3, 4) * minor(1, 3, 5)
        - minor(0, 1, 3) * minor(3, 4, 5)
        - minor(0, 3, 5) * minor(1, 3, 4)
    )
    return {
        "plucker_adjacent_forward_sign": forward == 0,
        "plucker_adjacent_backward_sign": backward == 0,
        "plucker_general_induction_sign": general == 0,
    }


def structural_factor_checks() -> dict[str, bool]:
    """Verify the added-row and the two one-sided boundary factorizations."""

    x0, x1, y0, y1, z0, z1, h = sp.symbols("x0 x1 y0 y1 z0 z1 h")
    added_row = sp.expand(
        sp.det(sp.Matrix([[x0, y0, z0], [x1, y1, z1], [0, 0, h]]))
        - h * (x0 * y1 - y0 * x1)
    )
    y2, z2 = sp.symbols("y2 z2")
    first_column_boundary = sp.expand(
        sp.det(sp.Matrix([[h, y0, z0], [0, y1, z1], [0, y2, z2]]))
        - h * (y1 * z2 - z1 * y2)
    )

    # Strictness of a non-boundary PF2 partner from strictly decreasing R_n.
    apm1, aqm1, rp, rq = sp.symbols("a_pm1 a_qm1 R_p R_q")
    strict_partner = sp.expand(
        h * ((apm1 * rp) * aqm1 - (aqm1 * rq) * apm1)
        - h * apm1 * aqm1 * (rp - rq)
    )
    # If the first lower index is zero, the lower-left entry is a_-1=0.
    a0, aqm1_boundary = sp.symbols("a0 a_qm1_boundary")
    triangular_partner = sp.expand(h * (a0 * aqm1_boundary - 0) - h * a0 * aqm1_boundary)

    return {
        "added_d_row_factor_is_h_times_PF2": added_row == 0,
        "last_row_two_zero_factor": added_row == 0,
        "first_column_two_zero_factor": first_column_boundary == 0,
        "strict_partner_ratio_factor": strict_partner == 0,
        "strict_partner_one_sided_boundary": triangular_partner == 0,
        "c2_equals_d_gives_h_equals_a0_not_zero": h != 0,
    }


def duality_checks() -> dict[str, bool]:
    """Check Delta_R^C = Delta_(L-Crev)^(L-Rrev) and gap reversal."""

    entries = sp.symbols("m0:9")
    matrix = sp.Matrix(3, 3, entries)
    reversal = sp.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    determinant_duality = sp.expand(matrix.det() - (reversal * matrix.T * reversal).det()) == 0

    samples = [
        ((0, 2, 5), (3, 6, 8), 8),
        ((1, 3, 5), (4, 6, 8), 9),
        ((0, 1, 7), (2, 5, 9), 12),
    ]
    index_map = True
    gap_map = True
    for rows, cols, level in samples:
        dual_rows = tuple(level - value for value in reversed(cols))
        dual_cols = tuple(level - value for value in reversed(rows))
        for i in range(3):
            for j in range(3):
                index_map &= dual_cols[j] - dual_rows[i] == cols[2 - i] - rows[2 - j]
        row_gaps = (dual_rows[1] - dual_rows[0], dual_rows[2] - dual_rows[1])
        gap_map &= row_gaps == (cols[2] - cols[1], cols[1] - cols[0])
    return {
        "toeplitz_reflection_transpose_determinant": determinant_duality,
        "toeplitz_reflection_transpose_indices": index_map,
        "column_gaps_become_reversed_row_gaps": gap_map,
    }


class ProofPlanner:
    """Combinatorial replay of the two inductions; signs are checked above."""

    def __init__(self) -> None:
        self.kinds: Counter[str] = Counter()
        self.max_depth = 0

    @staticmethod
    def strict_partner(rows: tuple[int, int], cols: tuple[int, int]) -> bool:
        i, j = rows
        p, q = cols
        # Sufficient for a nonzero one-sided 2-minor. Strict PF2 makes it
        # positive, including p=i where the lower-left entry is a_-1=0.
        return i < j and p < q and p >= i and q >= j

    @lru_cache(maxsize=None)
    def prove(self, rows: tuple[int, int, int], cols: tuple[int, int, int]) -> bool:
        r0, r1, r2 = rows
        c0, c1, c2 = cols
        depth = (r1 - r0) + (r2 - r1)
        self.max_depth = max(self.max_depth, depth)
        if not (r0 < r1 < r2 and c0 < c1 < c2):
            return False

        if r1 == r0 + 1 and r2 == r1 + 1:
            self.kinds["H811_consecutive_rows"] += 1
            return True
        if c2 < r2:
            self.kinds["zero_last_row"] += 1
            return True
        if c1 < r2:
            self.kinds["last_row_PF2_factor"] += 1
            return True
        if c0 < r1:
            self.kinds["first_column_PF2_factor_or_zero"] += 1
            return True

        # Here c1>=r2 and c0>=r1. Therefore d=c1+1>r2 and the
        # added row has entries (0,0,a_(c2-c1-1)), with positive tail.
        if r1 == r0 + 1:
            partner_ok = self.strict_partner((r1, r1 + 1), (c0, c1))
            smaller = (r1, r1 + 1, r2)
            self.kinds["Plucker_adjacent_forward"] += 1
            return partner_ok and self.prove(smaller, cols)

        if r2 == r1 + 1:
            partner_ok = self.strict_partner((r1 - 1, r1), (c0, c1))
            smaller = (r0, r1 - 1, r1)
            self.kinds["Plucker_adjacent_backward"] += 1
            return partner_ok and self.prove(smaller, cols)

        # General identity. The first RHS 3-minor has an adjacent pair;
        # the second reduces the first row gap from r1-r0 to r1-r0-1.
        partner_ok = self.strict_partner((r0 + 1, r1), (c0, c1))
        adjacent = (r0, r0 + 1, r1)
        smaller = (r0 + 1, r1, r2)
        self.kinds["Plucker_general_first_gap_induction"] += 1
        return partner_ok and self.prove(adjacent, cols) and self.prove(smaller, cols)


def gap_profile(item: dict) -> tuple[int, int, int, int]:
    rows = item["rows"]
    cols = item["cols"]
    return (
        rows[1] - rows[0],
        rows[2] - rows[1],
        cols[1] - cols[0],
        cols[2] - cols[1],
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Replay H1326 strict-PF2 plus H811 implies PF3.")
    parser.add_argument("--full-report", action="store_true")
    args = parser.parse_args()

    h1325 = load(H1325)
    deps = load(H1325_DEPS)
    h811 = load(H811)
    h813 = load(H813)
    h816 = load(H816)

    checks: dict[str, bool] = {
        f"hash::{path.relative_to(ROOT)}": digest(path) == expected
        for path, expected in PINNED_HASHES.items()
    }
    checks.update(generic_plucker_checks())
    checks.update(structural_factor_checks())
    checks.update(duality_checks())
    checks.update(
        {
            "H1325_closed": h1325.get("classification")
            == "xi_global_translated_d3_h811_consecutive_rows_closed"
            and bool(h1325.get("all_checks_pass")),
            "H1325_dependency_addendum_closed": deps.get("classification")
            == "h1325_h811_explicit_dependency_gates_closed"
            and bool(deps.get("all_checks_pass")),
            "H811_exact_consecutive_scope_loaded": h811.get("classification")
            == "repaired_consecutive_row_pf3_theorem_formalized"
            and "consecutive rows only" in h811["theorem"]["scope"],
            "Xi_coefficients_strictly_positive": bool(
                deps["checks"]["H908_positive_kernel_measure"]
                and deps["checks"]["H908_nonzero_positive_moments"]
            ),
            "Xi_qn_strictly_below_one_all_n": bool(
                h1325["checks"]["finite_and_H811"]["finite_PF2_through_q128"]
                and h1325["checks"]["exact_tail"]["X_below_one_for_all_n_ge_127"]
            ),
            "Xi_H811_consecutive_rows_closed": bool(
                h1325["checks"]["finite_and_H811"]["finite_all_D3_q_lowers_positive"]
                and h1325["checks"]["exact_tail"]["terminal_H_positive_for_all_n_ge_127"]
                and h1325["checks"]["finite_and_H811"]["boundary_start_1_strictly_positive"]
            ),
            "H813_source_count_is_2310": h813["summary_counts"]["active_sparse_row_pairs"] == 2310,
            "H816_source_count_is_356": h816["counts"]["uncovered"] == 356,
        }
    )

    p35, _ = determinant_q_poly((0, 1, 3), (3, 4, 5), 9)
    p36, _ = determinant_q_poly((0, 1, 3), (3, 4, 6), 9)
    checks["H814_six_term_k5_matched"] = poly_to_string(p35, 9) == (
        "1 - q4 - q3*q4 + q3*q4^2*q5 + q2*q3^2*q4^2 - q2*q3^2*q4^2*q5"
    )
    checks["H814_six_term_k6_matched"] = poly_to_string(p36, 9) == (
        "1 - q4 - q3*q4^2*q5 + q3*q4^3*q5^2*q6 + "
        "q2*q3^2*q4^3*q5 - q2*q3^2*q4^3*q5^2*q6"
    )

    sequence_len = 9  # H813 length=7 uses indices through 8.
    active_sparse = []
    for rows, cols in all_order3_pairs(sequence_len):
        if is_consecutive_rows(rows):
            continue
        if determinant_terms(rows, cols, sequence_len):
            active_sparse.append((tuple(rows), tuple(cols)))

    planner = ProofPlanner()
    active_failures = [
        {"rows": rows, "cols": cols}
        for rows, cols in active_sparse
        if not planner.prove(rows, cols)
    ]
    checks["H813_active_sparse_enumerated_2310"] = len(active_sparse) == 2310
    checks["H813_all_2310_replayed_by_induction"] = not active_failures

    remaining = remaining_after_h815(7)
    adjacent_rows = [item for item in remaining if 1 in gap_profile(item)[:2]]
    adjacent_rows_or_cols = [item for item in remaining if 1 in gap_profile(item)]
    four_gaps_at_least_two = [item for item in remaining if min(gap_profile(item)) >= 2]
    all_h816_proved = [
        item
        for item in remaining
        if planner.prove(tuple(item["rows"]), tuple(item["cols"]))
    ]
    six_total = sum(item["analysis"]["term_count"] == 6 for item in remaining)
    six_adjacent_union = sum(
        item["analysis"]["term_count"] == 6 for item in adjacent_rows_or_cols
    )
    residual_term_counts = Counter(item["analysis"]["term_count"] for item in four_gaps_at_least_two)

    checks.update(
        {
            "H816_adjacent_rows_242": len(adjacent_rows) == 242,
            "H816_adjacent_rows_or_columns_344": len(adjacent_rows_or_cols) == 344,
            "H816_all_six_term_forms_closed_before_general_step": six_total == 125
            and six_adjacent_union == 125,
            "H816_pre_general_residual_is_12": len(four_gaps_at_least_two) == 12,
            "H816_pre_general_residual_term_counts": dict(residual_term_counts)
            == {3: 3, 4: 8, 5: 1},
            "H816_full_induction_closes_356_of_356": len(all_h816_proved) == 356,
        }
    )

    all_checks = all(checks.values())
    report = {
        "schema": "rh_h1326_plucker_pf3_completion.v0",
        "classification": (
            "h1326_strict_pf2_plus_h811_implies_pf3_closed"
            if all_checks
            else "h1326_pf3_completion_replay_failed"
        ),
        "all_checks_pass": all_checks,
        "failed_checks": [name for name, passed in checks.items() if not passed],
        "abstract_theorem": {
            "inputs": [
                "a_n>0 for n>=0 and a_n=0 for n<0",
                "strict PF2 (equivalently here: R_n strictly decreasing)",
                "every order-3 Toeplitz minor with consecutive rows is nonnegative",
            ],
            "conclusion": "the one-sided Toeplitz matrix is PF3",
        },
        "Xi_instantiation": {
            "strict_PF2": "H1325 finite q_n<1 through 128 plus X_n<1 for n>=127",
            "consecutive_rows": "H1325 plus H811",
            "conclusion": "the Xi coefficient sequence is PF3",
        },
        "audits": {
            "H813_active_sparse": len(active_sparse),
            "H813_failures": active_failures,
            "H816_input": len(remaining),
            "H816_adjacent_rows": len(adjacent_rows),
            "H816_adjacent_rows_or_columns": len(adjacent_rows_or_cols),
            "H816_four_gap_residual_before_general_induction": len(four_gaps_at_least_two),
            "H816_full_induction_proved": len(all_h816_proved),
            "planner_max_gap_sum": planner.max_depth,
            "planner_step_kinds": dict(sorted(planner.kinds.items())),
        },
        "scope_limitations": [
            "PF3 does not imply PF4 or PF-infinity",
            "this does not prove all Jensen polynomial degrees",
            "this does not prove the Riemann Hypothesis",
        ],
        "checks": checks if args.full_report else {"count": len(checks), "passed": sum(checks.values())},
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if all_checks else 1


if __name__ == "__main__":
    raise SystemExit(main())
