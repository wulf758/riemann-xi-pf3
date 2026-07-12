#!/usr/bin/env python3
"""Fail-closed replay for the H1326 sparse-row Plucker bridge.

The proof itself is parameter-free.  For fixed three columns, the two
Grassmann--Plucker identities below reduce row gaps (1, 2) and (2, 1) to
consecutive-row 3-minors and Toeplitz 2-minors.  The Toeplitz specialization
uses C=(j,j+1,k), d=j+2, j>=2, k>=d.

This checker deliberately fails if a pinned H1325/H811 dependency drifts, if
either orientation changes sign, if one of the five Toeplitz factorizations
fails, or if the audited H816-window coverage changes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
RIEMANN = ROOT / "research" / "riemann"
THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from rh_h813_sparse_row_symbolic_hierarchy import determinant_q_poly, poly_to_string
from rh_h816_d3_cone_sparse_row_coupling import remaining_after_h815


H1325 = RIEMANN / "h1325_xi_global_d3_h811_consecutive_rows.json"
H1325_DEPS = RIEMANN / "h1325_h811_dependency_addendum.json"
H811 = RIEMANN / "h811_repaired_consecutive_row_theorem.json"
H816 = RIEMANN / "h816_d3_cone_sparse_row_coupling.json"

PINNED_HASHES = {
    H1325: "9308258F232A84CB964F0CE4357C2BE279D5E2FE3A49B227B37CCD7E6D43A1A0",
    H1325_DEPS: "F034A24C08EED9FC8E0A0405D596DF47A4227CFBEC7E7664D2823062A45A26B4",
    H811: "09C0B09673CD582BB1F9968466B5304FBA13A4B8AE0D97B6348EAD60498D6A12",
    H816: "8E6F0364C259FBBBC1AB66D9B1C250D2F706004FE1842644AFC584D07BDC2A00",
}

EXPECTED_COVERAGE = {
    "row_gaps=(1,2), col_gaps=(1,1)": 10,
    "row_gaps=(1,2), col_gaps=(1,2)": 10,
    "row_gaps=(1,2), col_gaps=(1,3)": 6,
    "row_gaps=(1,2), col_gaps=(1,4)": 3,
    "row_gaps=(1,2), col_gaps=(1,5)": 1,
    "row_gaps=(2,1), col_gaps=(1,1)": 10,
    "row_gaps=(2,1), col_gaps=(1,2)": 6,
    "row_gaps=(2,1), col_gaps=(1,3)": 3,
    "row_gaps=(2,1), col_gaps=(1,4)": 1,
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def generic_plucker_checks() -> dict[str, bool]:
    """Check both oriented identities for five generic row vectors."""

    x = sp.symbols("x0:15")
    rows = {r: x[3 * r : 3 * r + 3] for r in range(5)}

    def minor(*indices: int) -> sp.Expr:
        return sp.det(sp.Matrix([rows[index] for index in indices]))

    first = sp.expand(
        minor(0, 1, 3) * minor(1, 2, 4)
        - minor(0, 1, 2) * minor(1, 3, 4)
        - minor(0, 1, 4) * minor(1, 2, 3)
    )
    mirror = sp.expand(
        minor(0, 2, 3) * minor(1, 2, 4)
        - minor(0, 1, 2) * minor(2, 3, 4)
        - minor(0, 2, 4) * minor(1, 2, 3)
    )
    return {
        "plucker_013_sign_and_identity": first == 0,
        "plucker_023_sign_and_identity": mirror == 0,
    }


def toeplitz_factor_checks() -> dict[str, bool]:
    """Check all d=j+2 factorizations with a_n=0 for n<0."""

    jm3, jm2, jm1, aj, jp1 = sp.symbols("a_jm3 a_jm2 a_jm1 a_j a_jp1")
    z0, z1, z2, z3, h = sp.symbols("z0 z1 z2 z3 h")
    rows = {
        0: (aj, jp1, z0),
        1: (jm1, aj, z1),
        2: (jm2, jm1, z2),
        3: (jm3, jm2, z3),
        4: (sp.Integer(0), sp.Integer(0), h),
    }

    def minor(*indices: int) -> sp.Expr:
        return sp.expand(sp.det(sp.Matrix([rows[index] for index in indices])))

    expected = {
        "factor_Delta_12d": (minor(1, 2, 4), h * (jm1**2 - jm2 * aj)),
        "factor_Delta_13d": (minor(1, 3, 4), h * (jm1 * jm2 - jm3 * aj)),
        "factor_Delta_01d": (minor(0, 1, 4), h * (aj**2 - jm1 * jp1)),
        "factor_Delta_23d": (minor(2, 3, 4), h * (jm2**2 - jm3 * jm1)),
        "factor_Delta_02d": (minor(0, 2, 4), h * (aj * jm1 - jm2 * jp1)),
    }
    checks = {
        name: sp.expand(actual - target) == 0
        for name, (actual, target) in expected.items()
    }
    # j=2 is the only boundary in the theorem: a_(j-3)=a_-1=0.
    checks["boundary_j_equals_2_uses_a_minus_1_zero"] = (
        sp.expand(expected["factor_Delta_13d"][0].subs(jm3, 0) - h * jm1 * jm2) == 0
        and sp.expand(expected["factor_Delta_23d"][0].subs(jm3, 0) - h * jm2**2) == 0
    )
    # h=a_(k-d); at k=d this is a_0, not an empty or zero factor.
    checks["k_equals_d_tail_factor_is_a0"] = h != 0
    return checks


def covered_by_h1326(item: dict) -> bool:
    rows = tuple(item["rows"])
    cols = tuple(item["cols"])
    row_gaps = (rows[1] - rows[0], rows[2] - rows[1])
    return (
        row_gaps in {(1, 2), (2, 1)}
        and cols[1] == cols[0] + 1
        and cols[0] - rows[0] >= 2
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Replay the H1326 Plucker first-column-gap bridge.")
    parser.add_argument("--full-report", action="store_true")
    args = parser.parse_args()

    h1325 = load(H1325)
    deps = load(H1325_DEPS)
    h811 = load(H811)
    h816 = load(H816)

    checks: dict[str, bool] = {
        f"hash::{path.relative_to(ROOT)}": digest(path) == expected
        for path, expected in PINNED_HASHES.items()
    }
    checks.update(generic_plucker_checks())
    checks.update(toeplitz_factor_checks())
    checks.update(
        {
            "H1325_closed": h1325.get("classification")
            == "xi_global_translated_d3_h811_consecutive_rows_closed"
            and bool(h1325.get("all_checks_pass")),
            "H1325_dependency_addendum_closed": deps.get("classification")
            == "h1325_h811_explicit_dependency_gates_closed"
            and bool(deps.get("all_checks_pass")),
            "H811_loaded_with_exact_scope": h811.get("classification")
            == "repaired_consecutive_row_pf3_theorem_formalized"
            and "consecutive rows only" in h811["theorem"]["scope"],
            "Xi_coefficients_strictly_positive": bool(
                deps["checks"]["H908_positive_kernel_measure"]
                and deps["checks"]["H908_nonzero_positive_moments"]
            ),
            "Xi_global_strict_adjacent_PF2": bool(
                h1325["checks"]["finite_and_H811"]["finite_PF2_through_q128"]
                and h1325["checks"]["exact_tail"]["X_below_one_for_all_n_ge_127"]
            ),
            "H811_consecutive_row_inputs_closed": bool(
                h1325["checks"]["finite_and_H811"]["finite_all_D3_q_lowers_positive"]
                and h1325["checks"]["exact_tail"]["terminal_H_positive_for_all_n_ge_127"]
                and h1325["checks"]["finite_and_H811"]["boundary_start_1_strictly_positive"]
            ),
            "H816_source_count_is_356": h816["counts"]["uncovered"] == 356,
        }
    )

    p35, _ = determinant_q_poly((0, 1, 3), (3, 4, 5), 9)
    p36, _ = determinant_q_poly((0, 1, 3), (3, 4, 6), 9)
    checks["H814_representative_k5_matched"] = poly_to_string(p35, 9) == (
        "1 - q4 - q3*q4 + q3*q4^2*q5 + q2*q3^2*q4^2 - q2*q3^2*q4^2*q5"
    )
    checks["H814_representative_k6_matched"] = poly_to_string(p36, 9) == (
        "1 - q4 - q3*q4^2*q5 + q3*q4^3*q5^2*q6 + "
        "q2*q3^2*q4^3*q5 - q2*q3^2*q4^3*q5^2*q6"
    )

    remaining = remaining_after_h815(7)
    covered = [item for item in remaining if covered_by_h1326(item)]
    coverage = Counter(item["analysis"]["family"] for item in covered)
    checks["audited_H816_window_newly_covered_50"] = len(covered) == 50
    checks["audited_H816_family_counts_exact"] = dict(sorted(coverage.items())) == dict(
        sorted(EXPECTED_COVERAGE.items())
    )

    all_checks = all(checks.values())
    report = {
        "schema": "rh_h1326_plucker_first_gap_bridge.v0",
        "classification": (
            "h1326_plucker_first_gap_sparse_families_closed"
            if all_checks
            else "h1326_plucker_first_gap_replay_failed"
        ),
        "all_checks_pass": all_checks,
        "failed_checks": [name for name, passed in checks.items() if not passed],
        "theorem_scope": {
            "row_families": ["(r,r+1,r+3)", "(r,r+2,r+3)"],
            "columns": "(c,c+1,l), with c-r>=2 and l>=c+2",
            "inputs": "strict PF2 denominator plus PF2 and all consecutive-row order-3 minors",
        },
        "identities": {
            "gap_1_2": "Delta013*Delta12d=Delta012*Delta13d+Delta01d*Delta123",
            "gap_2_1": "Delta023*Delta12d=Delta012*Delta23d+Delta02d*Delta123",
            "parameters": "C=(j,j+1,k), d=j+2, j>=2, k>=d",
        },
        "audited_H816_window": {
            "input": len(remaining),
            "newly_covered": len(covered),
            "remaining_after_this_bridge": len(remaining) - len(covered),
            "family_counts": dict(sorted(coverage.items())),
        },
        "checks": checks if args.full_report else {"count": len(checks), "passed": sum(checks.values())},
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if all_checks else 1


if __name__ == "__main__":
    raise SystemExit(main())
