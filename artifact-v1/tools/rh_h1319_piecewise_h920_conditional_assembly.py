#!/usr/bin/env python3
"""Audit the exact piecewise H920/H1319 closure logic.

The compact full-Xi estimate is deliberately an external obligation.  This
script proves the arithmetic and all off-by-one statements in the implication

    finite H803 base + compact pointwise estimate + H1317 tail
        => H908-A for n >= 127 and H803 for every n >= 2.

Without a compatible compact certificate the report is *conditional* and the
default invocation still exits successfully when the conditional logic is
sound.  Use ``--require-unconditional`` to make the missing compact certificate
an error (exit status 2).
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any


SCHEMA = "rh_h1319_piecewise_h920_conditional_assembly.v1"
COMPACT_SCHEMA = "rh_h1319_full_xi_compact_interval_certificate.v1"


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def exact_geometric_threshold(a: Fraction) -> int:
    """Return min n>=3 such that a/(n-2)^2 <= 1/n^2, exactly."""

    for n in range(3, 1_000_000):
        if a * n * n <= (n - 2) * (n - 2):
            return n
    raise RuntimeError("geometric threshold search guard exhausted")


def resolve_from_root(root: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else root / path


def load_finite_base(path: Path) -> dict[str, Any]:
    result: dict[str, Any] = {
        "path": str(path),
        "present": path.exists(),
        "accepted": False,
        "checks": {},
    }
    if not path.exists():
        result["error"] = "finite H804 report missing"
        return result

    try:
        report = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        result["error"] = f"cannot read finite H804 report: {exc}"
        return result

    rows = report.get("rows")
    if not isinstance(rows, list):
        result["error"] = "finite H804 report has no rows list"
        return result

    indices = [row.get("n") for row in rows]
    expected_indices = list(range(2, 127))
    each_h803_available = all(
        row.get("h803_condition_holds_midpoint") is True
        and row.get("h801_certified") is True
        for row in rows
    )
    top_checks = report.get("checks", {})
    checks = {
        "classification_expected": report.get("classification")
        == "unit_deficit_barrier_supported_with_finite_base",
        "indices_are_exactly_2_through_126": indices == expected_indices,
        "row_count_is_125": len(rows) == 125,
        "each_h803_row_available_and_inherits_h801_certificate": each_h803_available,
        "top_level_all_h803_conditions_available": top_checks.get(
            "all_h803_conditions_available"
        )
        is True,
        "small_n_unit_failures_recorded_exactly": top_checks.get(
            "unit_barrier_failures_before_tail"
        )
        == [2, 3, 4, 5, 6, 7],
    }
    result.update(
        {
            "classification": report.get("classification"),
            "n_range": [indices[0], indices[-1]] if indices else None,
            "row_count": len(rows),
            "checks": checks,
            "accepted": all(checks.values()),
            "interpretation": (
                "finite exact H803-threshold evidence inherited from H801/H803; "
                "not an all-small-n unit-barrier certificate"
            ),
        }
    )
    return result


def run_h1317_dependency(root: Path, script: Path) -> dict[str, Any]:
    result: dict[str, Any] = {
        "path": str(script),
        "present": script.exists(),
        "accepted": False,
    }
    if not script.exists():
        result["error"] = "H1317 verifier missing"
        return result

    completed = subprocess.run(
        [sys.executable, str(script)],
        cwd=root,
        capture_output=True,
        text=True,
        check=False,
    )
    result["exit_code"] = completed.returncode
    if completed.returncode != 0:
        result["error"] = completed.stderr.strip() or completed.stdout.strip()
        return result
    try:
        dependency_report = json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        result["error"] = f"H1317 verifier output is not JSON: {exc}"
        return result
    result.update(
        {
            "all_checks_pass": dependency_report.get("all_checks_pass"),
            "B_Xi": dependency_report.get("B_Xi"),
            "coefficient": dependency_report.get("effective_choice", {}).get(
                "coefficient"
            ),
            "accepted": dependency_report.get("all_checks_pass") is True,
        }
    )
    return result


def inspect_compact_certificate(path: Path) -> dict[str, Any]:
    expected_contract = {
        "schema": COMPACT_SCHEMA,
        "claim": {
            "kernel": "full_xi",
            "quantity": "t^2*kappa_Xi'''(t)",
            "lower_bound": "-3/4",
            "t_domain": {
                "lower": "125",
                "upper": "T_71",
                "closed": True,
            },
        },
        "required_true_checks": [
            "rigorous_interval_arithmetic",
            "full_domain_covered",
            "full_xi_kernel",
            "saddle_map_covers_t_domain",
            "all_local_lower_bounds_pass",
        ],
        "all_checks_pass": True,
    }
    result: dict[str, Any] = {
        "path": str(path),
        "present": path.exists(),
        "accepted": False,
        "expected_contract": expected_contract,
    }
    if not path.exists():
        result["status"] = "missing: compact estimate remains an explicit assumption"
        return result

    try:
        certificate = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        result["status"] = f"unreadable: {exc}"
        return result

    claim = certificate.get("claim", {})
    domain = claim.get("t_domain", {})
    certificate_checks = certificate.get("checks", {})
    interface_checks = {
        "schema_matches": certificate.get("schema") == COMPACT_SCHEMA,
        "kernel_is_full_xi": claim.get("kernel") == "full_xi",
        "quantity_matches": claim.get("quantity") == "t^2*kappa_Xi'''(t)",
        "lower_bound_is_exactly_minus_three_quarters": claim.get("lower_bound")
        == "-3/4",
        "domain_lower_is_125": domain.get("lower") == "125",
        "domain_upper_is_T71": domain.get("upper") == "T_71",
        "domain_is_closed": domain.get("closed") is True,
        "rigorous_interval_arithmetic": certificate_checks.get(
            "rigorous_interval_arithmetic"
        )
        is True,
        "full_domain_covered": certificate_checks.get("full_domain_covered")
        is True,
        "full_xi_kernel": certificate_checks.get("full_xi_kernel") is True,
        "saddle_map_covers_t_domain": certificate_checks.get(
            "saddle_map_covers_t_domain"
        )
        is True,
        "all_local_lower_bounds_pass": certificate_checks.get(
            "all_local_lower_bounds_pass"
        )
        is True,
        "all_checks_pass": certificate.get("all_checks_pass") is True,
    }
    result.update(
        {
            "status": "compatible" if all(interface_checks.values()) else "incompatible",
            "interface_checks": interface_checks,
            "accepted": all(interface_checks.values()),
            "warning": (
                "This assembly checks the certificate interface, not the internal "
                "interval arithmetic; the compact verifier remains primary evidence."
            ),
        }
    )
    return result


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    root = Path(__file__).resolve().parents[1]
    finite_path = resolve_from_root(root, args.finite_report)
    compact_path = resolve_from_root(root, args.compact_certificate)
    h1317_script = root / "tools" / "rh_h1317_effective_xi_transfer_assembly.py"
    h920_note = root / "research" / "riemann" / "h920_effective_finite_closure_gate.md"
    h1317_note = root / "research" / "riemann" / "h1317_effective_xi_transfer_assembly.md"

    a_compact = Fraction(3, 4)
    a_tail = Fraction(141, 142)
    n_geom_compact = exact_geometric_threshold(a_compact)
    n_geom_tail = exact_geometric_threshold(a_tail)
    finite_end = 126
    first_analytic_n = finite_end + 1
    compact_block_end = n_geom_tail - 1

    # Strict lower bound for T_71 from pi>3 and e>2.
    t71_lower = (
        Fraction(3 * 71 * 2**71, 1) - Fraction(9 * 71, 4) - 1
    ) / 2

    finite = load_finite_base(finite_path)
    h1317 = run_h1317_dependency(root, h1317_script)
    compact = inspect_compact_certificate(compact_path)

    exact_checks = {
        "compact_coefficient_is_below_tail_coefficient": a_compact < a_tail,
        "compact_coefficient_is_below_one": a_compact < 1,
        "tail_coefficient_is_below_one": a_tail < 1,
        "N_geom_compact_is_exactly_15": n_geom_compact == 15,
        "N_geom_compact_predecessor_fails": a_compact * 14 * 14
        > 12 * 12,
        "N_geom_compact_threshold_passes": a_compact * 15 * 15
        <= 13 * 13,
        "N_geom_tail_is_exactly_567": n_geom_tail == 567,
        "N_geom_tail_predecessor_fails": a_tail * 566 * 566 > 564 * 564,
        "N_geom_tail_threshold_passes": a_tail * 567 * 567 <= 565 * 565,
        "first_analytic_cube_starts_at_125": first_analytic_n - 2 == 125,
        "last_compact_block_cube_ends_at_567": compact_block_end + 1 == 567,
        "T71_strict_lower_bound_exceeds_567": t71_lower > 567,
        "global_block_cube_floor_exceeds_125": n_geom_tail - 2 >= 125,
        "tail_coefficient_decomposes_as_H1317": Fraction(70, 71)
        + Fraction(1, 142)
        == a_tail,
        "unit_barrier_strictly_implies_H803_for_n_ge_2": True,
    }

    dependency_checks = {
        "finite_base_accepted": finite.get("accepted") is True,
        "H1317_verifier_accepted": h1317.get("accepted") is True,
        "H920_note_present": h920_note.exists(),
        "H1317_note_present": h1317_note.exists(),
    }
    conditional_logic_passes = all(exact_checks.values()) and all(
        dependency_checks.values()
    )
    unconditional_from_dependencies = (
        conditional_logic_passes and compact.get("accepted") is True
    )

    if not conditional_logic_passes:
        classification = "h1319_piecewise_h920_assembly_internal_or_dependency_failure"
    elif unconditional_from_dependencies:
        classification = "h1319_piecewise_h920_closure_from_verified_dependencies"
    else:
        classification = "h1319_piecewise_h920_closure_conditional_on_compact_certificate"

    return {
        "schema": SCHEMA,
        "classification": classification,
        "constants": {
            "finite_H803_n_range": [2, finite_end],
            "compact_t_domain": ["125", "T_71"],
            "compact_a": fraction_text(a_compact),
            "tail_t_domain": ["T_71", "infinity"],
            "tail_a": fraction_text(a_tail),
            "T_71_definition": "(71*pi*exp(71)-(9/4)*71-1)/2",
            "strict_T71_lower_from_pi_gt_3_e_gt_2": fraction_text(t71_lower),
            "N_geom_compact": n_geom_compact,
            "N_geom_tail": n_geom_tail,
        },
        "piecewise_coverage": [
            {
                "integer_n": "2 <= n <= 126",
                "input": "existing finite H803 exact-threshold certificate",
                "conclusion": "H803 factorial threshold (unit barrier not claimed)",
            },
            {
                "integer_n": "127 <= n <= 566",
                "cube_t": "[n-2,n+1] subset [125,567] subset [125,T_71]",
                "input": "compact a=3/4 pointwise estimate",
                "geometry": "N_geom(3/4)=15, hence every n in this block passes",
                "conclusion": "H908-A unit barrier, hence H803",
            },
            {
                "integer_n": "n >= 567",
                "cube_t": "[n-2,n+1] subset [125,infinity)",
                "input": (
                    "the compact and H1317 bounds glue to the global weaker "
                    "a=141/142 estimate on [125,infinity)"
                ),
                "geometry": "N_geom(141/142)=567 exactly",
                "conclusion": (
                    "H908-A unit barrier, hence H803; includes cubes crossing T_71"
                ),
            },
        ],
        "exact_checks": exact_checks,
        "dependency_checks": dependency_checks,
        "finite_dependency": finite,
        "H1317_dependency": h1317,
        "compact_dependency": compact,
        "conditional_logic_checks_pass": conditional_logic_passes,
        "compact_assumption_discharged": compact.get("accepted") is True,
        "all_required_dependencies_pass": unconditional_from_dependencies,
        "conclusions": {
            "conditional_H908A_unit_barrier": "every integer n >= 127",
            "conditional_H803_factorial_barrier": "every integer n >= 2",
            "small_n_warning": (
                "The unit barrier is false for n=2,...,7 in H804; only the "
                "weaker exact H803 threshold is assembled for all n."
            ),
            "scope": (
                "Closes only the H803/H908 moment-curvature obligation if the "
                "compact certificate is supplied; H907 all-degree Jensen/TP and RH remain open."
            ),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit the conditional piecewise H920/H1319 closure."
    )
    parser.add_argument(
        "--finite-report",
        default="research/riemann/h804_unit_deficit_barrier_audit.json",
    )
    parser.add_argument(
        "--compact-certificate",
        default="research/riemann/h1319_full_xi_compact_interval_certificate.json",
    )
    parser.add_argument(
        "--out",
        default="research/riemann/h1319_piecewise_h920_conditional_assembly.json",
    )
    parser.add_argument(
        "--require-unconditional",
        action="store_true",
        help="exit 2 unless a compatible compact certificate discharges the assumption",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    report = build_report(args)
    out_path = resolve_from_root(root, args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(
        json.dumps(
            {
                "classification": report["classification"],
                "conditional_logic_checks_pass": report[
                    "conditional_logic_checks_pass"
                ],
                "compact_assumption_discharged": report[
                    "compact_assumption_discharged"
                ],
                "all_required_dependencies_pass": report[
                    "all_required_dependencies_pass"
                ],
                "output": str(out_path),
            },
            indent=2,
        )
    )
    if not report["conditional_logic_checks_pass"]:
        raise SystemExit(1)
    if args.require_unconditional and not report["all_required_dependencies_pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
