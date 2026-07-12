#!/usr/bin/env python3
"""Exact dependency assembly for the global gamma-log moment box (H1315).

The verifier joins the computer-assisted low band from H941 to the audited
H1314 uniform pair cap and the H1313 global even-moment certificate.  Exact
rational arithmetic then closes the two H925 score budgets.  Finally, the
published H922a, H927, and H922 lemmas assemble the normalized moment box and
the one-sided third-cumulant estimate.

This script reruns H1313.  It attests the published H1314 theorem statement;
it does not recreate H1314's still-unpublished polynomial/Arb verifier.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import rh_h1313_h946_global_moment_assembly as h1313


ROOT = Path(__file__).resolve().parents[1]
RIEMANN = ROOT / "research" / "riemann"
DEFAULT_OUT = RIEMANN / "h1315_gamma_log_moment_box_global_assembly.json"


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def fraction_payload(value: Fraction) -> dict[str, Any]:
    return {
        "fraction": fraction_text(value),
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": format(float(value), ".17g"),
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def note_attestation(path: Path, required_fragments: list[str]) -> dict[str, Any]:
    exists = path.exists()
    text = path.read_text(encoding="utf-8") if exists else ""
    fragments = [
        {"fragment": fragment, "found": fragment in text}
        for fragment in required_fragments
    ]
    return {
        "path": str(path.relative_to(ROOT)),
        "exists": exists,
        "sha256": sha256(path) if exists else None,
        "required_fragments": fragments,
        "passed": exists and all(row["found"] for row in fragments),
    }


def dependency_certificate() -> dict[str, Any]:
    h1313_rerun = h1313.certificate()
    h1313_bounds = h1313_rerun.get("statement", {}).get("bounds", {})
    h1313_rerun_pass = (
        bool(h1313_rerun.get("all_checks_pass"))
        and h1313_rerun.get("statement", {}).get("range") == "r>=5/2"
        and h1313_bounds == {"M2": "26/25", "M4": "7/2", "M6": "21"}
    )

    notes = {
        "H1313": note_attestation(
            RIEMANN / "h1313_h946_global_even_moment_bounds.md",
            [
                "h1313_h946_global_even_moment_subproblem_proved",
                "every `r>=5/2`",
                "E_r[W^2] <= 26/25",
                "E_r[W^4] <= 7/2",
                "E_r[W^6] <= 21",
            ],
        ),
        "H1314": note_attestation(
            RIEMANN / "h1314_uniform_pair_cap_proved_audited.md",
            [
                "h1314_uniform_pair_cap_proved_audited_pending_verifier_persistence",
                "R_r(t) <= 51/100",
                "E[W^2 K]",
                "663/1250",
                "14229/5000",
            ],
        ),
        "H941": note_attestation(
            RIEMANN / "h941_gamma_log_score_adaptive_cover_to_2p50.md",
            [
                "h941_gamma_log_score_adaptive_cover_to_2p50_proved_local_computer_assisted",
                "2 <= r <= 2.50",
                "gamma-log H928 lower-sign cover",
                "Combining H925, H927, H922a, H922",
            ],
        ),
        "H925": note_attestation(
            RIEMANN / "h925_score_excess_budget_reduction.md",
            [
                "h925_score_excess_budget_reduction_proved_budget_open",
                "E[W] = -E[N(W)]",
                "E[W^3]",
                "E[N(W)] <= a",
                "E[(W^2+2)N(W)] <= 3a",
            ],
        ),
        "H922a": note_attestation(
            RIEMANN / "h922a_gamma_log_mode_asymmetry_signs.md",
            [
                "h922a_gamma_log_mode_asymmetry_signs_proved",
                "E[W] <= 0",
                "E[W^3] <= 0",
            ],
        ),
        "H927": note_attestation(
            RIEMANN / "h927_gamma_log_a_le_one_low_threshold.md",
            [
                "h927_gamma_log_a_le_one_r_ge_2_proved",
                "0<a<=1",
                "for every `r>=2`",
            ],
        ),
        "H922": note_attestation(
            RIEMANN / "h922_gamma_log_normalized_moment_box_gate.md",
            [
                "h922_normalized_moment_box_gate_proved_moment_box_open",
                "B/A^3 <= 7/(q^2 r_q)",
                "K_alpha'''(q) >= -35/(q^2 r_q)",
                "alpha=9/4",
                "r_q>=2",
            ],
        ),
    }

    return {
        "H1313": {
            "note": notes.pop("H1313"),
            "rerun": {
                "id": h1313_rerun.get("id"),
                "range": h1313_rerun.get("statement", {}).get("range"),
                "bounds": h1313_bounds,
                "dependency_pass": bool(h1313_rerun.get("dependency_pass")),
                "all_checks_pass": bool(h1313_rerun.get("all_checks_pass")),
                "passed": bool(h1313_rerun_pass),
            },
            "passed": bool(notes.get("H1313", {}).get("passed")),
        },
        **notes,
    }


def exact_arithmetic_certificate() -> dict[str, Any]:
    pair_cap = Fraction(51, 100)
    moment_2 = Fraction(26, 25)
    moment_4 = Fraction(7, 2)

    first_budget = pair_cap * moment_2
    third_budget = pair_cap * (moment_4 + 2 * moment_2)
    expected_first = Fraction(663, 1250)
    expected_third = Fraction(14229, 5000)

    first_pass = first_budget == expected_first and first_budget < 1
    third_pass = third_budget == expected_third and third_budget < 3

    return {
        "inputs": {
            "uniform_pair_cap": fraction_payload(pair_cap),
            "M2_upper": fraction_payload(moment_2),
            "M4_upper": fraction_payload(moment_4),
        },
        "first_score_budget": {
            "derivation": "(51/100)*(26/25)",
            "value": fraction_payload(first_budget),
            "expected": fraction_payload(expected_first),
            "target": "<1",
            "passed": bool(first_pass),
        },
        "cubic_score_budget": {
            "derivation": "(51/100)*(7/2+2*(26/25))",
            "value": fraction_payload(third_budget),
            "expected": fraction_payload(expected_third),
            "target": "<3",
            "passed": bool(third_pass),
        },
        "all_checks_pass": bool(first_pass and third_pass),
    }


def certificate() -> dict[str, Any]:
    dependencies = dependency_certificate()
    arithmetic = exact_arithmetic_certificate()

    # H1313 is special because its compact certificate is rerun here.
    h1313_note_pass = dependencies["H1313"]["note"]["passed"]
    h1313_rerun_pass = dependencies["H1313"]["rerun"]["passed"]
    dependencies["H1313"]["passed"] = bool(
        h1313_note_pass and h1313_rerun_pass
    )

    required_dependency_ids = (
        "H1313",
        "H1314",
        "H941",
        "H925",
        "H922a",
        "H927",
        "H922",
    )
    dependencies_pass = all(
        dependencies[name]["passed"] for name in required_dependency_ids
    )

    low_band_pass = dependencies["H941"]["passed"]
    high_band_pass = (
        dependencies["H1313"]["passed"]
        and dependencies["H1314"]["passed"]
        and arithmetic["all_checks_pass"]
    )
    coverage_pass = low_band_pass and high_band_pass

    moment_box_pass = (
        coverage_pass
        and dependencies["H925"]["passed"]
        and dependencies["H922a"]["passed"]
        and dependencies["H927"]["passed"]
    )

    centered_moment_constant = Fraction(3) + Fraction(2)
    coefficient_constant = Fraction(7)
    final_constant = centered_moment_constant * coefficient_constant
    cumulant_pass = (
        moment_box_pass
        and dependencies["H922"]["passed"]
        and centered_moment_constant == 5
        and final_constant == 35
    )

    all_checks_pass = dependencies_pass and cumulant_pass
    return {
        "id": "h1315_gamma_log_moment_box_global_assembly",
        "classification": "h1315_gamma_log_moment_box_global_assembly_proved_from_audited_dependencies",
        "statement": {
            "alpha": "9/4",
            "range": "r_q>=2",
            "moment_box": [
                "0<a<=1",
                "-a<=E[W]<=0",
                "E[W^3]>=-3a",
            ],
            "conclusion": "K_alpha'''(q)>=-35/(q^2*r_q)",
        },
        "dependencies": dependencies,
        "dependency_checks_pass": bool(dependencies_pass),
        "exact_arithmetic": arithmetic,
        "coverage": {
            "domain": "r>=2",
            "bands": [
                {
                    "left": "2",
                    "right": "5/2",
                    "dependency": "H941",
                    "role": "computer-assisted H925 lower-sign budgets",
                    "passed": bool(low_band_pass),
                },
                {
                    "left": "5/2",
                    "right": "infinity",
                    "dependency": "H1314+H1313+exact arithmetic",
                    "role": "uniform pair cap plus global even moments",
                    "passed": bool(high_band_pass),
                },
            ],
            "seam": "5/2",
            "no_gap": True,
            "passed": bool(coverage_pass),
        },
        "moment_box_assembly": {
            "H925_lower_signs": "score budgets imply E[W]>=-a and E[W^3]>=-3a",
            "H922a_upper_signs": "E[W]<=0 and E[W^3]<=0",
            "H927_cubic_coefficient": "0<a<=1",
            "passed": bool(moment_box_pass),
        },
        "cumulant_consequence": {
            "H922_centered_moment_factor": fraction_payload(
                centered_moment_constant
            ),
            "H922_H921_coefficient_factor": fraction_payload(
                coefficient_constant
            ),
            "final_constant": fraction_payload(final_constant),
            "conclusion": "K_alpha'''(q)>=-35/(q^2*r_q) for alpha=9/4 and r_q>=2",
            "passed": bool(cumulant_pass),
        },
        "verification_scope": {
            "rerun": "H1313 compact dependency certificate",
            "attested": "H1314, H941, H925, H922a, H927, H922 published statements",
            "limitation": "does not recreate H1314's polynomial/Arb proof",
        },
        "does_not_prove": [
            "effective full-Xi transfer constants",
            "finite H920 closure",
            "all-degree Jensen hyperbolicity or total positivity",
            "the Riemann Hypothesis",
        ],
        "all_checks_pass": bool(all_checks_pass),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="run all checks and print JSON without writing the canonical file",
    )
    args = parser.parse_args()

    report = certificate()
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if not args.check_only:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0 if report["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
