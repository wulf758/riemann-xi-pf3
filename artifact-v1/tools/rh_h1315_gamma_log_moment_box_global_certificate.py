#!/usr/bin/env python3
"""Canonical exact dependency certificate for H1315.

This standalone repository verifier reruns H1313, attests the published
H1314/H941/H925/H922a/H927/H922 dependency statements, performs the H1315
budget arithmetic with exact fractions, and writes the canonical H1315 JSON.

It deliberately does not recreate H1314's unpublished polynomial/Arb driver.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
from typing import Any

import rh_h1313_h946_global_moment_assembly as h1313


ROOT = Path(__file__).resolve().parents[1]
RIEMANN = ROOT / "research" / "riemann"
DEFAULT_OUT = RIEMANN / "h1315_gamma_log_moment_box_global_assembly.json"


def fraction_text(value: Fraction) -> str:
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def fraction_payload(value: Fraction) -> dict[str, Any]:
    with localcontext() as context:
        context.prec = 50
        decimal = format(
            (Decimal(value.numerator) / Decimal(value.denominator)).normalize(),
            "f",
        )
    return {
        "fraction": fraction_text(value),
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": decimal,
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def note_attestation(path: Path, fragments: list[str]) -> dict[str, Any]:
    exists = path.exists()
    text = path.read_text(encoding="utf-8") if exists else ""
    checks = [
        {"fragment": fragment, "found": fragment in text}
        for fragment in fragments
    ]
    return {
        "path": str(path.relative_to(ROOT)),
        "exists": exists,
        "sha256": sha256(path) if exists else None,
        "required_fragments": checks,
        "passed": exists and all(check["found"] for check in checks),
    }


def dependency_certificate() -> dict[str, Any]:
    h1313_run = h1313.certificate()
    h1313_bounds = h1313_run.get("statement", {}).get("bounds", {})
    h1313_rerun_pass = (
        bool(h1313_run.get("all_checks_pass"))
        and bool(h1313_run.get("dependency_pass"))
        and h1313_run.get("statement", {}).get("range") == "r>=5/2"
        and h1313_bounds == {"M2": "26/25", "M4": "7/2", "M6": "21"}
    )

    h1313_note = note_attestation(
        RIEMANN / "h1313_h946_global_even_moment_bounds.md",
        [
            "h1313_h946_global_even_moment_subproblem_proved",
            "every `r>=5/2`",
            "E_r[W^2] <= 26/25",
            "E_r[W^4] <= 7/2",
            "E_r[W^6] <= 21",
        ],
    )

    dependencies: dict[str, Any] = {
        "H1313": {
            "note": h1313_note,
            "rerun": {
                "id": h1313_run.get("id"),
                "range": h1313_run.get("statement", {}).get("range"),
                "bounds": h1313_bounds,
                "dependency_pass": bool(h1313_run.get("dependency_pass")),
                "all_checks_pass": bool(h1313_run.get("all_checks_pass")),
                "passed": bool(h1313_rerun_pass),
            },
            "passed": bool(h1313_note["passed"] and h1313_rerun_pass),
        },
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
                "computer-assisted gamma-log H928 lower-sign",
                "Combining H925, H927, H922a, H922",
            ],
        ),
        "H925": note_attestation(
            RIEMANN / "h925_score_excess_budget_reduction.md",
            [
                "h925_score_excess_budget_reduction_proved_budget_open",
                "E[W] = -E[N(W)]",
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
    return dependencies


def arithmetic_certificate() -> dict[str, Any]:
    cap = Fraction(51, 100)
    m2 = Fraction(26, 25)
    m4 = Fraction(7, 2)
    first = cap * m2
    cubic = cap * (m4 + 2 * m2)
    expected_first = Fraction(663, 1250)
    expected_cubic = Fraction(14229, 5000)

    first_pass = first == expected_first and first < 1
    cubic_pass = cubic == expected_cubic and cubic < 3
    return {
        "inputs": {
            "uniform_pair_cap": fraction_payload(cap),
            "M2_upper": fraction_payload(m2),
            "M4_upper": fraction_payload(m4),
        },
        "first_score_budget": {
            "identity": "E[N]/a=E[W^2*K]",
            "derivation": "(51/100)*(26/25)",
            "value": fraction_payload(first),
            "expected": fraction_payload(expected_first),
            "target": "<1",
            "passed": bool(first_pass),
        },
        "cubic_score_budget": {
            "identity": "E[(W^2+2)N]/a=E[(W^2+2)W^2*K]",
            "derivation": "(51/100)*(7/2+2*(26/25))",
            "value": fraction_payload(cubic),
            "expected": fraction_payload(expected_cubic),
            "target": "<3",
            "passed": bool(cubic_pass),
        },
        "all_checks_pass": bool(first_pass and cubic_pass),
    }


def certificate() -> dict[str, Any]:
    dependencies = dependency_certificate()
    arithmetic = arithmetic_certificate()
    dependency_ids = ("H1313", "H1314", "H941", "H925", "H922a", "H927", "H922")
    dependencies_pass = all(dependencies[name]["passed"] for name in dependency_ids)

    low_left, low_right = Fraction(2), Fraction(5, 2)
    high_left = Fraction(5, 2)
    no_gap = low_left == 2 and low_right == high_left
    low_pass = dependencies["H941"]["passed"]
    high_pass = (
        dependencies["H1313"]["passed"]
        and dependencies["H1314"]["passed"]
        and arithmetic["all_checks_pass"]
    )
    coverage_pass = no_gap and low_pass and high_pass

    moment_box_pass = (
        coverage_pass
        and dependencies["H925"]["passed"]
        and dependencies["H922a"]["passed"]
        and dependencies["H927"]["passed"]
    )
    centered_factor = Fraction(3) + Fraction(2)
    coefficient_factor = Fraction(7)
    final_factor = centered_factor * coefficient_factor
    cumulant_pass = (
        moment_box_pass
        and dependencies["H922"]["passed"]
        and centered_factor == 5
        and final_factor == 35
    )

    return {
        "id": "h1315_gamma_log_moment_box_global_assembly",
        "classification": "h1315_gamma_log_moment_box_global_assembly_proved_from_audited_dependencies",
        "canonical_verifier": "tools/rh_h1315_gamma_log_moment_box_global_certificate.py",
        "statement": {
            "alpha": "9/4",
            "range": "r_q>=2",
            "moment_box": ["0<a<=1", "-a<=E[W]<=0", "E[W^3]>=-3a"],
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
                    "passed": bool(low_pass),
                },
                {
                    "left": "5/2",
                    "right": "infinity",
                    "dependency": "H1314+H1313+exact arithmetic",
                    "passed": bool(high_pass),
                },
            ],
            "seam": fraction_text(low_right),
            "no_gap": bool(no_gap),
            "passed": bool(coverage_pass),
        },
        "moment_box_assembly": {
            "H925_lower_signs": "E[W]>=-a and E[W^3]>=-3a",
            "H922a_upper_signs": "E[W]<=0 and E[W^3]<=0",
            "H927_cubic_coefficient": "0<a<=1",
            "passed": bool(moment_box_pass),
        },
        "cumulant_consequence": {
            "centered_moment_factor": fraction_payload(centered_factor),
            "coefficient_factor": fraction_payload(coefficient_factor),
            "final_constant": fraction_payload(final_factor),
            "conclusion": "K_alpha'''(q)>=-35/(q^2*r_q) for alpha=9/4 and r_q>=2",
            "passed": bool(cumulant_pass),
        },
        "verification_scope": {
            "rerun": "H1313 compact Arb/dependency certificate",
            "attested": "H1314, H941, H925, H922a, H927, H922 theorem statements",
            "limitation": "does not recreate H1314's polynomial/Arb proof",
        },
        "does_not_prove": [
            "effective full-Xi transfer constants",
            "finite H920 closure",
            "all-degree Jensen hyperbolicity or total positivity",
            "the Riemann Hypothesis",
        ],
        "all_checks_pass": bool(dependencies_pass and cumulant_pass),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--check-only", action="store_true")
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
