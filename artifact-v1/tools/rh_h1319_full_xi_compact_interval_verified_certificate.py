#!/usr/bin/env python3
"""Final audited entry point for the H1319 compact interval certificate.

This wrapper has two jobs beyond the component manifest assembler:

1. certify saddle monotonicity without a dependency-destroying giant r-ball;
2. require the independent high-precision audit of the score mean-value
   engine, its H1272 dependency, and the three detailed middle-band reports.

Only this audited wrapper writes the final compact-schema artifact consumed by
the strict H920 assembly.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

from flint import arb, ctx

import rh_h1319_full_xi_compact_interval_manifest as implementation


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_AUDIT = (
    ROOT
    / "research"
    / "riemann"
    / "h1319_score_mean_value_independent_audit_canonical.json"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def qarb(value: Fraction) -> arb:
    return arb(value.numerator) / arb(value.denominator)


def normalized(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT.resolve())).replace("\\", "/")


def inspect_independent_audit(
    path: Path, compact_report: dict[str, Any]
) -> dict[str, Any]:
    audit = json.loads(path.read_text(encoding="utf-8"))
    audit_meta = audit.get("audit", {})
    engine = audit.get("engine", {})
    dependency = audit.get("H1272_dependency", {})
    union = audit.get("union", {})
    audited_reports = audit.get("reports", [])
    mean_segments = [
        segment
        for segment in compact_report["r_proof_chain"]["segments"]
        if segment["kind"] == "score_mean_value_arb"
    ]
    audited_by_path = {
        row.get("path"): row for row in audited_reports if isinstance(row, dict)
    }

    engine_path = ROOT / str(engine.get("path", ""))
    runner_path = ROOT / str(audit_meta.get("canonical_runner", ""))
    h1272_path = ROOT / str(dependency.get("path", ""))
    detailed_hashes_match = True
    for segment in mean_segments:
        row = audited_by_path.get(segment["path"])
        detailed_hashes_match = detailed_hashes_match and bool(row)
        if row:
            detailed_hashes_match = (
                detailed_hashes_match
                and row.get("sha256") == segment["sha256"]
                and row.get("all_pass") is True
            )

    checks = {
        "audit_id": audit.get("id")
        == "h1319_score_mean_value_independent_audit_canonical",
        "audit_classification": audit.get("classification")
        == "independent_audit_pass_not_a_new_curvature_certificate",
        "audit_declared_canonical": audit.get("canonical") is True,
        "audit_all_checks_pass": audit.get("all_checks_pass") is True,
        "audit_runner_exit_zero": audit_meta.get("exit_code") == 0,
        "audit_runner_all_checks_pass": audit_meta.get("all_checks_pass") is True,
        "oracle_mpmath_dps_at_least_140": audit_meta.get("mpmath_dps", 0) >= 140,
        "oracle_digits_at_least_115": audit_meta.get(
            "formatted_significant_digits", 0
        )
        >= 115,
        "oracle_has_explicit_formatting_radius": audit_meta.get(
            "explicit_decimal_formatting_radius"
        )
        is True,
        "dual_oracle_all_pass": audit.get(
            "dual_against_independent_mpmath", {}
        ).get("all_pass")
        is True,
        "dual_oracle_failures_empty": audit.get(
            "dual_against_independent_mpmath", {}
        ).get("failures")
        == [],
        "symbolic_checks_all_true": all(
            audit.get("symbolic_identity_checks", {}).values()
        ),
        "proof_structure_checks_all_true": all(
            audit.get("proof_structure_checks", {}).values()
        ),
        "H1272_audit_pass": dependency.get("all_pass") is True,
        "H1272_range_is_r_ge_5_over_2": dependency.get("declared_range")
        == "r>=5/2",
        "H1272_live_hash_matches": h1272_path.is_file()
        and sha256(h1272_path) == dependency.get("sha256"),
        "engine_live_hash_matches": engine_path.is_file()
        and sha256(engine_path) == engine.get("sha256"),
        "runner_live_hash_matches": runner_path.is_file()
        and sha256(runner_path) == audit_meta.get("canonical_runner_sha256"),
        "audited_union_is_4_to_22": union.get("covered_r_interval")
        == ["4", "22"],
        "audited_union_gap_free": union.get("exact_gap_free_chain") is True,
        "audited_union_has_360_boxes": union.get("accepted_box_count") == 360,
        "all_three_mean_reports_audited": len(mean_segments) == 3
        and len(audited_reports) == 3,
        "detailed_report_hashes_match_audit": detailed_hashes_match,
    }
    return {
        "path": normalized(path),
        "sha256": sha256(path),
        "engine": {
            "path": engine.get("path"),
            "sha256": engine.get("sha256"),
        },
        "H1272": {
            "path": dependency.get("path"),
            "sha256": dependency.get("sha256"),
        },
        "oracle": {
            "mpmath_dps": audit_meta.get("mpmath_dps"),
            "formatted_significant_digits": audit_meta.get(
                "formatted_significant_digits"
            ),
            "explicit_decimal_formatting_radius": audit_meta.get(
                "explicit_decimal_formatting_radius"
            ),
            "comparison_count": audit.get(
                "dual_against_independent_mpmath", {}
            ).get("value_and_r_derivative_comparison_count"),
        },
        "checks": checks,
        "passes": all(checks.values()),
    }


def build_verified_report(args: argparse.Namespace) -> dict[str, Any]:
    report = implementation.build_report(args)

    # q''(r)=pi*exp(r)*(r+2)>0 for r>0, so the endpoint suffices.  This
    # avoids the severe wrapping seen in a single ball spanning [3.23,71].
    left = qarb(Fraction(323, 100))
    q_prime_left = arb.pi() * left.exp() * (left + 1) - arb(9) / 4
    report["saddle_checks"].pop("q_prime_on_323_over_100_to_71", None)
    report["saddle_checks"].update(
        {
            "q_prime_at_323_over_100": implementation.interval_payload(
                q_prime_left
            ),
            "q_second_derivative_sign": "pi*exp(r)*(r+2)>0 for r>0",
            "monotonicity_conclusion": (
                "q'(r)>=q'(323/100)>0 for every r>=323/100"
            ),
        }
    )
    report["checks"]["saddle_map_covers_t_domain"] = bool(
        q_prime_left.lower() > 0
    )

    independent = inspect_independent_audit(Path(args.audit), report)
    report["independent_mean_value_audit"] = independent
    report["checks"]["independent_mean_value_audit_pass"] = independent[
        "passes"
    ]
    report["canonical_generator"] = normalized(Path(__file__))
    report["all_checks_pass"] = all(report["checks"].values())
    report["classification"] = (
        "h1319_full_xi_compact_interval_closed"
        if report["all_checks_pass"]
        else "h1319_full_xi_compact_interval_manifest_failure"
    )
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--first", default=str(implementation.DEFAULT_FIRST))
    parser.add_argument("--mean-4-5", default=str(implementation.DEFAULT_MEAN_4_5))
    parser.add_argument("--mean-5-10", default=str(implementation.DEFAULT_MEAN_5_10))
    parser.add_argument("--mean-10-22", default=str(implementation.DEFAULT_MEAN_10_22))
    parser.add_argument("--high", default=str(implementation.DEFAULT_HIGH))
    parser.add_argument("--audit", default=str(DEFAULT_AUDIT))
    parser.add_argument("--out", default=str(implementation.DEFAULT_OUT))
    parser.add_argument("--markdown-out", default=str(implementation.DEFAULT_MARKDOWN))
    parser.add_argument("--dps", type=int, default=90)
    args = parser.parse_args()
    ctx.dps = args.dps

    report = build_verified_report(args)
    out = Path(args.out)
    markdown = Path(args.markdown_out)
    out.parent.mkdir(parents=True, exist_ok=True)
    markdown.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    note = implementation.build_markdown(report)
    note += (
        "\n## Independent audit gate\n\n"
        f"Canonical audit: `{report['independent_mean_value_audit']['path']}`. "
        f"PASS=`{report['independent_mean_value_audit']['passes']}`. The audit "
        "uses a 140-dps independent mpmath oracle, 115 formatted digits, and "
        "an explicit decimal-formatting radius; it also pins H1272 and the "
        "mean-value engine by SHA-256.\n"
    )
    markdown.write_text(note, encoding="utf-8")
    print(
        json.dumps(
            {
                "schema": report["schema"],
                "classification": report["classification"],
                "all_checks_pass": report["all_checks_pass"],
                "independent_audit_pass": report[
                    "independent_mean_value_audit"
                ]["passes"],
                "segments": len(report["r_proof_chain"]["segments"]),
                "out": str(out),
                "markdown": str(markdown),
            },
            indent=2,
        )
    )
    return 0 if report["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
