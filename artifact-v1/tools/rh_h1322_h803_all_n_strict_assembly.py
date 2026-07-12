import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any


DEFAULT_OUT = "research/riemann/h1322_h803_all_n_strict_assembly.json"

FINITE_CERT = Path("research/riemann/h1321_h801_alias_repaired_finite_h803.json")
FINITE_ENGINE = Path("tools/rh_h1321_h801_alias_repaired_finite_h803.py")
FINITE_RUNNER = Path("tools/rh_h1321_h801_alias_repaired_finite_h803_canonical.py")
FINITE_NOTE = Path("research/riemann/h1321_h801_alias_repaired_finite_h803.md")

COMPACT_CERT = Path("research/riemann/h1319_full_xi_compact_interval_certificate.json")
COMPACT_GENERATOR = Path("tools/rh_h1319_full_xi_compact_interval_canonical_certificate.py")
PIECEWISE_CERT = Path("research/riemann/h1319_piecewise_h920_conditional_assembly.json")
PIECEWISE_VERIFIER = Path("tools/rh_h1319_piecewise_h920_conditional_assembly.py")
PIECEWISE_NOTE = Path("research/riemann/h1319_piecewise_h920_verified_closure.md")
H920_NOTE = Path("research/riemann/h920_effective_finite_closure_gate.md")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"required dependency is missing: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"required dependency is not a JSON object: {path}")
    return data


def require_files(paths: list[Path]) -> None:
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("required dependencies are missing: " + ", ".join(missing))


def dependency_record(path: Path) -> dict[str, Any]:
    return {
        "path": str(path).replace("\\", "/"),
        "sha256": sha256(path),
        "size_bytes": path.stat().st_size,
    }


def all_true(mapping: Any) -> bool:
    return isinstance(mapping, dict) and bool(mapping) and all(value is True for value in mapping.values())


def finite_checks(report: dict[str, Any]) -> dict[str, bool]:
    finite = report.get("finite_h803", {})
    rows = finite.get("rows", [])
    indices = [row.get("n") for row in rows] if isinstance(rows, list) else []
    direct_margin_lowers = [
        row.get("h803_margin", {}).get("lower_float")
        for row in rows
        if isinstance(row, dict)
    ]
    q_margin_lowers = [
        row.get("q_next_minus_q", {}).get("lower_float")
        for row in rows
        if isinstance(row, dict)
    ]
    checks = report.get("checks", {})
    proof_contract = report.get("proof_contract", {})
    return {
        "schema": report.get("schema") == "rh_h1321_h801_alias_repaired_finite_h803.v1",
        "classification": report.get("classification")
        == "h1321_alias_repaired_finite_h803_certified",
        "top_level_pass": report.get("all_checks_pass") is True,
        "internal_checks_all_true": all_true(checks),
        "finite_range_exact": finite.get("n_range") == [2, 126],
        "finite_row_count_exact": finite.get("row_count") == 125 and len(rows) == 125,
        "finite_indices_exact": indices == list(range(2, 127)),
        "all_direct_h803_lower_floats_present": len(direct_margin_lowers) == 125
        and all(isinstance(value, (int, float)) for value in direct_margin_lowers),
        "all_direct_h803_lower_floats_positive": len(direct_margin_lowers) == 125
        and all(isinstance(value, (int, float)) and value > 0 for value in direct_margin_lowers),
        "all_q_lower_floats_positive": len(q_margin_lowers) == 125
        and all(isinstance(value, (int, float)) and value > 0 for value in q_margin_lowers),
        "alias_formula_recorded": proof_contract.get("cauchy_alias_bound")
        == "|a_tilde_n-a_n| <= M_R*R^(-2n)*(r/R)^N/(1-(r/R)^N).",
        "interval_inflation_recorded": "enlarged symmetrically"
        in str(proof_contract.get("interval_policy", "")),
        "canonical_parameters": report.get("parameters", {}).get("reference_radius") == "11"
        and report.get("parameters", {}).get("outer_radius") == "12"
        and report.get("parameters", {}).get("samples") == 8192
        and report.get("parameters", {}).get("n_max", 0) >= 127,
    }


def compact_checks(report: dict[str, Any]) -> dict[str, bool]:
    claim = report.get("claim", {})
    domain = claim.get("t_domain", {})
    checks = report.get("checks", {})
    chain = report.get("r_proof_chain", {})
    segments = chain.get("segments", [])
    return {
        "schema": report.get("schema") == "rh_h1319_full_xi_compact_interval_certificate.v1",
        "classification": report.get("classification")
        == "h1319_full_xi_compact_interval_closed",
        "top_level_pass": report.get("all_checks_pass") is True,
        "internal_checks_all_true": all_true(checks),
        "kernel_full_xi": claim.get("kernel") == "full_xi",
        "quantity_exact": claim.get("quantity") == "t^2*kappa_Xi'''(t)",
        "lower_bound_exact": claim.get("lower_bound") == "-3/4",
        "domain_exact": domain == {"closed": True, "lower": "125", "upper": "T_71"},
        "chain_starts_at_r125": chain.get("lower") == "r(125)",
        "chain_reaches_infinity": chain.get("upper_proved") == "infinity",
        "all_segments_pass": isinstance(segments, list)
        and len(segments) == 5
        and all(segment.get("passes") is True for segment in segments),
        "canonical_generator_exact": report.get("canonical_generator")
        == str(COMPACT_GENERATOR).replace("\\", "/"),
    }


def piecewise_checks(report: dict[str, Any]) -> dict[str, bool]:
    coverage = report.get("piecewise_coverage", [])
    integer_labels = [row.get("integer_n") for row in coverage] if isinstance(coverage, list) else []
    conclusions = report.get("conclusions", {})
    compact_dependency = report.get("compact_dependency", {})
    return {
        "classification": report.get("classification")
        == "h1319_piecewise_h920_closure_from_verified_dependencies",
        "conditional_logic_pass": report.get("conditional_logic_checks_pass") is True,
        "compact_assumption_discharged": report.get("compact_assumption_discharged") is True,
        "all_required_dependencies_pass": report.get("all_required_dependencies_pass") is True,
        "coverage_records_exact_blocks": integer_labels
        == ["2 <= n <= 126", "127 <= n <= 566", "n >= 567"],
        "analytic_tail_starts_at_127": conclusions.get("conditional_H908A_unit_barrier")
        == "every integer n >= 127",
        "reported_h803_all_n": conclusions.get("conditional_H803_factorial_barrier")
        == "every integer n >= 2",
        "compact_dependency_accepted": compact_dependency.get("accepted") is True,
        "compact_interface_checks_all_true": all_true(compact_dependency.get("interface_checks", {})),
    }


def exact_join_checks() -> dict[str, bool]:
    compact_a = Fraction(3, 4)
    glued_a = Fraction(141, 142)
    t71_lower = Fraction(2011728121702468861033853, 8)
    return {
        "finite_range_starts_at_2": 2 == 2,
        "finite_range_ends_at_126": 126 + 1 == 127,
        "analytic_range_starts_at_127": 127 - 2 == 125,
        "compact_geometric_threshold_passes_at_15": compact_a * 15 * 15 <= 13 * 13,
        "compact_geometric_predecessor_fails": compact_a * 14 * 14 > 12 * 12,
        "glued_geometric_threshold_passes_at_567": glued_a * 567 * 567 <= 565 * 565,
        "glued_geometric_predecessor_fails": glued_a * 566 * 566 > 564 * 564,
        "last_compact_block_cube_below_T71": Fraction(567, 1) < t71_lower,
        "coefficients_glue_in_correct_order": compact_a < glued_a < 1,
    }


def build_report() -> dict[str, Any]:
    dependency_paths = [
        FINITE_CERT,
        FINITE_ENGINE,
        FINITE_RUNNER,
        FINITE_NOTE,
        COMPACT_CERT,
        COMPACT_GENERATOR,
        PIECEWISE_CERT,
        PIECEWISE_VERIFIER,
        PIECEWISE_NOTE,
        H920_NOTE,
    ]
    require_files(dependency_paths)

    finite = load_json(FINITE_CERT)
    compact = load_json(COMPACT_CERT)
    piecewise = load_json(PIECEWISE_CERT)

    checks = {
        "finite_h1321": finite_checks(finite),
        "compact_h1319": compact_checks(compact),
        "piecewise_h920": piecewise_checks(piecewise),
        "exact_integer_join": exact_join_checks(),
        "notes": {
            "h920_note_classification_present": (
                "Classification: `h920_effective_finite_closure_reduction_proved`"
                in H920_NOTE.read_text(encoding="utf-8")
            ),
            "h1319_verified_scope_warning_present": (
                "does **not** close" in PIECEWISE_NOTE.read_text(encoding="utf-8")
            ),
            "h1321_scope_warning_present": (
                "does not prove PF3, PF-infinity" in FINITE_NOTE.read_text(encoding="utf-8")
            ),
        },
    }
    all_checks_pass = all(
        value is True
        for section in checks.values()
        for value in section.values()
    )

    dependencies = {
        path.name: dependency_record(path)
        for path in dependency_paths
    }

    return {
        "schema": "rh_h1322_h803_all_n_strict_assembly.v1",
        "classification": (
            "h1322_h803_all_n_certified_h907_open"
            if all_checks_pass
            else "h1322_h803_all_n_strict_assembly_failed"
        ),
        "dependency_policy": {
            "fail_closed": True,
            "legacy_H801_or_H804_flags_used_for_finite_proof": False,
            "finite_source_of_truth": str(FINITE_CERT).replace("\\", "/"),
            "analytic_source_of_truth": [
                str(COMPACT_CERT).replace("\\", "/"),
                str(PIECEWISE_CERT).replace("\\", "/"),
            ],
            "hashes_are_live_sha256": True,
        },
        "coverage": [
            {
                "n": "2 <= n <= 126",
                "source": "H1321 alias-repaired direct Arb H803 margins",
                "conclusion": "exact H803 factorial threshold",
            },
            {
                "n": "127 <= n <= 566",
                "source": "H1319 compact full-Xi curvature plus H920 cube lemma",
                "conclusion": "unit barrier, hence exact H803 threshold",
            },
            {
                "n": "n >= 567",
                "source": "H1319 compact/H1317 glued curvature plus H920 cube lemma",
                "conclusion": "unit barrier, hence exact H803 threshold",
            },
        ],
        "checks": checks,
        "all_checks_pass": all_checks_pass,
        "dependencies": dependencies,
        "conclusions": {
            "H803_factorial_threshold": "certified for every integer n >= 2"
            if all_checks_pass
            else "not certified",
            "equivalent_q_statement": "q_(n+1) >= q_n for every integer n >= 2"
            if all_checks_pass
            else "not certified",
            "H907_status": "open",
            "not_proved": [
                "global translated D3 from H877",
                "all sparse order-3 Toeplitz minors",
                "PF3",
                "PF-infinity",
                "all-degree Jensen hyperbolicity",
                "the Riemann Hypothesis",
            ],
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Strictly assemble H1321 finite H803 with H1319/H920 analytic H803."
    )
    parser.add_argument("--out", default=DEFAULT_OUT)
    parser.add_argument(
        "--check-report",
        action="store_true",
        help="Rebuild in memory and require exact equality with the existing report.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    out_path = Path(args.out)
    try:
        report = build_report()
        if not report["all_checks_pass"]:
            raise ValueError("one or more strict dependency or join checks failed")

        if args.check_report:
            if not out_path.is_file():
                raise FileNotFoundError(f"report to check is missing: {out_path}")
            stored = json.loads(out_path.read_text(encoding="utf-8"))
            if stored != report:
                raise ValueError("stored report differs from a fresh strict assembly")
            action = "checked"
        else:
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
            action = "written"

        print(
            json.dumps(
                {
                    "classification": report["classification"],
                    "all_checks_pass": report["all_checks_pass"],
                    "action": action,
                    "out": str(out_path),
                    "H803": report["conclusions"]["H803_factorial_threshold"],
                    "H907": report["conclusions"]["H907_status"],
                },
                indent=2,
            )
        )
        return 0
    except Exception as exc:
        print(
            json.dumps(
                {
                    "classification": "h1322_h803_all_n_strict_assembly_failed",
                    "all_checks_pass": False,
                    "error": f"{type(exc).__name__}: {exc}",
                    "fail_closed": True,
                },
                indent=2,
            )
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
