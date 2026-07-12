#!/usr/bin/env python3
"""Canonical end-to-end gate for H921/H941/H1314/H1315 and H1317."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
RIEMANN = ROOT / "research" / "riemann"
OUTPUT = RIEMANN / "h1314_h1315_end_to_end_certificate.json"


COMPONENTS = [
    {
        "name": "H921_main_coefficient",
        "script": "tools/rh_h921_explicit_saddle_main_coefficient_certificate.py",
        "json": "research/riemann/h921_explicit_saddle_main_coefficient_certificate.json",
        "id": "h921_explicit_saddle_main_coefficient_certificate",
    },
    {
        "name": "H941_low_band",
        "script": "tools/rh_h941_gamma_log_score_cover_assembly.py",
        "json": "research/riemann/h941_gamma_log_score_cover_assembly.json",
        "id": "h941_gamma_log_score_cover_assembly",
    },
    {
        "name": "H943C_U4",
        "script": "tools/rh_h943c_u4_five_interval_certificate.py",
        "json": "research/riemann/h943c_u4_five_interval_certificate.json",
        "id": "h943c_u4_five_interval_certificate",
    },
    {
        "name": "H1314_parameters",
        "script": "tools/rh_h1314_h943c_parameter_polynomial_verifier.py",
        "json": "research/riemann/h1314_h943c_parameter_polynomial_certificate.json",
        "id": "h1314_h943c_parameter_polynomial_certificate",
    },
    {
        "name": "H1314_pairing_tail_bridge",
        "script": "tools/rh_h1314_h943c_pairing_tail_bridge.py",
        "json": "research/riemann/h1314_h943c_pairing_tail_bridge.json",
        "id": "h1314_h943c_pairing_tail_bridge",
    },
    {
        "name": "H1314_uniform_pair",
        "script": "tools/rh_h1314_h943c_uniform_pair_verifier.py",
        "json": "research/riemann/h1314_h943c_uniform_pair_certificate.json",
        "id": "h1314_h943c_uniform_pair_certificate",
    },
    {
        "name": "H1315_moment_box",
        "script": "tools/rh_h1315_gamma_log_moment_box_global_certificate.py",
        "json": "research/riemann/h1315_gamma_log_moment_box_global_assembly.json",
        "id": "h1315_gamma_log_moment_box_global_assembly",
    },
]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def run_json_component(component: dict[str, str]) -> dict[str, Any]:
    script = ROOT / component["script"]
    output = ROOT / component["json"]
    rows = []
    for run_number in (1, 2):
        completed = subprocess.run(
            [sys.executable, str(script)],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        parsed = (
            json.loads(output.read_text(encoding="utf-8"))
            if completed.returncode == 0 and output.exists()
            else {}
        )
        rows.append(
            {
                "run": run_number,
                "returncode": completed.returncode,
                "json_exists": output.exists(),
                "json_sha256": sha256_file(output) if output.exists() else None,
                "certificate_id": parsed.get("id"),
                "expected_id": component["id"],
                "id_matches": parsed.get("id") == component["id"],
                "all_checks_pass": bool(parsed.get("all_checks_pass")),
                "stderr_tail": completed.stderr[-1000:],
            }
        )
    deterministic = (
        rows[0]["json_sha256"] is not None
        and rows[0]["json_sha256"] == rows[1]["json_sha256"]
    )
    passed = deterministic and all(
        row["returncode"] == 0
        and row["id_matches"]
        and row["all_checks_pass"]
        for row in rows
    )
    return {
        "name": component["name"],
        "script": component["script"],
        "json": component["json"],
        "runs": rows,
        "deterministic_json": deterministic,
        "passed": passed,
    }


def load(name: str) -> dict[str, Any]:
    component = next(row for row in COMPONENTS if row["name"] == name)
    return json.loads((ROOT / component["json"]).read_text(encoding="utf-8"))


def dependency_cross_links() -> dict[str, Any]:
    h921 = load("H921_main_coefficient")
    h941 = load("H941_low_band")
    u4 = load("H943C_U4")
    parameters = load("H1314_parameters")
    bridge = load("H1314_pairing_tail_bridge")
    uniform = load("H1314_uniform_pair")
    h1315 = load("H1315_moment_box")

    parameter_path = "research/riemann/h1314_h943c_parameter_polynomial_certificate.json"
    parameter_hash = sha256_file(ROOT / parameter_path)
    h1314_first = uniform["H943C_budgets"]["first_budget"]
    h1314_second = uniform["H943C_budgets"]["second_budget"]
    h1315_first = h1315["exact_arithmetic"]["first_score_budget"]["value"][
        "fraction"
    ]
    h1315_second = h1315["exact_arithmetic"]["cubic_score_budget"]["value"][
        "fraction"
    ]

    coefficient_7 = h1315["cumulant_consequence"]["coefficient_factor"][
        "numerator"
    ]
    centered_5 = h1315["cumulant_consequence"]["centered_moment_factor"][
        "numerator"
    ]
    final_35 = h1315["cumulant_consequence"]["final_constant"]["numerator"]

    checks = {
        "H921_range": h921["range"] == "alpha=9/4, r>=2",
        "H921_conclusion": h921["conclusion"] == "0<B/A^3<=7/(q^2*r)",
        "H921_factor_matches_H1315": coefficient_7 == 7,
        "H941_exact_domain": h941["statement"]["domain"] == "2<=r<=5/2",
        "H941_coverage_passes": bool(h941["coverage"]["passed"]),
        "H941_artifact_count_23": len(h941["artifacts"]) == 23,
        "H941_box_count_724": h941["coverage"]["box_count"] == 724,
        "H1315_low_band_matches_H941": (
            h1315["coverage"]["bands"][0]["left"] == "2"
            and h1315["coverage"]["bands"][0]["right"] == "5/2"
            and h1315["dependencies"]["H941"]["passed"]
        ),
        "bridge_parameter_path": bridge["parameter_dependency"]["path"] == parameter_path,
        "bridge_parameter_hash": bridge["parameter_dependency"]["sha256"] == parameter_hash,
        "bridge_pairing_residual_zero": (
            bridge["pairing_algebra"]["symbolic_residual"] == "0"
            and bridge["pairing_algebra"]["even_expectation_symbolic_residual"] == "0"
        ),
        "uniform_parameter_id": (
            uniform["parameter_certificate"]["id"] == parameters["id"]
        ),
        "uniform_parameter_path": (
            uniform["parameter_certificate"]["canonical_json"] == parameter_path
        ),
        "uniform_U4_id": uniform["compact_certificate"]["u4_certificate_id"] == u4["id"],
        "uniform_tail_matches_bridge_H1": (
            uniform["tail_certificate"]["junction"]["H1_upper"]
            == bridge["tail_constants"]["H1_upper"]
            == "-17/4"
        ),
        "uniform_tail_matches_bridge_derivative": (
            bridge["tail_constants"]["Hprime_negative_margin"] == "11/4"
        ),
        "H1314_first_budget": h1314_first == "663/1250",
        "H1314_second_budget": h1314_second == "14229/5000",
        "H1315_first_budget_matches_H1314": h1315_first == h1314_first,
        "H1315_second_budget_matches_H1314": h1315_second == h1314_second,
        "H1315_high_band_uses_H1314": (
            h1315["coverage"]["bands"][1]["left"] == "5/2"
            and h1315["coverage"]["bands"][1]["right"] == "infinity"
            and h1315["dependencies"]["H1314"]["passed"]
        ),
        "H1315_centered_factor": centered_5 == 5,
        "H1315_final_factor": final_35 == centered_5 * coefficient_7 == 35,
        "H1315_global_coverage": bool(h1315["coverage"]["passed"]),
        "H1315_moment_box": bool(h1315["moment_box_assembly"]["passed"]),
    }
    return {
        "H1314_budgets": [h1314_first, h1314_second],
        "H1315_budgets": [h1315_first, h1315_second],
        "H1315_factors": {
            "centered_moment": centered_5,
            "main_coefficient": coefficient_7,
            "gamma_log_constant": final_35,
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def h1317_link(final_35: int) -> dict[str, Any]:
    script = ROOT / "tools" / "rh_h1317_effective_xi_transfer_assembly.py"
    rows = []
    for run_number in (1, 2):
        completed = subprocess.run(
            [sys.executable, str(script)],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        try:
            parsed = json.loads(completed.stdout)
        except json.JSONDecodeError:
            parsed = {}
        rows.append(
            {
                "run": run_number,
                "returncode": completed.returncode,
                "stdout_sha256": sha256_bytes(completed.stdout.encode("utf-8")),
                "certificate_id": parsed.get("id"),
                "all_checks_pass": bool(parsed.get("all_checks_pass")),
                "full_bound": parsed.get("full_bound"),
                "B_Xi": parsed.get("B_Xi"),
            }
        )

    expected_id = "h1317_effective_xi_transfer_assembly"
    expected_bound = "kappa_Xi'''(t) >= -(70/r(t)+B_Xi/t)/t^2"
    converted_constant = 8 * final_35 // 4
    checks = {
        "both_runs_succeed": all(row["returncode"] == 0 for row in rows),
        "both_ids_match": all(row["certificate_id"] == expected_id for row in rows),
        "both_certificates_pass": all(row["all_checks_pass"] for row in rows),
        "deterministic_stdout": rows[0]["stdout_sha256"] == rows[1]["stdout_sha256"],
        "q_2t_third_derivative_conversion": converted_constant == 70,
        "full_bound_consumes_70": all(row["full_bound"] == expected_bound for row in rows),
        "B_Xi_matches": all(row["B_Xi"] == 45_450_578_214 for row in rows),
    }
    return {
        "script": str(script.relative_to(ROOT)).replace("\\", "/"),
        "conversion": "8*35/(2^2)=70 because q=2t and d_t^3=8*d_q^3",
        "converted_gamma_log_constant": converted_constant,
        "runs": rows,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    components = [run_json_component(component) for component in COMPONENTS]
    links = dependency_cross_links()
    final_35 = links["H1315_factors"]["gamma_log_constant"]
    transfer = h1317_link(final_35)
    all_checks = (
        all(component["passed"] for component in components)
        and links["all_checks_pass"]
        and transfer["all_checks_pass"]
    )
    report = {
        "id": "h1314_h1315_end_to_end_certificate",
        "components": components,
        "dependency_cross_links": links,
        "H1315_to_H1317_link": transfer,
        "analytic_dependency_boundary": (
            "polynomial, interval, rational, artifact-coverage and constant "
            "assemblies are regenerated; cited analytic identities are persisted "
            "in the H1314 pairing bridge and their source notes"
        ),
        "redundant_open_reproducibility_item": (
            "the stronger H946 R<=1/2+C_r*t^2 polynomial is absent from AXMEM; "
            "it is not used by this H1314-H1315-H1317 route"
        ),
        "all_checks_pass": all_checks,
    }
    OUTPUT.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    if not all_checks:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
