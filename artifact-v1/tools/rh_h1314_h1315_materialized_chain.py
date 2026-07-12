#!/usr/bin/env python3
"""End-to-end reproducibility gate for the materialized H1314-H1315 chain."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
RIEMANN = ROOT / "research" / "riemann"
OUTPUT = RIEMANN / "h1314_h1315_materialized_chain.json"


COMPONENTS = [
    {
        "name": "H943C_U4",
        "script": ROOT / "tools" / "rh_h943c_u4_five_interval_certificate.py",
        "json": RIEMANN / "h943c_u4_five_interval_certificate.json",
    },
    {
        "name": "H1314_parameters",
        "script": ROOT / "tools" / "rh_h1314_h943c_parameter_polynomial_verifier.py",
        "json": RIEMANN / "h1314_h943c_parameter_polynomial_certificate.json",
    },
    {
        "name": "H1314_uniform_pair",
        "script": ROOT / "tools" / "rh_h1314_h943c_uniform_pair_verifier.py",
        "json": RIEMANN / "h1314_h943c_uniform_pair_certificate.json",
    },
    {
        "name": "H1315_moment_box",
        "script": ROOT / "tools" / "rh_h1315_gamma_log_moment_box_global_certificate.py",
        "json": RIEMANN / "h1315_gamma_log_moment_box_global_assembly.json",
    },
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_component(component: dict[str, Any]) -> dict[str, Any]:
    runs = []
    for run_number in (1, 2):
        completed = subprocess.run(
            [sys.executable, str(component["script"])],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        json_path: Path = component["json"]
        parsed = (
            json.loads(json_path.read_text(encoding="utf-8"))
            if completed.returncode == 0 and json_path.exists()
            else {}
        )
        runs.append(
            {
                "run": run_number,
                "returncode": completed.returncode,
                "json_exists": json_path.exists(),
                "json_sha256": sha256(json_path) if json_path.exists() else None,
                "certificate_id": parsed.get("id"),
                "all_checks_pass": bool(parsed.get("all_checks_pass")),
                "stderr_tail": completed.stderr[-1000:],
            }
        )

    deterministic = (
        runs[0]["json_sha256"] is not None
        and runs[0]["json_sha256"] == runs[1]["json_sha256"]
    )
    passed = deterministic and all(
        row["returncode"] == 0 and row["all_checks_pass"] for row in runs
    )
    return {
        "name": component["name"],
        "script": str(component["script"].relative_to(ROOT)).replace("\\", "/"),
        "json": str(component["json"].relative_to(ROOT)).replace("\\", "/"),
        "runs": runs,
        "deterministic_json": deterministic,
        "passed": passed,
    }


def load_component(name: str) -> dict[str, Any]:
    component = next(item for item in COMPONENTS if item["name"] == name)
    return json.loads(component["json"].read_text(encoding="utf-8"))


def cross_link_certificate() -> dict[str, Any]:
    parameters = load_component("H1314_parameters")
    uniform = load_component("H1314_uniform_pair")
    h1315 = load_component("H1315_moment_box")

    h1314_first = uniform["H943C_budgets"]["first_budget"]
    h1314_second = uniform["H943C_budgets"]["second_budget"]
    h1315_first = h1315["exact_arithmetic"]["first_score_budget"]["value"][
        "fraction"
    ]
    h1315_second = h1315["exact_arithmetic"]["cubic_score_budget"]["value"][
        "fraction"
    ]

    checks = {
        "parameter_json_consumed": bool(parameters["all_checks_pass"]),
        "uniform_pair_json_consumes_parameter": bool(
            uniform["parameter_certificate"]["all_checks_pass"]
        ),
        "uniform_pair_statement": (
            uniform["statement"] == "R_r(t)<=51/100 for r>=5/2 and t>=0"
        ),
        "h1314_first_budget_is_663_over_1250": h1314_first == "663/1250",
        "h1314_second_budget_is_14229_over_5000": (
            h1314_second == "14229/5000"
        ),
        "h1315_first_budget_matches_h1314_json": h1315_first == h1314_first,
        "h1315_second_budget_matches_h1314_json": h1315_second == h1314_second,
        "h1315_coverage_passes": bool(h1315["coverage"]["passed"]),
        "h1315_moment_box_passes": bool(h1315["moment_box_assembly"]["passed"]),
        "h1315_final_constant_is_35": (
            h1315["cumulant_consequence"]["final_constant"]["fraction"] == "35"
        ),
    }
    return {
        "H1314_budget_values": [h1314_first, h1314_second],
        "H1315_budget_values": [h1315_first, h1315_second],
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def downstream_transfer_smoke() -> dict[str, Any]:
    script = ROOT / "tools" / "rh_h1317_effective_xi_transfer_assembly.py"
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
    passed = completed.returncode == 0 and bool(parsed.get("all_checks_pass"))
    return {
        "script": str(script.relative_to(ROOT)).replace("\\", "/"),
        "returncode": completed.returncode,
        "certificate_id": parsed.get("id"),
        "all_checks_pass": bool(parsed.get("all_checks_pass")),
        "passed": passed,
        "scope": "downstream smoke test, not a replacement for analytic dependencies",
    }


def main() -> None:
    component_runs = [run_component(component) for component in COMPONENTS]
    cross_link = cross_link_certificate()
    downstream = downstream_transfer_smoke()
    all_checks = (
        all(component["passed"] for component in component_runs)
        and cross_link["all_checks_pass"]
        and downstream["passed"]
    )
    report = {
        "id": "h1314_h1315_materialized_chain",
        "components": component_runs,
        "H1314_to_H1315_cross_link": cross_link,
        "downstream_H1317_smoke": downstream,
        "analytic_dependency_boundary": (
            "the scripts regenerate polynomial, rational, Arb and assembly "
            "calculations; proved analytic pairing identities remain cited inputs"
        ),
        "redundant_unmaterialized_branch": (
            "the stronger H946 R<=1/2+C_r*t^2 polynomial is absent from the "
            "cards and is not used by this H1314-H1315-Xi transfer chain"
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
