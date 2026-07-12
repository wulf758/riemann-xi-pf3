#!/usr/bin/env python3
"""Audit the sign inputs used by H1325's high-band upper cumulant bound."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "research" / "riemann"
SOURCES = {
    "H925": R / "h925_score_excess_budget_reduction.md",
    "H927": R / "h927_gamma_log_a_le_one_low_threshold.md",
    "H1315": R / "h1315_gamma_log_moment_box_global_assembly.md",
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def main() -> int:
    texts = {name: path.read_text(encoding="utf-8") for name, path in SOURCES.items()}
    checks = {
        "H927_a_strictly_positive": "and hence `a>0`." in texts["H927"],
        "H925_score_excess_nonnegative": "N(w)>=0." in texts["H925"],
        "H1315_X_identity": "E[N]/a = E[W^2 K]" in texts["H1315"],
        "H1315_Y_identity": "E[(W^2+2)N]/a" in texts["H1315"],
        "H1315_X_upper": "663/1250" in texts["H1315"],
        "H1315_E2_upper": "E[W^2] <= 26/25" in texts["H1315"],
    }
    passed = all(checks.values())
    report = {
        "schema": "rh_h1325_high_band_sign_dependency_audit.v0",
        "classification": (
            "h1325_high_band_sign_inputs_explicitly_audited"
            if passed
            else "h1325_high_band_sign_input_failure"
        ),
        "all_checks_pass": passed,
        "checks": checks,
        "deduction": (
            "a>0 and N(w)>=0 imply X=E[N]/a>=0 and "
            "Y=E[(W^2+2)N]/a>=0; hence -Y and -2*a^2*X^3 "
            "may be discarded in an upper bound"
        ),
        "source_hashes": {
            str(path.relative_to(ROOT)).replace("\\", "/"): digest(path)
            for path in SOURCES.values()
        },
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
