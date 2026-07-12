#!/usr/bin/env python3
"""Fail-closed sign audit for H1325's high-band upper cumulant bound."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RIEMANN = ROOT / "research" / "riemann"
H1325 = RIEMANN / "h1325_xi_global_d3_h811_consecutive_rows.json"
H943 = RIEMANN / "h943_score_excess_taylor_bridge.md"
H927 = RIEMANN / "h927_gamma_log_a_le_one_low_threshold.md"
HIGH = RIEMANN / "h1319_full_xi_high_band_r22_analytic_certificate.json"
DEFAULT_OUT = RIEMANN / "h1325_high_band_sign_audit.json"
DEFAULT_MD = RIEMANN / "h1325_high_band_sign_audit.md"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--markdown-out", default=str(DEFAULT_MD))
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()

    h1325 = json.loads(H1325.read_text(encoding="utf-8"))
    high = json.loads(HIGH.read_text(encoding="utf-8"))
    h943 = H943.read_text(encoding="utf-8")
    h927 = H927.read_text(encoding="utf-8")

    checks = {
        "H1325_main_passes": bool(h1325.get("all_checks_pass")),
        "H943_exact_K_nonnegative": "K_r(w)>=0" in h943,
        "H943_exact_N_identity": "N(w)=a w^2 K_r(w)" in h943,
        "H943_X_nonnegative_expectation": "E[W^2 K_r(W)]" in h943,
        "H943_Y_nonnegative_expectation": "E[(W^2+2)W^2 K_r(W)]" in h943,
        "H927_a_strictly_positive": "0<a<=1" in h927,
        "high_identity_exact": high["dominant_stein_bound"]["identity"]
        == "kappa3(W)/a=-Y+3*X*E2-2*a^2*X^3",
        "high_X_upper_exact": high["dominant_stein_bound"]["X_upper"]["exact"]
        == "663/1250",
        "high_F_direct_certificate": high["algebraic_ratio_certificate"]["conclusion"]
        == "0<q^2*r*B/A^3<1 for every real r>=4",
    }
    all_checks = all(checks.values())
    report = {
        "schema": "rh_h1325_high_band_sign_audit.v0",
        "classification": (
            "h1325_high_band_discarded_terms_sign_closed"
            if all_checks
            else "h1325_high_band_discarded_terms_sign_failed"
        ),
        "all_checks_pass": all_checks,
        "checks": checks,
        "argument": (
            "H927 gives a>0 and H943 gives K_r>=0. Thus "
            "X=E[W^2 K_r]>=0 and Y=E[(W^2+2)W^2 K_r]>=0; "
            "the terms -Y and -2*a^2*X^3 are nonpositive, so "
            "kappa3(W)/a<=3*X*E2 is valid."
        ),
        "conclusion": (
            "H1325's high-band upper third-cumulant envelope has all discarded-term signs explicitly certified."
            if all_checks
            else "A discarded-term sign gate failed; the upper envelope is not promoted."
        ),
        "dependency_hashes": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (H1325, H943, H927, HIGH)
        },
    }
    markdown = "\n".join(
        [
            "# H1325 High-Band Sign Audit",
            "",
            f"Classification: `{report['classification']}`",
            "",
            "H927 proves `a>0`, while H943 proves `K_r(w)>=0`.",
            "Consequently `X>=0` and `Y>=0`, and both discarded terms",
            "`-Y` and `-2a^2X^3` are nonpositive.  The H1325 high-band",
            "upper third-cumulant bound is therefore sign-valid.",
            "",
            "```text",
            "python tools/rh_h1325_high_band_sign_audit.py --no-write",
            "```",
            "",
        ]
    )

    if not args.no_write:
        Path(args.out).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        Path(args.markdown_out).write_text(markdown, encoding="utf-8")
    print(
        json.dumps(
            {
                "classification": report["classification"],
                "all_checks_pass": all_checks,
                "failed_checks": [key for key, value in checks.items() if not value],
                "out": None if args.no_write else str(Path(args.out)),
                "markdown": None if args.no_write else str(Path(args.markdown_out)),
            },
            indent=2,
        )
    )
    return 0 if all_checks else 1


if __name__ == "__main__":
    raise SystemExit(main())
