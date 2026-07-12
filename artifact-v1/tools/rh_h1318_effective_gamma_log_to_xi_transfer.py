#!/usr/bin/env python3
"""Mechanical checks for the H1318 effective Xi transfer assembly."""

from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    gamma_path = root / "research" / "riemann" / "h922m_qbox_analytic_certificate.md"
    beta5_path = root / "research" / "riemann" / "h1317_beta5_h923_derivative_certificate.md"
    mtail_path = root / "research" / "riemann" / "h1316_xi_mtail_explicit_transfer.md"
    audit_path = root / "research" / "riemann" / "h1316_xi_transfer_chain_independent_audit.md"

    gamma = gamma_path.read_text(encoding="utf-8")
    beta5 = beta5_path.read_text(encoding="utf-8")
    mtail = mtail_path.read_text(encoding="utf-8")
    audit = audit_path.read_text(encoding="utf-8")

    dependencies = {
        "gamma_threshold_present": "q >= exp(100)" in gamma,
        "gamma_bound_present": "K_alpha'''(q) >= -35/(q^2 r_q)" in gamma,
        "beta5_B1_present": "B_1=1415" in beta5,
        "beta5_B2_present": "B_2=7073" in beta5,
        "beta5_B3_present": "B_3=52337" in beta5,
        "beta5_BXi_present": "45,450,578,214" in beta5,
        "mtail_audited_result_present": "48507/t^3" in mtail,
        "normalization_audit_present": "M_n=8*4^n" in audit,
        "h918_invalidation_present": "H918 does **not** prove" in audit,
    }

    b_xi = 45_450_578_214
    transfer_threshold = 4 * b_xi
    qbox_t = math.exp(100.0) / 2.0
    coefficient = Fraction(70, 94) + Fraction(1, 4)
    geometric_threshold = 2.0 / (1.0 - math.sqrt(187.0 / 188.0))
    r0 = 100.0 - math.log(100.0) - math.log(math.pi)

    arithmetic = {
        "quarter_threshold_correct": transfer_threshold == 181_802_312_856,
        "qbox_dominates_transfer_threshold": qbox_t > transfer_threshold,
        "r0_gt_94": r0 > 94.0,
        "coefficient_is_187_over_188": coefficient == Fraction(187, 188),
        "coefficient_lt_1": coefficient < 1,
        "qbox_dominates_geometric_threshold": qbox_t + 2 > geometric_threshold,
    }

    checks = {**dependencies, **arithmetic}
    result = {
        "id": "h1318_effective_gamma_log_to_xi_transfer",
        "a": "187/188",
        "T": "exp(100)/2",
        "B_Xi": b_xi,
        "quarter_transfer_threshold": transfer_threshold,
        "qbox_t_float": qbox_t,
        "r0_float": r0,
        "geometric_threshold_float": geometric_threshold,
        "dependencies": dependencies,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    print(json.dumps(result, indent=2, sort_keys=True))

    if not result["all_checks_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

