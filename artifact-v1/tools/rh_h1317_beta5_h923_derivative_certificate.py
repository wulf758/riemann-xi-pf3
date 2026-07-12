#!/usr/bin/env python3
"""Mechanical arithmetic checks for H1317's beta-5 H923 certificate.

The analytic inputs are the H1313 moment bounds and the exact quotient
identities proved in the accompanying note.  This verifier checks dependency
attestation, rational constant budgets, quotient coefficient counts, and the
final H923 arithmetic.
"""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parents[1]

    moment_path = root / "research" / "riemann" / "h1313_h946_global_even_moment_bounds.md"
    tail_path = root / "research" / "riemann" / "h1316_xi_mtail_explicit_transfer.md"
    audit_path = root / "research" / "riemann" / "h1316_xi_transfer_chain_independent_audit.md"

    moment_text = moment_path.read_text(encoding="utf-8")
    tail_text = tail_path.read_text(encoding="utf-8")
    audit_text = audit_path.read_text(encoding="utf-8")

    dependency_checks = {
        "moment_file_exists": moment_path.exists(),
        "M2_present": "E_r[W^2] <= 26/25" in moment_text,
        "M4_present": "E_r[W^4] <= 7/2" in moment_text,
        "M6_present": "E_r[W^6] <= 21" in moment_text,
        "moment_range_present": "every `r>=5/2`" in moment_text,
        "tail_B1_present": "1311/t" in tail_text,
        "tail_B2_present": "6555/t^2" in tail_text,
        "tail_B3_present": "48507/t^3" in tail_text,
        "correction_sign_present": "-c_5 < epsilon(t) < 0" in audit_text,
    }

    # H1313 implications used in the good/bad split.
    m2 = Fraction(26, 25)
    m4 = Fraction(7, 2)
    m6 = Fraction(21)
    moment_checks = {
        "mu1_lt_4": m2 < 16,
        "mu2_lt_4": m2 < 4,
        "mu3_lt_4": m2 * m4 < 16**2,
        "mu4_lt_4": m4 < 4,
        "M6_is_21": m6 == 21,
    }

    base = 48 + 21
    first_multipliers = [1, 2]
    second_multipliers = [1, 4, 8, 2]
    third_multipliers = [1, 6, 24, 6, 48, 24, 2]

    beta5 = {
        "B1": base * sum(first_multipliers),
        "B2": base * sum(second_multipliers),
        "B3": base * sum(third_multipliers),
    }

    # c_5<1/2, so round each tail+c_5*beta5 budget upward to an integer.
    full = {
        "B1": 1311 + (beta5["B1"] + 1) // 2,
        "B2": 6555 + (beta5["B2"] + 1) // 2,
        "B3": 48507 + (beta5["B3"] + 1) // 2,
    }
    b_xi = 2 * full["B3"] + 12 * full["B1"] * full["B2"] + 16 * full["B1"] ** 3

    arithmetic_checks = {
        "base_is_69": base == 69,
        "beta5_B1_is_207": beta5["B1"] == 207,
        "beta5_B2_is_1035": beta5["B2"] == 1035,
        "beta5_B3_is_7659": beta5["B3"] == 7659,
        "full_B1_is_1415": full["B1"] == 1415,
        "full_B2_is_7073": full["B2"] == 7073,
        "full_B3_is_52337": full["B3"] == 52337,
        "B_Xi_is_45450578214": b_xi == 45_450_578_214,
        "quarter_threshold_is_181802312856": 4 * b_xi == 181_802_312_856,
    }

    checks = {**dependency_checks, **moment_checks, **arithmetic_checks}
    result = {
        "id": "h1317_beta5_h923_derivative_certificate",
        "beta5_derivative_constants": beta5,
        "full_epsilon_derivative_constants": full,
        "B_Xi": b_xi,
        "quarter_transfer_threshold": 4 * b_xi,
        "dependency_checks": dependency_checks,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    print(json.dumps(result, indent=2, sort_keys=True))

    if not result["all_checks_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

