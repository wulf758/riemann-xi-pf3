#!/usr/bin/env python3
"""Exact scalar and dependency checks for the H1316 full Xi assembly."""

from __future__ import annotations

from decimal import Decimal, getcontext
from fractions import Fraction as Q
from math import factorial
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def taylor_exp(x: int, degree: int) -> Q:
    return sum((Q(x) ** k / factorial(k) for k in range(degree + 1)), Q(0))


def main() -> None:
    beta = (207, 1035, 7659)
    tail = (1311, 6555, 48507)
    full_exact = tuple(Q(tail[j]) + Q(beta[j], 2) for j in range(3))
    full_ceil = (1415, 7073, 52337)
    assert full_exact == (Q(2829, 2), Q(14145, 2), Q(104673, 2))
    assert all(full_exact[j] < full_ceil[j] for j in range(3))

    b1, b2, b3 = full_ceil
    b_xi = 2 * b3 + 12 * b1 * b2 + 16 * b1**3
    assert b_xi == 45450578214

    # Denominator-sign scalar estimates.
    assert taylor_exp(9, 12) > 7000
    assert taylor_exp(12, 7) > 1000
    geometric_upper = Q(16, 7000) / (1 - Q(16, 1000))
    assert geometric_upper < Q(1, 400)
    assert Q(1, 400) < Q(3, 8)  # 3/(2*pi)>3/8 because pi<4.

    # r>94 when q>=exp(100): exp(6)>94*pi via pi<22/7.
    exp6_lower = taylor_exp(6, 7)
    assert exp6_lower == Q(2101, 7)
    assert exp6_lower > Q(94 * 22, 7)

    # Full correction uses less than the exact 1/188 coefficient slack.
    assert 2**100 > 376 * b_xi
    assert Q(35, 47) + Q(1, 188) == Q(3, 4)

    # H920's geometric term is below 15; T+2 is overwhelmingly larger.
    # 2/(1-sqrt(3)/2)<15 follows from sqrt(3)<26/15.
    assert Q(3) < Q(26, 15) ** 2
    assert 2**99 + 2 > 15

    dependencies = {
        "beta_note": ROOT / "research/riemann/h1316_beta5_epsilon_derivative_box.md",
        "beta_tool": ROOT / "tools/rh_h1316_beta5_epsilon_box_certificate.py",
        "tail_note": ROOT / "research/riemann/h1316_xi_mtail_explicit_transfer.md",
        "tail_tool": ROOT / "tools/rh_h1316_xi_mtail_explicit_transfer.py",
        "audit_note": ROOT / "research/riemann/h1316_xi_transfer_chain_independent_audit.md",
        "h922m": ROOT / "research/riemann/h922m_qbox_analytic_certificate.md",
        "h923": ROOT / "research/riemann/h923_xi_log_correction_transfer_gate.md",
    }
    assert all(path.exists() for path in dependencies.values())

    getcontext().prec = 50
    n_eff_scale = Decimal(100).exp() / 2 + 2

    print("beta5_constants=", beta)
    print("tail_constants=", tail)
    print("full_exact_before_ceiling=", tuple(str(x) for x in full_exact))
    print("full_constants=", full_ceil)
    print("B_Xi=", b_xi)
    print("exp6_taylor7=", exp6_lower)
    print("N_eff_scale=", n_eff_scale)
    print("dependencies=", {key: path.exists() for key, path in dependencies.items()})
    print("all_checks_pass=true")


if __name__ == "__main__":
    main()
