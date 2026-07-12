#!/usr/bin/env python3
"""Exact arithmetic checks for the H1317 effective Xi transfer assembly."""

from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path

from flint import arb, ctx


def main() -> None:
    # Centered quotient bounds for R_5, and the independently proved m>=2 tail.
    r5 = (207, 1035, 7659)
    mtail = (1311, 6555, 48507)
    total = tuple(math.ceil(Fraction(a, 2) + b) for a, b in zip(r5, mtail))
    b1, b2, b3 = total
    b_xi = 2 * b3 + 12 * b1 * b2 + 16 * b1**3

    # Choose r>=71 and t>=142*B_Xi.  Then 70/r+B_Xi/t<=141/142.
    coefficient = Fraction(70, 71) + Fraction(1, 142)
    threshold_multiple = 142 * b_xi

    # The saddle threshold T_71=q(71)/2 is vastly larger than 142*B_Xi.
    # Use only pi>3 and e>2, so this comparison is exact integer arithmetic.
    t71_lower = (
        Fraction(3 * 71 * 2**71, 1)
        - Fraction(9 * 71, 4)
        - 1
    ) / 2

    # H920's geometric threshold for a=141/142 is below 567.
    geom_square_guard = Fraction(141, 142) < Fraction(565, 567) ** 2

    # Rigorous Arb check of the global sign estimate for the full correction.
    # sum_{m>=2} m^4 exp(-pi(m^2-1)) < 1/400 < 3/(2*pi).
    ctx.prec = 256
    pi = arb.pi()
    theta_tail = sum(
        arb(m) ** 4 * (-pi * (m * m - 1)).exp()
        for m in range(2, 20)
    )
    # The omitted m>=20 terms are below a deliberately enormous geometric cap.
    omitted_cap = arb(20) ** 4 * (-pi * (20**2 - 1)).exp() / (
        1 - (-4 * pi).exp()
    )
    sign_upper = theta_tail + omitted_cap
    c5 = arb(3) / (2 * pi)

    root = Path(__file__).resolve().parents[1]
    dependencies = {
        "moments": root / "research" / "riemann" / "h1313_h946_global_even_moment_bounds.md",
        "beta5": root / "research" / "riemann" / "h1316_beta5_direct_transfer_bound.md",
        "mtail": root / "research" / "riemann" / "h1316_xi_mtail_explicit_transfer.md",
        "audit": root / "research" / "riemann" / "h1316_xi_transfer_chain_independent_audit.md",
    }
    dependency_status = {name: path.exists() for name, path in dependencies.items()}

    checks = {
        "B1_is_1415": b1 == 1415,
        "B2_is_7073": b2 == 7073,
        "B3_is_52337": b3 == 52337,
        "B_Xi_is_45450578214": b_xi == 45_450_578_214,
        "coefficient_is_141_over_142": coefficient == Fraction(141, 142),
        "coefficient_lt_1": coefficient < 1,
        "T71_lower_exceeds_142_BXi": t71_lower > threshold_multiple,
        "H920_geometric_threshold_below_567": geom_square_guard,
        "theta_tail_lt_1_over_400": sign_upper.upper() < arb(1) / 400,
        "one_over_400_lt_c5": arb(1) / 400 < c5.lower(),
        "dependencies_exist": all(dependency_status.values()),
    }
    report = {
        "id": "h1317_effective_xi_transfer_assembly",
        "R5_derivative_constants": {"B1": r5[0], "B2": r5[1], "B3": r5[2]},
        "mtail_derivative_constants": {"B1": mtail[0], "B2": mtail[1], "B3": mtail[2]},
        "total_epsilon_constants": {"B1": b1, "B2": b2, "B3": b3},
        "B_Xi": b_xi,
        "full_bound": "kappa_Xi'''(t) >= -(70/r(t)+B_Xi/t)/t^2",
        "effective_choice": {
            "r_floor": 71,
            "t_floor": "T_71=(71*pi*exp(71)-(9/4)*71-1)/2",
            "coefficient": str(coefficient),
            "H920_geometric_threshold_upper": 567,
        },
        "exact_T71_lower": str(t71_lower),
        "required_142_BXi": threshold_multiple,
        "pointwise_sign_tail_upper": sign_upper.str(30),
        "c5": c5.str(30),
        "dependency_status": dependency_status,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "scope": "effective eventual Xi transfer; finite H803/Jensen/RH gaps remain",
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    if not report["all_checks_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
