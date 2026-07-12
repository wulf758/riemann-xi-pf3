#!/usr/bin/env python3
"""High-precision, explicitly rounded runner for the H1319 audit."""

from __future__ import annotations

from typing import Any

import mpmath as mp
from flint import arb

import rh_h1319_score_mean_value_independent_audit as audit


def enclosed_mp_ball(value: mp.mpf) -> arb:
    """Enclose the decimal formatting error of a 115-digit oracle value."""

    if value == 0:
        return arb(0, 1) * arb("1e-130")
    significant_digits = 115
    text = mp.nstr(value, significant_digits, strip_zeros=False)
    decimal_exponent = int(mp.floor(mp.log10(abs(value))))
    formatting_radius = arb(f"1e{decimal_exponent-significant_digits+1}")
    return arb(text) + arb(0, 1) * formatting_radius


def dual_probe_checks_high_precision() -> dict[str, Any]:
    """Repeat every Dual probe at 140 dps with explicit oracle radii."""

    mp.mp.dps = 140
    probes = (
        (mp.mpf("4.025"), mp.mpf("-24")),
        (mp.mpf("4.025"), mp.mpf("0")),
        (mp.mpf("4.025"), mp.mpf("24")),
        (mp.mpf("10.025"), mp.mpf("-24")),
        (mp.mpf("10.025"), mp.mpf("5")),
        (mp.mpf("10.025"), mp.mpf("24")),
        (mp.mpf("21.975"), mp.mpf("-24")),
        (mp.mpf("21.975"), mp.mpf("0")),
        (mp.mpf("21.975"), mp.mpf("24")),
    )
    names = (
        "psi",
        "density",
        "dominant_score_excess",
        "finite_kernel",
        "finite_kernel_w_derivative",
        "integrand_Z",
        "integrand_M2",
        "integrand_SA",
        "integrand_SC",
        "scale_8t2_over_A32",
    )
    failures: list[dict[str, str]] = []
    comparisons = 0
    for r_value, w_value in probes:
        model = audit.engine.DualScoreModel(
            arb(mp.nstr(r_value, 40)), arb(1), 8, 48
        )
        balls = audit.arb_fields(model, arb(mp.nstr(w_value, 40)))
        exact = audit.exact_fields(r_value, w_value)
        exact_derivatives = [
            mp.diff(
                lambda variable, index=index: audit.exact_fields(
                    variable, w_value
                )[index],
                r_value,
            )
            for index in range(len(names))
        ]
        for name, ball, value, derivative in zip(
            names, balls, exact, exact_derivatives
        ):
            comparisons += 2
            if not ball.val.overlaps(enclosed_mp_ball(value)):
                failures.append(
                    {
                        "r": mp.nstr(r_value),
                        "w": mp.nstr(w_value),
                        "field": name,
                        "component": "value",
                    }
                )
            if not ball.der.overlaps(enclosed_mp_ball(derivative)):
                failures.append(
                    {
                        "r": mp.nstr(r_value),
                        "w": mp.nstr(w_value),
                        "field": name,
                        "component": "r_derivative",
                    }
                )
    return {
        "mpmath_dps": 140,
        "formatted_significant_digits": 115,
        "explicit_decimal_formatting_radius": True,
        "probe_count": len(probes),
        "field_count_per_probe": len(names),
        "comparison_count": comparisons,
        "failures": failures,
        "all_pass": not failures,
    }


audit.dual_probe_checks = dual_probe_checks_high_precision


if __name__ == "__main__":
    raise SystemExit(audit.main())
