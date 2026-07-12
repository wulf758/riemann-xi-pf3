#!/usr/bin/env python3
"""Independent stdlib/Fraction audit of the H1324 transfer and H790 canary."""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PRIMARY_PATH = (
    ROOT
    / "research"
    / "riemann"
    / "h1324_moment_jensen_laguerre_transfer.json"
)
AUDIT_PATH = (
    ROOT
    / "research"
    / "riemann"
    / "h1324_moment_jensen_laguerre_transfer_independent_audit.json"
)


def text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def half_gamma_without_sqrt_pi(k: int) -> Fraction:
    """Return Gamma(k+1/2)/sqrt(pi) exactly."""

    return Fraction(math.factorial(2 * k), 4**k * math.factorial(k))


def coefficient_from_laguerre_transfer(d: int, n: int, j: int) -> Fraction:
    """Coefficient multiplying M_(n+j), independently of SymPy."""

    # K_(d,n), with sqrt(pi) cancelled against Gamma(n+d+1/2).
    K = Fraction(
        math.factorial(d),
        4**n,
    ) / half_gamma_without_sqrt_pi(n + d)

    # Coefficient of (u^2 X)^j in L_d^(n-1/2)(-u^2 X/4).
    laguerre = (
        half_gamma_without_sqrt_pi(n + d)
        / half_gamma_without_sqrt_pi(n + j)
        / math.factorial(d - j)
        / math.factorial(j)
        / 4**j
    )
    return K * laguerre


def coefficient_from_jensen(d: int, n: int, j: int) -> Fraction:
    m = n + j
    return Fraction(math.comb(d, j) * math.factorial(m), math.factorial(2 * m))


def build_audit() -> dict[str, object]:
    coefficient_checks: list[bool] = []
    checked = 0
    for d in range(13):
        for n in range(13):
            for j in range(d + 1):
                coefficient_checks.append(
                    coefficient_from_laguerre_transfer(d, n, j)
                    == coefficient_from_jensen(d, n, j)
                )
                checked += 1

    weights = (Fraction(9, 10), Fraction(1, 10))
    x_atoms = (Fraction(1), Fraction(16))
    moments = [sum((weights[k] * x_atoms[k] ** m for k in range(2)), Fraction()) for m in range(3)]
    gamma = [Fraction(math.factorial(m), math.factorial(2 * m)) * moments[m] for m in range(3)]
    jensen = [gamma[0], 2 * gamma[1], gamma[2]]
    discriminant = jensen[1] ** 2 - 4 * jensen[0] * jensen[2]

    # For d=2,n=0, K*L_2^(-1/2)(-xX/4)=1+xX+x^2X^2/12.
    integral_jensen = [
        moments[0],
        moments[1],
        moments[2] / 12,
    ]

    primary = json.loads(PRIMARY_PATH.read_text(encoding="utf-8"))
    primary_checks = {
        "classification": primary["classification"]
        == "exact_transfer_proved_generic_mixture_shortcut_invalidated",
        "primary_all_checks_pass": primary["all_checks_pass"] is True,
        "transfer_constant": primary["transfer_theorem"]["constant"]
        == "K_(d,n) = sqrt(pi)*d!/(4^n*Gamma(n+d+1/2))",
        "transfer_identity": primary["transfer_theorem"]["identity"]
        == (
            "J_gamma^(d,n)(X) = K_(d,n) integral u^(2n) "
            "L_d^(n-1/2)(-u^2*X/4) dmu(u)"
        ),
        "h790_coordinate_dictionary": primary["h790_canary"]["x_measure"]
        == "nu=(9/10)delta_1+(1/10)delta_16 in x=u^2"
        and primary["h790_canary"]["equivalent_u_measure"]
        == "mu=(9/10)delta_1+(1/10)delta_4 in u",
        "h790_discriminant": primary["h790_canary"]["discriminant"] == "-31/12",
        "remaining_gap_is_open": primary["remaining_gap"]["status"] == "open_not_claimed",
    }

    checks = {
        "fraction_coefficient_grid": all(coefficient_checks),
        "h790_moments": moments == [Fraction(1), Fraction(5, 2), Fraction(53, 2)],
        "h790_gamma": gamma == [Fraction(1), Fraction(5, 4), Fraction(53, 24)],
        "h790_jensen_coefficients": jensen
        == [Fraction(1), Fraction(5, 2), Fraction(53, 24)],
        "h790_integral_coefficients": integral_jensen == jensen,
        "h790_discriminant": discriminant == Fraction(-31, 12),
        "primary_fields": all(primary_checks.values()),
    }

    return {
        "schema": "rh_h1324_moment_jensen_laguerre_transfer_independent_audit.v0",
        "classification": "independent_fraction_audit_pass",
        "implementation": "python_stdlib_fraction_no_sympy_no_primary_import",
        "coefficient_grid": {
            "d_min": 0,
            "d_max": 12,
            "n_min": 0,
            "n_max": 12,
            "checks": checked,
        },
        "h790_reconstruction": {
            "x_atoms": [text(value) for value in x_atoms],
            "u_atoms": ["1", "4"],
            "weights": [text(value) for value in weights],
            "moments_M0_M1_M2": [text(value) for value in moments],
            "gamma_0_1_2": [text(value) for value in gamma],
            "jensen_ascending_coefficients": [text(value) for value in jensen],
            "discriminant": text(discriminant),
        },
        "primary_field_checks": primary_checks,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-report", action="store_true")
    args = parser.parse_args()
    audit = build_audit()
    if args.check_report:
        pinned = json.loads(AUDIT_PATH.read_text(encoding="utf-8"))
        if pinned != audit:
            raise SystemExit(f"audit mismatch: {AUDIT_PATH}")
    print(json.dumps(audit, indent=2, sort_keys=True))
    return 0 if audit["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
