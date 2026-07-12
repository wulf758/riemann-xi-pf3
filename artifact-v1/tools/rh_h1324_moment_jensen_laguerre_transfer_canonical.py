#!/usr/bin/env python3
"""Canonical SymPy certificate for the H1324 moment-to-Jensen transfer."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
REPORT_PATH = (
    ROOT
    / "research"
    / "riemann"
    / "h1324_moment_jensen_laguerre_transfer.json"
)


def _q(value: sp.Expr) -> str:
    return str(sp.factor(value))


def build_certificate() -> dict[str, object]:
    d, n, j = sp.symbols("d n j", integer=True, nonnegative=True)
    m = n + j

    # Coefficient of M_(n+j) after expanding the generalized Laguerre
    # polynomial, multiplying by K_(d,n), and integrating u^(2n+2j).
    raw_coefficient = (
        sp.sqrt(sp.pi)
        * sp.gamma(d + 1)
        / (
            4 ** (n + j)
            * sp.gamma(n + j + sp.Rational(1, 2))
            * sp.gamma(d - j + 1)
            * sp.gamma(j + 1)
        )
    )

    # Gamma duplication at m=n+j:
    # sqrt(pi)/(4^m Gamma(m+1/2)) = Gamma(m+1)/Gamma(2m+1).
    after_duplication = (
        sp.gamma(d + 1)
        * sp.gamma(m + 1)
        / (
            sp.gamma(d - j + 1)
            * sp.gamma(j + 1)
            * sp.gamma(2 * m + 1)
        )
    )
    target_coefficient = sp.binomial(d, j) * sp.gamma(m + 1) / sp.gamma(2 * m + 1)
    symbolic_ratio = sp.simplify(after_duplication / target_coefficient)

    # Check SymPy's actual assoc_laguerre convention on a nontrivial exact
    # grid.  Y denotes u^2, so u^(2n)Y^j integrates to M_(n+j).
    X, Y = sp.symbols("X Y")
    grid_checks: list[bool] = []
    grid_count = 0
    for d_value in range(0, 9):
        for n_value in range(0, 7):
            alpha = sp.Rational(2 * n_value - 1, 2)
            kernel = sp.assoc_laguerre(d_value, alpha, -Y * X / 4)
            K = (
                sp.sqrt(sp.pi)
                * sp.factorial(d_value)
                / (4**n_value * sp.gamma(n_value + d_value + sp.Rational(1, 2)))
            )
            for j_value in range(d_value + 1):
                actual = sp.simplify(K * sp.expand(kernel).coeff(X, j_value) / Y**j_value)
                expected = (
                    sp.binomial(d_value, j_value)
                    * sp.factorial(n_value + j_value)
                    / sp.factorial(2 * (n_value + j_value))
                )
                grid_checks.append(sp.simplify(actual - expected) == 0)
                grid_count += 1

    # H790 canary.  Its atoms 1 and 16 are in x=u^2.  Equivalently, the
    # measure in the present u-variable has atoms 1 and 4.
    weight_1 = sp.Rational(9, 10)
    weight_16 = sp.Rational(1, 10)
    x_atoms = (sp.Integer(1), sp.Integer(16))
    moments = [
        sp.simplify(weight_1 * x_atoms[0] ** k + weight_16 * x_atoms[1] ** k)
        for k in range(3)
    ]
    gamma = [sp.factorial(k) * moments[k] / sp.factorial(2 * k) for k in range(3)]
    canary_jensen = sp.expand(gamma[0] + 2 * gamma[1] * X + gamma[2] * X**2)
    canary_discriminant = sp.discriminant(canary_jensen, X)

    K_20 = sp.sqrt(sp.pi) * sp.factorial(2) / sp.gamma(sp.Rational(5, 2))
    atomic_kernel_20 = sp.assoc_laguerre(2, sp.Rational(-1, 2), -Y * X / 4)
    canary_integral = sp.expand(
        weight_1 * K_20 * atomic_kernel_20.subs(Y, x_atoms[0])
        + weight_16 * K_20 * atomic_kernel_20.subs(Y, x_atoms[1])
    )

    # A symbolic check of the atomic root scaling: if lambda is a positive
    # root of L_d^(alpha), X=-4 lambda/u^2 maps the Laguerre argument to lambda.
    lam, u = sp.symbols("lambda u", positive=True)
    root_argument = sp.simplify(-u**2 * (-4 * lam / u**2) / 4)

    checks = {
        "symbolic_duplication_ratio_is_one": symbolic_ratio == 1,
        "sympy_assoc_laguerre_grid": all(grid_checks),
        "h790_moments_are_1_5over2_53over2": moments
        == [sp.Integer(1), sp.Rational(5, 2), sp.Rational(53, 2)],
        "h790_jensen_polynomial": canary_jensen
        == 1 + sp.Rational(5, 2) * X + sp.Rational(53, 24) * X**2,
        "h790_integral_equals_jensen": sp.simplify(canary_integral - canary_jensen) == 0,
        "h790_discriminant_is_minus_31_over_12": canary_discriminant
        == sp.Rational(-31, 12),
        "atomic_root_argument_maps_to_lambda": root_argument == lam,
    }

    return {
        "schema": "rh_h1324_moment_jensen_laguerre_transfer.v0",
        "classification": "exact_transfer_proved_generic_mixture_shortcut_invalidated",
        "definitions": {
            "moments": "M_m = integral u^(2m) dmu(u)",
            "ordinary_coefficients": "a_m = M_m/(2m)!",
            "jensen_coefficients": "gamma_m = m!*a_m = m!*M_m/(2m)!",
            "jensen_polynomial": "J_gamma^(d,n)(X) = sum_(j=0)^d binom(d,j)*gamma_(n+j)*X^j",
        },
        "transfer_theorem": {
            "range": "integers d,n >= 0, whenever M_(n),...,M_(n+d) are finite",
            "constant": "K_(d,n) = sqrt(pi)*d!/(4^n*Gamma(n+d+1/2))",
            "identity": (
                "J_gamma^(d,n)(X) = K_(d,n) integral u^(2n) "
                "L_d^(n-1/2)(-u^2*X/4) dmu(u)"
            ),
        },
        "coefficient_proof": {
            "laguerre_expansion": (
                "L_d^(alpha)(z) = sum_(j=0)^d "
                "Gamma(d+alpha+1)/[Gamma(alpha+j+1)(d-j)!j!]*(-z)^j"
            ),
            "raw_M_n_plus_j_coefficient": _q(raw_coefficient),
            "duplication_identity": (
                "sqrt(pi)/(4^m*Gamma(m+1/2)) = Gamma(m+1)/Gamma(2m+1), m=n+j"
            ),
            "after_duplication": _q(after_duplication),
            "target": "binomial(d,j)*(n+j)!/(2n+2j)!",
            "target_times_moment": "binomial(d,j)*gamma_(n+j)",
            "symbolic_ratio": _q(symbolic_ratio),
            "sympy_exact_grid": {"d_min": 0, "d_max": 8, "n_min": 0, "n_max": 6, "coefficient_checks": grid_count},
        },
        "xi_specialization": {
            "abstract_positive_measure": (
                "Psi(z)=xi(1/2+z)=integral cosh(z*u)dmu_Psi(u); "
                "Xi(t)=Psi(i*t)=integral cos(t*u)dmu_Psi(u)"
            ),
            "dbn_convention": (
                "xi(1/2+z)=8*integral_0^infinity Phi_DBN(U)*cosh(2zU)dU"
            ),
            "pushforward_measure": "dmu_Psi(u)=4*Phi_DBN(u/2)du",
            "moment_normalization": (
                "M_m=integral u^(2m)dmu_Psi(u)="
                "8*4^m*integral_0^infinity U^(2m)Phi_DBN(U)dU"
            ),
            "warning": (
                "Writing dmu=Phi(u)du is exact only if Phi denotes the rescaled Psi-kernel "
                "4*Phi_DBN(u/2); this avoids losing the factors 8 and 4^m."
            ),
        },
        "atomic_hyperbolicity": {
            "statement": (
                "For one atom u0>0, J is a positive constant times "
                "L_d^(n-1/2)(-u0^2 X/4), hence is hyperbolic."
            ),
            "roots": "X_k(u0)=-4*lambda_(d,k)^(n-1/2)/u0^2, k=1,...,d",
            "reason": "n-1/2>-1, so all generalized-Laguerre roots lambda_(d,k) are positive.",
        },
        "common_interlacing_audit": {
            "observation": (
                "The atomic root sets are u^(-2) dilations.  On continuous support in "
                "u>0, each root curve sweeps toward -infinity as u tends to 0 and toward "
                "0 from below as u tends to infinity."
            ),
            "conclusion": (
                "No fixed common interlacer is supplied by the representation itself, so the "
                "standard common-interlacing sufficient lemma for positive mixtures cannot be "
                "invoked automatically.  This is not a theorem that common interlacing or every "
                "other preservation mechanism is impossible."
            ),
        },
        "h790_canary": {
            "x_measure": "nu=(9/10)delta_1+(1/10)delta_16 in x=u^2",
            "equivalent_u_measure": "mu=(9/10)delta_1+(1/10)delta_4 in u",
            "moments_M0_M1_M2": [_q(value) for value in moments],
            "gamma_0_1_2": [_q(value) for value in gamma],
            "J_gamma_2_0": str(canary_jensen),
            "discriminant": _q(canary_discriminant),
            "verdict": (
                "The discriminant is negative, so a positive mixture of individually "
                "hyperbolic Laguerre atoms need not be hyperbolic."
            ),
        },
        "remaining_gap": {
            "statement": (
                "Prove, using structure special to the Xi measure, that the displayed Laguerre "
                "mixture is hyperbolic uniformly for every d,n>=0 (or prove an equivalent "
                "all-degree criterion)."
            ),
            "status": "open_not_claimed",
            "scope": (
                "H1324 is an exact transfer plus an exact obstruction to the generic mixing "
                "shortcut.  It does not prove Jensen hyperbolicity, PF-infinity, or RH."
            ),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-report", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate()
    if args.check_report:
        pinned = json.loads(REPORT_PATH.read_text(encoding="utf-8"))
        if pinned != certificate:
            raise SystemExit(f"certificate mismatch: {REPORT_PATH}")
    print(json.dumps(certificate, indent=2, sort_keys=True))
    return 0 if certificate["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
