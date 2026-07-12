"""Exact symbolic audit of the Xi translated-D3 / kappa_2 reduction.

This checker deliberately proves only algebraic and finite-difference
identities.  It does not assert D_n >= 0 for the Xi coefficients.
"""

from __future__ import annotations

import json

import sympy as sp


def require_zero(name: str, expression: sp.Expr) -> None:
    reduced = sp.factor(sp.cancel(sp.combsimp(expression)))
    if reduced != 0:
        raise AssertionError(f"{name}: expected zero, got {reduced}")


def main() -> None:
    n = sp.symbols("n", integer=True, positive=True)

    # a_n=M(2n)/(2n)! gives this exact factorial multiplier in q_n.
    factorial_multiplier = sp.combsimp(
        sp.factorial(2 * n - 2) ** 2
        / (sp.factorial(2 * n) * sp.factorial(2 * n - 4))
    )
    f_n = (2 * n - 2) * (2 * n - 3) / ((2 * n) * (2 * n - 1))
    require_zero("factorial_multiplier", factorial_multiplier - f_n)

    # Cross-check the third-difference factorial ratio used for q monotonicity.
    C = lambda k: (2 * k) * (2 * k - 1)
    b_ratio = sp.factor(C(n + 1) * C(n - 1) / C(n) ** 2)
    b_ratio_expected = (
        (2 * n + 2)
        * (2 * n + 1)
        * (2 * n - 2)
        * (2 * n - 3)
        / ((2 * n) * (2 * n - 1)) ** 2
    )
    require_zero("third_difference_factorial_ratio", b_ratio - b_ratio_expected)

    # The factor 4 belongs to Delta_2^2 Lambda; factor 8 belongs to
    # Delta_2^3 Lambda.  Polynomial bases make an exact symbolic oracle.
    t, base, x, y, z = sp.symbols("t base x y z", real=True)
    for degree in range(9):
        F = t**degree
        double_difference = (
            F.subs(t, base + 4) - 2 * F.subs(t, base + 2) + F.subs(t, base)
        )
        double_average = 4 * sp.integrate(
            sp.diff(F, t, 2).subs(t, base + 2 * (x + y)),
            (x, 0, 1),
            (y, 0, 1),
        )
        require_zero(f"double_difference_degree_{degree}", double_difference - double_average)

        triple_difference = (
            F.subs(t, base + 6)
            - 3 * F.subs(t, base + 4)
            + 3 * F.subs(t, base + 2)
            - F.subs(t, base)
        )
        triple_average = 8 * sp.integrate(
            sp.diff(F, t, 3).subs(t, base + 2 * (x + y + z)),
            (x, 0, 1),
            (y, 0, 1),
            (z, 0, 1),
        )
        require_zero(f"triple_difference_degree_{degree}", triple_difference - triple_average)

    # Interior contiguous Toeplitz determinant: start n, central quotient q_{n+1}.
    a_nm2, a_nm1, a_n, a_np1, a_np2 = sp.symbols(
        "a_nm2 a_nm1 a_n a_np1 a_np2", nonzero=True
    )
    q_n = a_n * a_nm2 / a_nm1**2
    q_np1 = a_np1 * a_nm1 / a_n**2
    q_np2 = a_np2 * a_n / a_np1**2
    determinant = sp.det(
        sp.Matrix(
            [
                [a_n, a_np1, a_np2],
                [a_nm1, a_n, a_np1],
                [a_nm2, a_nm1, a_n],
            ]
        )
    )
    d_deficit = (1 - q_np1) ** 2 - q_np1**2 * (1 - q_n) * (1 - q_np2)
    require_zero("toeplitz_determinant_normalization", determinant / a_n**3 - d_deficit)

    # Expanded quotient polynomial, useful for comparison with the H8xx tools.
    q0, q1, q2 = sp.symbols("q0 q1 q2", nonzero=True)
    d_abstract = (1 - q1) ** 2 - q1**2 * (1 - q0) * (1 - q2)
    d_expanded = 1 - 2 * q1 + q1**2 * (q0 + q2 - q0 * q2)
    require_zero("expanded_D", d_abstract - d_expanded)

    # Stable t=-log(q) form.  Algebraically set q_j=exp(-t_j) and replace
    # exp(t_{n+1}) by 1/q_{n+1}; no analytic approximation is involved.
    stable_bracket = (1 / q1 - 1) ** 2 - (1 - q0) * (1 - q2)
    require_zero("stable_t_form", d_abstract - q1**2 * stable_bracket)

    # The start-1 Toeplitz block is a boundary case because a_{-1}=0.
    a0, a1, a2, a3 = sp.symbols("a0 a1 a2 a3", nonzero=True)
    q_2 = a2 * a0 / a1**2
    q_3 = a3 * a1 / a2**2
    boundary_det = sp.det(sp.Matrix([[a1, a2, a3], [a0, a1, a2], [0, a0, a1]]))
    boundary_poly = 1 - 2 * q_2 + q_2**2 * q_3
    require_zero("start_1_boundary", boundary_det / a1**3 - boundary_poly)

    report = {
        "schema": "claude_d3_kappa2_exact_reduction_check.v0",
        "status": "ok",
        "proved_scope": "exact algebra and finite-difference factors only",
        "factorial_multiplier_f_n": str(sp.factor(f_n)),
        "inverse_multiplier_rho_n": str(sp.factor(1 / f_n)),
        "third_difference_factorial_ratio": str(b_ratio),
        "third_difference_ratio_minus_one": str(sp.factor(b_ratio - 1)),
        "difference_oracle_degrees": [0, 8],
        "double_difference_factor": 4,
        "triple_difference_factor": 8,
        "interior_toeplitz_normalization": "det(start=n)/a_n^3 = D_n, n>=2",
        "boundary_start_1": "det(start=1)/a_1^3 = 1-2*q_2+q_2^2*q_3",
        "global_D3_positivity_proved": False,
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
