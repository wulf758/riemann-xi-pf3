"""Independent exact audit of the proposed two-sided-cumulant D3 tail lemma.

Conditional hypotheses (not proved by this checker):

    kappa(t) = Lambda(2t),
    kappa''(t)  <= 4/(5t),
    kappa'''(t) <= 11/(20t^2)        for t >= 125,
    q_n <= q_(n+1) <= q_(n+2).

The checker certifies the normalization, rational envelopes, monotonicity
identities, and shifted-polynomial signs for every n >= 127.  It does not
certify the two derivative hypotheses themselves.
"""

from __future__ import annotations

import json

import sympy as sp


N0 = 127


def require_zero(name: str, expression: sp.Expr) -> None:
    reduced = sp.factor(sp.cancel(sp.combsimp(expression)))
    if reduced != 0:
        raise AssertionError(f"{name}: expected zero, got {reduced}")


def shifted_coefficients(expression: sp.Expr, n: sp.Symbol, m: sp.Symbol) -> list[int]:
    polynomial = sp.Poly(sp.expand(expression.subs(n, m + N0)), m)
    coefficients = [int(value) for value in polynomial.all_coeffs()]
    if not coefficients or not all(value > 0 for value in coefficients):
        raise AssertionError(
            f"shifted coefficient positivity failed for {expression}: {coefficients}"
        )
    return coefficients


def main() -> None:
    n = sp.symbols("n", integer=True, positive=True)
    m = sp.symbols("m", integer=True, nonnegative=True)

    # Chain-rule audit for kappa(t)=Lambda(2t): factors 4 and 8.
    t = sp.symbols("t", real=True)
    for degree in range(9):
        Lambda = t**degree
        kappa = Lambda.subs(t, 2 * t)
        require_zero(
            f"kappa_second_chain_rule_degree_{degree}",
            sp.diff(kappa, t, 2) - 4 * sp.diff(Lambda, t, 2).subs(t, 2 * t),
        )
        require_zero(
            f"kappa_third_chain_rule_degree_{degree}",
            sp.diff(kappa, t, 3) - 8 * sp.diff(Lambda, t, 3).subs(t, 2 * t),
        )

    # Unit-step finite-difference averages in the kappa variable.
    base, r, s, w = sp.symbols("base r s w", real=True)
    for degree in range(9):
        F = t**degree
        delta2 = F.subs(t, base + 2) - 2 * F.subs(t, base + 1) + F.subs(t, base)
        average2 = sp.integrate(
            sp.diff(F, t, 2).subs(t, base + r + s),
            (r, 0, 1),
            (s, 0, 1),
        )
        require_zero(f"unit_delta2_degree_{degree}", delta2 - average2)

        delta3 = (
            F.subs(t, base + 3)
            - 3 * F.subs(t, base + 2)
            + 3 * F.subs(t, base + 1)
            - F.subs(t, base)
        )
        average3 = sp.integrate(
            sp.diff(F, t, 3).subs(t, base + r + s + w),
            (r, 0, 1),
            (s, 0, 1),
            (w, 0, 1),
        )
        require_zero(f"unit_delta3_degree_{degree}", delta3 - average3)

    def f(index: sp.Expr) -> sp.Expr:
        return sp.factor(
            (2 * index - 2)
            * (2 * index - 3)
            / ((2 * index) * (2 * index - 1))
        )

    # Check the factorial quotient directly from a_n=M(2n)/(2n)!.
    factorial_quotient = sp.combsimp(
        sp.factorial(2 * n - 2) ** 2
        / (sp.factorial(2 * n) * sp.factorial(2 * n - 4))
    )
    require_zero("factorial_quotient", factorial_quotient - f(n))

    u = sp.Rational(4, 5) / (n - 2)
    v = sp.Rational(11, 20) / (n - 2) ** 2
    X = sp.factor(f(n) / (1 - u))
    S = sp.factor((f(n + 1) / f(n)) / (1 - v))
    Y = sp.factor(S * X)
    H = sp.factor(1 - Y - Y**2 * (1 - X))

    expected_X = 5 * (n - 2) * (n - 1) * (2 * n - 3) / (
        n * (2 * n - 1) * (5 * n - 14)
    )
    expected_S = (
        20 * n**2 * (n - 2) ** 2 * (2 * n - 1) ** 2
        / (
            (n - 1)
            * (n + 1)
            * (2 * n - 3)
            * (2 * n + 1)
            * (20 * n**2 - 80 * n + 69)
        )
    )
    expected_Y = (
        100 * n * (n - 2) ** 3 * (2 * n - 1)
        / (
            (n + 1)
            * (2 * n + 1)
            * (5 * n - 14)
            * (20 * n**2 - 80 * n + 69)
        )
    )
    require_zero("X_formula", X - expected_X)
    require_zero("S_formula", S - expected_S)
    require_zero("Y_formula", Y - expected_Y)

    # Exact D3 lower step.  Set x=q_n, y=q_(n+1), z=q_(n+2).
    x, y, z = sp.symbols("x y z", positive=True)
    Hxy = 1 - y - y**2 * (1 - x)
    Dxyz = (1 - y) ** 2 - y**2 * (1 - x) * (1 - z)
    require_zero(
        "D_minus_monotone_lower_bound",
        Dxyz - (1 - y) * Hxy - y**2 * (1 - x) * (z - y),
    )

    # Coupling y<=Sx repairs the otherwise wrong one-variable monotonicity.
    Svar = sp.symbols("S", positive=True)
    h = sp.expand(Hxy.subs(y, Svar * x))
    h_derivative = Svar * (-1 + Svar * x * (3 * x - 2))
    require_zero("coupled_h_derivative", sp.diff(h, x) - h_derivative)
    require_zero("terminal_H", h.subs({x: X, Svar: S}) - H)

    # Denominators are positive for n>=127.  The only non-linear factor is
    # 20(n-2)^2-11, visibly positive already at n=127.
    denominator_facts = {
        "n": n,
        "2n-1": 2 * n - 1,
        "5n-14": 5 * n - 14,
        "n-1": n - 1,
        "n+1": n + 1,
        "2n-3": 2 * n - 3,
        "2n+1": 2 * n + 1,
        "20(n-2)^2-11": 20 * (n - 2) ** 2 - 11,
    }
    denominator_checks = {
        name: bool(expr.subs(n, N0) > 0) for name, expr in denominator_facts.items()
    }
    if not all(denominator_checks.values()):
        raise AssertionError(f"nonpositive denominator factor: {denominator_checks}")

    # Exact shifted-polynomial certificates n=m+127, m>=0.
    sign_targets = {
        "1-X": sp.factor(sp.together(1 - X).as_numer_denom()[0]),
        "S-1": sp.factor(sp.together(S - 1).as_numer_denom()[0]),
        "1-Y": sp.factor(sp.together(1 - Y).as_numer_denom()[0]),
        "H": sp.factor(sp.together(H).as_numer_denom()[0]),
    }
    shifted = {
        name: shifted_coefficients(expression, n, m)
        for name, expression in sign_targets.items()
    }

    report = {
        "schema": "claude_d3_two_sided_tail_lemma_independent_check.v0",
        "status": "ok",
        "conditional_hypotheses_proved_here": False,
        "conditional_hypotheses": [
            "kappa(t)=Lambda(2t)",
            "kappa''(t)<=4/(5t) for t>=125",
            "kappa'''(t)<=11/(20t^2) for t>=125",
            "q_n<=q_(n+1)<=q_(n+2)",
        ],
        "chain_rule": {"kappa_second_factor": 4, "kappa_third_factor": 8},
        "n_min": N0,
        "u": str(u),
        "v": str(v),
        "X": str(X),
        "S": str(S),
        "Y": str(Y),
        "H": str(H),
        "values_at_n_min": {
            "X": str(sp.N(X.subs(n, N0), 25)),
            "S": str(sp.N(S.subs(n, N0), 25)),
            "Y": str(sp.N(Y.subs(n, N0), 25)),
            "H": str(sp.N(H.subs(n, N0), 25)),
        },
        "denominator_checks": denominator_checks,
        "shift_n_equals_m_plus_127": shifted,
        "logical_chain": [
            "Delta^2 kappa(n-2)<=u and exp(u)<=1/(1-u) give q_n<=X",
            "Delta^3 kappa(n-2)<=v and exp(v)<=1/(1-v) give q_(n+1)/q_n<=S",
            "q monotone gives D_n>=(1-y)H(x,y)",
            "H decreases in y, so y<=Sx gives H(x,y)>=H(x,Sx)",
            "Y=SX<1 makes h_S'(x)<0 on 0<x<=X",
            "therefore H(x,y)>=H(X,Y)>0 and D_n>0",
        ],
        "conclusion": "conditional hypotheses imply D_n>0 for every n>=127",
        "global_Xi_D3_proved": False,
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
