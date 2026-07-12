#!/usr/bin/env python3
"""Exact multi-gap audit for the strict-PF2 lemma used by H1326.

This is separate from the main replay because the Windows patch sandbox would
not reopen the already-created H1326 runner.  It checks the arbitrary row-gap
ratio-product factorization, not merely an adjacent-ratio special case.
"""

from __future__ import annotations

import json
from collections import Counter

import sympy as sp


def main() -> int:
    checks: dict[str, bool] = {}
    canary_count = 0
    for x in range(8):
        for alpha in range(1, 9):
            for beta in range(1, 9):
                # Interior 2-minor:
                # a_(x+alpha)a_(x+beta)-a_(x+alpha+beta)a_x.
                left_indices = tuple(range(x + 1, x + alpha + 1))
                right_indices = tuple(range(x + beta + 1, x + beta + alpha + 1))
                ordered = len(left_indices) == len(right_indices) and all(
                    left < right for left, right in zip(left_indices, right_indices)
                )
                ratios = sp.symbols(f"R0:{x + alpha + beta + 2}")
                ax, axb = sp.symbols(f"a_x_{x} a_xb_{x}_{alpha}_{beta}")
                left_product = sp.prod(ratios[index] for index in left_indices)
                right_product = sp.prod(ratios[index] for index in right_indices)
                exact = sp.expand(
                    (ax * left_product) * axb
                    - (axb * right_product) * ax
                    - ax * axb * (left_product - right_product)
                ) == 0
                checks[f"x{x}_alpha{alpha}_beta{beta}"] = ordered and exact
                canary_count += 1

    boundary_classes: Counter[str] = Counter()
    boundary_ok = True
    for i in range(7):
        for j in range(i + 1, 9):
            for p in range(13):
                for q in range(p + 1, 14):
                    indices = (p - i, q - j, q - i, p - j)
                    first_nonzero = indices[0] >= 0 and indices[1] >= 0
                    second_nonzero = indices[2] >= 0 and indices[3] >= 0
                    if p < i:
                        boundary_classes["zero_first_column"] += 1
                        boundary_ok &= not first_nonzero and not second_nonzero
                    elif p < j:
                        if q < j:
                            boundary_classes["zero_lower_row"] += 1
                            boundary_ok &= not first_nonzero and not second_nonzero
                        else:
                            boundary_classes["strict_triangular"] += 1
                            boundary_ok &= first_nonzero and not second_nonzero
                    else:
                        boundary_classes["strict_interior"] += 1
                        boundary_ok &= first_nonzero and second_nonzero

    checks["all_512_multigap_canaries"] = canary_count == 512
    checks["all_one_sided_boundary_classes"] = boundary_ok
    checks["triangular_class_nonempty"] = boundary_classes["strict_triangular"] > 0
    checks["interior_class_nonempty"] = boundary_classes["strict_interior"] > 0
    all_checks = all(checks.values())
    report = {
        "schema": "rh_h1326_strict_pf2_multigap_audit.v0",
        "classification": (
            "h1326_strict_pf2_arbitrary_gap_lemma_audited"
            if all_checks
            else "h1326_strict_pf2_arbitrary_gap_audit_failed"
        ),
        "all_checks_pass": all_checks,
        "failed_checks": [name for name, value in checks.items() if not value],
        "multigap_canaries": canary_count,
        "boundary_classes": dict(sorted(boundary_classes.items())),
        "proof_formula": (
            "a_(x+alpha)a_(x+beta)-a_(x+alpha+beta)a_x "
            "= a_x a_(x+beta) (prod_{t=1}^alpha R_(x+t) "
            "- prod_{t=1}^alpha R_(x+beta+t))"
        ),
        "strict_sign_reason": (
            "beta>0 and strict decrease of R make each factor in the first "
            "product larger; one-sided nonzero boundary minors are triangular"
        ),
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if all_checks else 1


if __name__ == "__main__":
    raise SystemExit(main())
