"""Canonical runner for the alias-repaired finite translated-D3 certificate.

The materialized engine used ``arb.parent()`` to construct the scalar one,
but python-flint 0.8 does not expose that method.  This runner replaces only
that constructor detail; coefficient extraction, alias enlargement, D3
arithmetic, determinant cross-checks, and canaries remain unchanged.
"""

from typing import Any

import claude_d3_alias_repaired_finite_certificate as engine


def translated_d3(q_values: list[Any | None], n: int) -> Any:
    q_left = q_values[n]
    q_center = q_values[n + 1]
    q_right = q_values[n + 2]
    one = q_center * 0 + 1
    return (one - q_center) ** 2 - q_center**2 * (one - q_left) * (one - q_right)


engine.translated_d3 = translated_d3


if __name__ == "__main__":
    raise SystemExit(engine.main())
