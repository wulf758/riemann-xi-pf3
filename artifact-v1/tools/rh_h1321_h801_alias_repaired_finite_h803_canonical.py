"""Canonical runner for the H1321 finite-H803 alias repair.

The first materialized engine serialized interval records before selecting the
worst row.  This runner repairs that reporting-only defect without changing the
coefficient extraction, alias bound, or H803 arithmetic.
"""

from typing import Any

import rh_h1321_h801_alias_repaired_finite_h803 as engine


def arb_record(value: Any) -> dict[str, Any]:
    return {
        "repr": str(value),
        "lower": str(value.lower()),
        "upper": str(value.upper()),
        "lower_float": float(value.lower()),
        "upper_float": float(value.upper()),
        "mid": str(value.mid()),
        "rad": str(value.rad()),
        "positive_lower_bound": bool(value.lower() > 0),
        "contains_zero": bool(value.lower() <= 0 <= value.upper()),
    }


def minimum_row(rows: list[dict[str, Any]], key: str) -> dict[str, Any]:
    return min(rows, key=lambda row: row[key]["lower_float"])


engine.arb_record = arb_record
engine.minimum_row = minimum_row


if __name__ == "__main__":
    raise SystemExit(engine.main())
