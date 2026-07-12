#!/usr/bin/env python3
"""Canonical JSON runner for the exact H1320 counterexample certificate."""

from __future__ import annotations

import argparse
import json

import sympy as sp

from rh_h1320_factorial_barrier_all_degree_counterexample import (
    DEFAULT_REPORT,
    build_report,
)


def normalize_json(value: object) -> object:
    """Convert exact SymPy booleans while rejecting non-JSON arithmetic."""
    if value is sp.true:
        return True
    if value is sp.false:
        return False
    if isinstance(value, dict):
        return {str(key): normalize_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [normalize_json(item) for item in value]
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--check-report",
        action="store_true",
        help="Require the checked-in JSON report to equal the recomputed report.",
    )
    args = parser.parse_args()

    report = normalize_json(build_report())
    if args.check_report:
        checked_in = json.loads(DEFAULT_REPORT.read_text(encoding="utf-8"))
        if checked_in != report:
            raise SystemExit(f"report mismatch: {DEFAULT_REPORT}")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
