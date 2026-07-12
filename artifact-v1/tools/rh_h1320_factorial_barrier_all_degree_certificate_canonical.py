#!/usr/bin/env python3
"""Canonical H1320 certificate runner with exact rational theta_2 rebuild."""

from __future__ import annotations

import argparse
import json

import sympy as sp

from rh_h1320_factorial_barrier_all_degree_certificate import (
    REPORT_PATH,
    build_certificate,
)


def c(index: int) -> int:
    return 2 * index * (2 * index - 1)


def build_canonical_certificate() -> dict[str, object]:
    certificate = build_certificate()
    theta2 = sp.Rational(c(1) * c(3), c(2) ** 2)
    reconstructed_ratio = sp.Rational(
        certificate["special_n_2_reconstructed"]["T_3_over_T_2"]
    )
    certificate["special_n_2_reconstructed"]["theta_2"] = str(theta2)
    certificate["checks"]["n2_factorial_ratio_reconstructed"] = bool(
        reconstructed_ratio == theta2
    )
    certificate["all_checks_pass"] = all(certificate["checks"].values())
    return certificate


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-report", action="store_true")
    args = parser.parse_args()
    certificate = build_canonical_certificate()
    if args.check_report:
        pinned = json.loads(REPORT_PATH.read_text(encoding="utf-8"))
        if pinned != certificate:
            raise SystemExit(f"certificate mismatch: {REPORT_PATH}")
    print(json.dumps(certificate, indent=2, sort_keys=True))
    return 0 if certificate["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
