#!/usr/bin/env python3
"""Precision-safe runner for the independent H1319 mean-value audit."""

from __future__ import annotations

import mpmath as mp
from flint import arb

import rh_h1319_score_mean_value_independent_audit as audit


def precise_mp_ball(value: mp.mpf) -> arb:
    """Convert the 100-digit oracle without an artificial 80-digit cut."""

    return arb(mp.nstr(value, 100))


audit.mp_ball = precise_mp_ball


if __name__ == "__main__":
    raise SystemExit(audit.main())
