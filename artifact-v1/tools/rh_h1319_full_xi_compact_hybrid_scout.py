#!/usr/bin/env python3
"""Canonical hybrid-precision H1319 full-Xi compact scout.

The low and moderate saddle range uses the stable float64 quadrature from the
preliminary H1319 scout.  From a configurable threshold onward, every sample
is recomputed with direct mpmath Gauss-Legendre moments.  This avoids the
float64 odd-moment cancellation that appears when r is large.  No finite
differences are used anywhere.

Diagnostic only: finite windows, point samples, and non-interval arithmetic do
not constitute a proof on the compact r interval.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import mpmath as mp

import rh_h1319_full_xi_compact_scaled_cumulant_scout as base


DEFAULT_OUT = Path(
    "research/riemann/h1319_full_xi_compact_scaled_cumulant_scout.json"
)


def mp_to_dense_row(r: float, result: dict[str, Any]) -> dict[str, float]:
    saddle = base.saddle_quantities_float(r)
    keys = (
        "epsilon",
        "dominant_mean_W",
        "dominant_variance_W",
        "dominant_kappa3_W",
        "full_mean_W",
        "full_variance_W",
        "full_kappa3_W",
        "scaled_dominant_kappa3_t2",
        "scaled_correction_kappa3_t2",
        "scaled_full_kappa3_t2",
        "target_margin",
    )
    return {
        "r": r,
        **saddle,
        **{key: float(mp.mpf(result[key])) for key in keys},
    }


def endpoint_convergence(r_text: str) -> dict[str, Any]:
    configs = (
        (50, "18", 8, 20, "endpoint_dps50_W18_m8_o20"),
        (70, "24", 12, 28, "endpoint_dps70_W24_m12_o28"),
        (90, "28", 16, 32, "endpoint_dps90_W28_m16_o32_reference"),
    )
    rows = []
    for dps, width, cutoff, order, label in configs:
        row = base.mp_sample(r_text, dps, width, cutoff, order)
        row["label"] = label
        rows.append(row)
    reference = mp.mpf(rows[-1]["scaled_full_kappa3_t2"])
    for row in rows:
        row["abs_difference_from_reference"] = mp.nstr(
            abs(mp.mpf(row["scaled_full_kappa3_t2"]) - reference), 24
        )
    return {"r": r_text, "reference": rows[-1]["label"], "rows": rows}


def hybrid_report(args: argparse.Namespace) -> dict[str, Any]:
    r_min_mp = base.saddle_r_for_t(mp.mpf(125), 90)
    r_min = float(r_min_mp)
    float_rows, zoom_history = base.adaptive_scan(
        r_min,
        args.r_max,
        args.window,
        args.m_cutoff,
        args.order,
        args.log_points,
        args.linear_points,
        args.zoom_rounds,
        args.zoom_points,
    )

    rows = []
    replacement_count = 0
    first_replacement_comparison = None
    for float_row in float_rows:
        r = float_row["r"]
        if r >= args.mp_threshold:
            mp_result = base.mp_sample(
                format(r, ".17g"),
                args.mp_dps,
                str(args.mp_window),
                args.mp_m_cutoff,
                args.mp_order,
            )
            row = mp_to_dense_row(r, mp_result)
            replacement_count += 1
            if first_replacement_comparison is None:
                first_replacement_comparison = {
                    "r": r,
                    "float64_scaled_full": float_row["scaled_full_kappa3_t2"],
                    "mpmath_scaled_full": row["scaled_full_kappa3_t2"],
                    "absolute_difference": abs(
                        float_row["scaled_full_kappa3_t2"]
                        - row["scaled_full_kappa3_t2"]
                    ),
                }
            rows.append(row)
        else:
            rows.append(float_row)

    rows.sort(key=lambda row: row["r"])
    minimum = min(rows, key=lambda row: row["scaled_full_kappa3_t2"])
    dominant_minimum = min(rows, key=lambda row: row["scaled_dominant_kappa3_t2"])
    correction_minimum = min(rows, key=lambda row: row["scaled_correction_kappa3_t2"])
    sampled_non_decreasing = all(
        right["scaled_full_kappa3_t2"]
        >= left["scaled_full_kappa3_t2"] - args.monotone_tolerance
        for left, right in zip(rows, rows[1:])
    )

    minimum_convergence = base.convergence_at(mp.nstr(r_min_mp, 70))
    high_convergence = endpoint_convergence(format(args.r_max, ".17g"))
    return {
        "classification": "h1319_full_xi_compact_hybrid_dense_scout_not_proof",
        "statement_tested": {
            "range": "r(t=125)<=r<=71",
            "target": "t^2*kappa_Xi'''(t)>=-3/4",
            "finite_window_diagnostic_only": True,
        },
        "configuration": {
            "r_min": mp.nstr(r_min_mp, 70),
            "r_max": args.r_max,
            "float64_engine": {
                "range": f"r<{args.mp_threshold}",
                "window": args.window,
                "m_cutoff": args.m_cutoff,
                "order": args.order,
            },
            "mpmath_engine": {
                "range": f"r>={args.mp_threshold}",
                "dps": args.mp_dps,
                "window": args.mp_window,
                "m_cutoff": args.mp_m_cutoff,
                "order": args.mp_order,
                "replacement_count": replacement_count,
            },
            "log_points": args.log_points,
            "linear_points_on_low_band": args.linear_points,
            "zoom_rounds": args.zoom_rounds,
            "zoom_points_per_round": args.zoom_points,
            "finite_differences_used": False,
        },
        "summary": {
            "sample_count": len(rows),
            "target_holds_on_all_samples": minimum["target_margin"] > 0,
            "minimum_full": base.compact_row(minimum),
            "minimum_dominant": base.compact_row(dominant_minimum),
            "minimum_correction": base.compact_row(correction_minimum),
            "full_scaled_cumulant_sampled_non_decreasing": sampled_non_decreasing,
            "all_sampled_corrections_nonnegative": all(
                row["scaled_correction_kappa3_t2"] >= 0 for row in rows
            ),
            "first_multiprecision_replacement": first_replacement_comparison,
        },
        "zoom_history": zoom_history,
        "convergence_at_sampled_minimum": minimum_convergence,
        "convergence_at_high_endpoint": high_convergence,
        "rows": [base.compact_row(row) for row in rows],
        "preliminary_float64_diagnostic": {
            "script": "tools/rh_h1319_full_xi_compact_scaled_cumulant_scout.py",
            "status": "retained for formulas and low-r engine; its all-float64 high-r odd cumulants are noncanonical",
        },
        "limitations": [
            "finite-window point quadrature, not interval arithmetic",
            "no interval coverage between sampled r values",
            "no rigorous analytic tail remainder beyond the W windows",
            "the float64-to-mpmath engine seam is numerically checked, not certified",
            "sampling cannot exclude an unsampled interior dip",
        ],
        "next_rigorous_step": "adaptive Arb r-box cover with interval W quadrature and analytic tails on [r(t=125),71]",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--r-max", type=float, default=71.0)
    parser.add_argument("--window", type=float, default=24.0)
    parser.add_argument("--m-cutoff", type=int, default=12)
    parser.add_argument("--order", type=int, default=32)
    parser.add_argument("--log-points", type=int, default=161)
    parser.add_argument("--linear-points", type=int, default=161)
    parser.add_argument("--zoom-rounds", type=int, default=4)
    parser.add_argument("--zoom-points", type=int, default=24)
    parser.add_argument("--mp-threshold", type=float, default=30.0)
    parser.add_argument("--mp-dps", type=int, default=60)
    parser.add_argument("--mp-window", type=int, default=20)
    parser.add_argument("--mp-m-cutoff", type=int, default=8)
    parser.add_argument("--mp-order", type=int, default=24)
    parser.add_argument("--monotone-tolerance", type=float, default=5.0e-12)
    args = parser.parse_args()

    report = hybrid_report(args)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "classification": report["classification"],
                "summary": report["summary"],
                "out": str(args.out),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
