#!/usr/bin/env python3
"""Dense finite-window scout for the full-Xi scaled third cumulant.

This is a numerical diagnostic, not an interval certificate.  It evaluates
dominant and full-Xi cumulants directly in the normalized saddle coordinate W,
without finite differences.  A fast stable float64 Gauss-Legendre pass scans
the compact H1319 r-range, while selected mpmath Gauss-Legendre passes test
precision, window, kernel cutoff, and quadrature-order convergence at the
sampled minimum.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import mpmath as mp
import numpy as np


ALPHA = 9.0 / 4.0
TARGET = -3.0 / 4.0
DEFAULT_OUT = Path(
    "research/riemann/h1319_full_xi_compact_scaled_cumulant_scout.json"
)
BASE_KNOTS = (-24, -18, -14, -10, -7, -5, -3, -2, -1, 0, 1, 2, 3, 5, 7, 10, 14, 18, 24)


def expm1_minus_x_np(x: np.ndarray) -> np.ndarray:
    """Stable exp(x)-1-x, including x near machine zero."""

    x = np.asarray(x, dtype=float)
    result = np.empty_like(x)
    small = np.abs(x) < 1.0e-4
    xs = x[small]
    result[small] = xs * xs * (
        0.5
        + xs
        * (
            1.0 / 6.0
            + xs
            * (
                1.0 / 24.0
                + xs * (1.0 / 120.0 + xs * (1.0 / 720.0 + xs / 5040.0))
            )
        )
    )
    xl = x[~small]
    result[~small] = np.expm1(xl) - xl
    return result


def quadrature_grid(width: float, order: int) -> tuple[np.ndarray, np.ndarray]:
    knots = sorted({-width, width, *(float(k) for k in BASE_KNOTS if abs(k) < width)})
    nodes, weights = np.polynomial.legendre.leggauss(order)
    all_nodes: list[np.ndarray] = []
    all_weights: list[np.ndarray] = []
    for left, right in zip(knots, knots[1:]):
        midpoint = (left + right) / 2.0
        half = (right - left) / 2.0
        all_nodes.append(midpoint + half * nodes)
        all_weights.append(half * weights)
    return np.concatenate(all_nodes), np.concatenate(all_weights)


def saddle_quantities_float(r: float) -> dict[str, float]:
    p = math.pi * math.exp(r)
    q = p * r - ALPHA * r - 1.0
    a2 = r * (p * (r + 1.0) - ALPHA)
    b3 = r * (p * (r * r + 3.0 * r + 1.0) - ALPHA)
    return {"P": p, "q": q, "t": q / 2.0, "A": a2, "B": b3}


def stable_potential_and_delta_np(
    r: float,
    w: np.ndarray,
    m_cutoff: int,
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    data = saddle_quantities_float(r)
    p = data["P"]
    sqrt_a = math.sqrt(data["A"])
    s = w / sqrt_a
    dr = r * np.expm1(s)

    # Exact identity with both saddle-cancelling differences evaluated stably:
    # Psi = P(expm1(dr)-dr) + (P-alpha)r(expm1(s)-s).
    psi = p * expm1_minus_x_np(dr) + (p - ALPHA) * r * expm1_minus_x_np(s)

    local_r = r + dr
    inv_x = np.exp(-local_r) / math.pi
    x_value = 1.0 / inv_x

    # H=Phi/full dominant ratio = 1+delta.  Keeping delta separate retains
    # corrections far below the unit roundoff of 1+delta.
    delta = -1.5 * inv_x
    for m in range(2, m_cutoff + 1):
        m2 = float(m * m)
        delta += (m2 * m2 - 1.5 * m2 * inv_x) * np.exp(
            -x_value * (m2 - 1.0)
        )
    return psi, delta, data


def cumulants_from_base_and_correction(
    base: list[float],
    correction: list[float],
    a2: float,
    t: float,
) -> dict[str, float]:
    raw = [value / base[0] for value in base]
    mean = raw[1]
    variance = raw[2] - mean * mean
    third = raw[3] - 3.0 * mean * raw[2] + 2.0 * mean**3

    epsilon = correction[0] / base[0]
    correction_raw = [value / base[0] for value in correction]
    delta_raw = [
        (correction_raw[j] - epsilon * raw[j]) / (1.0 + epsilon)
        for j in range(4)
    ]
    d1, d2, d3 = delta_raw[1], delta_raw[2], delta_raw[3]
    delta_variance = d2 - 2.0 * mean * d1 - d1 * d1
    delta_third = (
        d3
        - 3.0 * (mean * d2 + raw[2] * d1 + d1 * d2)
        + 6.0 * mean * mean * d1
        + 6.0 * mean * d1 * d1
        + 2.0 * d1**3
    )

    scale = 8.0 * t * t / (a2 ** 1.5)
    scaled_dominant = scale * third
    scaled_correction = scale * delta_third
    scaled_full = scaled_dominant + scaled_correction
    return {
        "epsilon": epsilon,
        "dominant_mean_W": mean,
        "dominant_variance_W": variance,
        "dominant_kappa3_W": third,
        "full_mean_W": mean + d1,
        "full_variance_W": variance + delta_variance,
        "full_kappa3_W": third + delta_third,
        "delta_kappa3_W": delta_third,
        "scaled_dominant_kappa3_t2": scaled_dominant,
        "scaled_correction_kappa3_t2": scaled_correction,
        "scaled_full_kappa3_t2": scaled_full,
        "target_margin": scaled_full - TARGET,
    }


def sample_float(
    r: float,
    width: float,
    m_cutoff: int,
    order: int,
    grid: tuple[np.ndarray, np.ndarray] | None = None,
) -> dict[str, float]:
    w, weights = grid if grid is not None else quadrature_grid(width, order)
    psi, delta, data = stable_potential_and_delta_np(r, w, m_cutoff)
    weighted_density = weights * np.exp(-psi)
    base = [float(np.sum(weighted_density * w**j)) for j in range(4)]
    correction = [
        float(np.sum(weighted_density * delta * w**j)) for j in range(4)
    ]
    moments = cumulants_from_base_and_correction(base, correction, data["A"], data["t"])
    return {"r": r, **data, **moments}


def saddle_r_for_t(t: mp.mpf, dps: int = 80) -> mp.mpf:
    mp.mp.dps = dps
    alpha = mp.mpf(9) / 4
    target_q = 2 * t
    equation = lambda r: mp.pi * r * mp.exp(r) - alpha * r - 1 - target_q
    return mp.findroot(equation, (mp.mpf(3), mp.mpf(4)))


def initial_grid(
    r_min: float,
    r_max: float,
    log_points: int,
    linear_points: int,
) -> list[float]:
    values = set(float(x) for x in np.geomspace(r_min, r_max, log_points))
    linear_right = min(r_max, 8.0)
    values.update(float(x) for x in np.linspace(r_min, linear_right, linear_points))
    values.update(
        x
        for x in (4.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, r_max)
        if r_min <= x <= r_max
    )
    values.add(r_min)
    values.add(r_max)
    return sorted(values)


def adaptive_scan(
    r_min: float,
    r_max: float,
    width: float,
    m_cutoff: int,
    order: int,
    log_points: int,
    linear_points: int,
    zoom_rounds: int,
    zoom_points: int,
) -> tuple[list[dict[str, float]], list[dict[str, Any]]]:
    grid = quadrature_grid(width, order)
    cache: dict[str, dict[str, float]] = {}

    def evaluate(points: list[float]) -> None:
        for r in points:
            key = float(r).hex()
            if key not in cache:
                cache[key] = sample_float(r, width, m_cutoff, order, grid)

    evaluate(initial_grid(r_min, r_max, log_points, linear_points))
    history: list[dict[str, Any]] = []
    for round_index in range(zoom_rounds):
        rows = sorted(cache.values(), key=lambda row: row["r"])
        minimum_index = min(
            range(len(rows)), key=lambda index: rows[index]["scaled_full_kappa3_t2"]
        )
        left_index = max(0, minimum_index - 1)
        right_index = min(len(rows) - 1, minimum_index + 1)
        if left_index == right_index:
            break
        left = rows[left_index]["r"]
        right = rows[right_index]["r"]
        history.append(
            {
                "round": round_index + 1,
                "minimum_r_before_zoom": rows[minimum_index]["r"],
                "minimum_value_before_zoom": rows[minimum_index][
                    "scaled_full_kappa3_t2"
                ],
                "zoom_left": left,
                "zoom_right": right,
            }
        )
        candidates = [
            float(x) for x in np.linspace(left, right, zoom_points + 2)[1:-1]
        ]
        evaluate(candidates)
    return sorted(cache.values(), key=lambda row: row["r"]), history


def mp_knots(width: mp.mpf) -> list[mp.mpf]:
    values = {-width, width}
    values.update(mp.mpf(k) for k in BASE_KNOTS if abs(k) < width)
    return sorted(values)


def mp_sample(
    r_text: str,
    dps: int,
    width_text: str,
    m_cutoff: int,
    order: int,
) -> dict[str, Any]:
    mp.mp.dps = dps
    r = mp.mpf(r_text)
    width = mp.mpf(width_text)
    alpha = mp.mpf(9) / 4
    p = mp.pi * mp.exp(r)
    q = p * r - alpha * r - 1
    t = q / 2
    a2 = r * (p * (r + 1) - alpha)
    root_a = mp.sqrt(a2)
    nodes, weights = mp.gauss_quadrature(order, "legendre")
    base = [mp.mpf(0) for _ in range(4)]
    correction = [mp.mpf(0) for _ in range(4)]

    for left, right in zip(mp_knots(width), mp_knots(width)[1:]):
        midpoint = (left + right) / 2
        half = (right - left) / 2
        for node, weight in zip(nodes, weights):
            w = midpoint + half * node
            quadrature_weight = half * weight
            s = w / root_a
            dr = r * mp.expm1(s)
            psi = p * (mp.expm1(dr) - dr) + (p - alpha) * r * (
                mp.expm1(s) - s
            )
            density = mp.exp(-psi)
            local_r = r + dr
            inv_x = mp.exp(-local_r) / mp.pi
            x_value = 1 / inv_x
            delta = -mp.mpf(3) * inv_x / 2
            for m in range(2, m_cutoff + 1):
                m2 = mp.mpf(m * m)
                delta += (m2 * m2 - mp.mpf(3) * m2 * inv_x / 2) * mp.exp(
                    -x_value * (m2 - 1)
                )
            mass = quadrature_weight * density
            power = mp.mpf(1)
            for j in range(4):
                base[j] += mass * power
                correction[j] += mass * delta * power
                power *= w

    raw = [value / base[0] for value in base]
    mean = raw[1]
    variance = raw[2] - mean * mean
    third = raw[3] - 3 * mean * raw[2] + 2 * mean**3
    epsilon = correction[0] / base[0]
    correction_raw = [value / base[0] for value in correction]
    delta_raw = [
        (correction_raw[j] - epsilon * raw[j]) / (1 + epsilon)
        for j in range(4)
    ]
    d1, d2, d3 = delta_raw[1], delta_raw[2], delta_raw[3]
    delta_variance = d2 - 2 * mean * d1 - d1 * d1
    delta_third = (
        d3
        - 3 * (mean * d2 + raw[2] * d1 + d1 * d2)
        + 6 * mean * mean * d1
        + 6 * mean * d1 * d1
        + 2 * d1**3
    )
    scale = 8 * t * t / a2 ** mp.mpf("1.5")
    scaled_dominant = scale * third
    scaled_correction = scale * delta_third
    scaled_full = scaled_dominant + scaled_correction
    digits = min(dps - 8, 60)
    return {
        "dps": dps,
        "window": width_text,
        "m_cutoff": m_cutoff,
        "order": order,
        "r": mp.nstr(r, digits),
        "t": mp.nstr(t, digits),
        "epsilon": mp.nstr(epsilon, digits),
        "dominant_mean_W": mp.nstr(mean, digits),
        "dominant_variance_W": mp.nstr(variance, digits),
        "dominant_kappa3_W": mp.nstr(third, digits),
        "full_mean_W": mp.nstr(mean + d1, digits),
        "full_variance_W": mp.nstr(variance + delta_variance, digits),
        "full_kappa3_W": mp.nstr(third + delta_third, digits),
        "scaled_dominant_kappa3_t2": mp.nstr(scaled_dominant, digits),
        "scaled_correction_kappa3_t2": mp.nstr(scaled_correction, digits),
        "scaled_full_kappa3_t2": mp.nstr(scaled_full, digits),
        "target_margin": mp.nstr(scaled_full + mp.mpf(3) / 4, digits),
    }


def convergence_at(r_text: str) -> dict[str, Any]:
    configs = [
        (50, "18", 8, 24, "precision_50"),
        (70, "18", 8, 24, "precision_70"),
        (70, "24", 8, 24, "window_24"),
        (70, "24", 12, 24, "cutoff_12"),
        (70, "24", 12, 32, "order_32"),
        (90, "24", 12, 32, "precision_90"),
        (90, "24", 16, 32, "cutoff_16"),
        (90, "28", 16, 32, "window_28_reference"),
    ]
    rows = []
    for dps, width, cutoff, order, label in configs:
        row = mp_sample(r_text, dps, width, cutoff, order)
        row["label"] = label
        rows.append(row)
    reference = mp.mpf(rows[-1]["scaled_full_kappa3_t2"])
    for row in rows:
        row["abs_difference_from_reference"] = mp.nstr(
            abs(mp.mpf(row["scaled_full_kappa3_t2"]) - reference), 20
        )
    return {
        "critical_r": r_text,
        "reference_label": rows[-1]["label"],
        "rows": rows,
    }


def compact_row(row: dict[str, float]) -> dict[str, float]:
    keys = (
        "r",
        "q",
        "t",
        "A",
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
    return {key: row[key] for key in keys}


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    r_min_mp = saddle_r_for_t(mp.mpf(125), 90)
    r_min = float(r_min_mp)
    rows, zoom_history = adaptive_scan(
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
    minimum = min(rows, key=lambda row: row["scaled_full_kappa3_t2"])
    dominant_minimum = min(rows, key=lambda row: row["scaled_dominant_kappa3_t2"])
    correction_minimum = min(rows, key=lambda row: row["scaled_correction_kappa3_t2"])
    monotone_tolerance = 5.0e-12
    sampled_non_decreasing = all(
        right["scaled_full_kappa3_t2"]
        >= left["scaled_full_kappa3_t2"] - monotone_tolerance
        for left, right in zip(rows, rows[1:])
    )
    convergence = None if args.skip_convergence else convergence_at(mp.nstr(r_min_mp, 70))
    return {
        "classification": "h1319_full_xi_compact_scaled_cumulant_dense_scout_not_proof",
        "statement_tested": {
            "range": "r(t=125)<=r<=71",
            "target": "t^2*kappa_Xi'''(t)>=-3/4",
            "finite_window_diagnostic_only": True,
        },
        "configuration": {
            "r_min": mp.nstr(r_min_mp, 70),
            "r_max": args.r_max,
            "window": args.window,
            "m_cutoff": args.m_cutoff,
            "gauss_legendre_order_per_segment": args.order,
            "log_points": args.log_points,
            "linear_points_on_low_band": args.linear_points,
            "zoom_rounds": args.zoom_rounds,
            "zoom_points_per_round": args.zoom_points,
            "finite_differences_used": False,
        },
        "summary": {
            "sample_count": len(rows),
            "target_holds_on_all_samples": minimum["target_margin"] > 0,
            "minimum_full": compact_row(minimum),
            "minimum_dominant": compact_row(dominant_minimum),
            "minimum_correction": compact_row(correction_minimum),
            "full_scaled_cumulant_sampled_non_decreasing": sampled_non_decreasing,
            "all_sampled_corrections_nonnegative": all(
                row["scaled_correction_kappa3_t2"] >= -1.0e-15 for row in rows
            ),
        },
        "zoom_history": zoom_history,
        "convergence_at_sampled_minimum": convergence,
        "rows": [compact_row(row) for row in rows],
        "limitations": [
            "float64 finite-window quadrature for the dense scan",
            "mpmath fixed quadrature only at the sampled critical point",
            "no interval arithmetic in r or W",
            "no rigorous analytic remainder for |W| beyond the finite window",
            "sampling cannot exclude an unsampled interior dip",
        ],
        "next_rigorous_step": "replace the sampled grid by adaptive Arb r-boxes, interval W quadrature, and analytic tails",
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
    parser.add_argument("--skip-convergence", action="store_true")
    args = parser.parse_args()
    report = build_report(args)
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
