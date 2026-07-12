import argparse
import json
import math
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np


THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from rh_h807_essential_minor_hard_solver_audit import (
    build_sequence_fraction,
    fraction_text,
    toeplitz_minor_fraction,
)
from rh_h815_four_term_sparse_bridge import (
    certify_four_term_bridge,
    mixed_unclassified_pairs,
)


DEFAULT_OUT = "research/riemann/h817b_h811_sparse_nonlinear_countercheck.json"


def remaining_pairs(length: int) -> list[dict[str, Any]]:
    return [
        item
        for item in mixed_unclassified_pairs(length)
        if certify_four_term_bridge(item["poly"]) is None
    ]


def edge_biased_unit(rng: np.random.Generator, shape: tuple[int, ...]) -> np.ndarray:
    """Sample [0,1], mixing bulk and both boundary layers."""

    raw = rng.random(shape)
    mode = rng.integers(0, 5, size=shape)
    out = raw.copy()
    out[mode == 1] = raw[mode == 1] ** 2
    out[mode == 2] = raw[mode == 2] ** 5
    out[mode == 3] = 1.0 - (1.0 - raw[mode == 3]) ** 2
    out[mode == 4] = 1.0 - (1.0 - raw[mode == 4]) ** 5
    return out


def q_from_slacks_float(slacks: np.ndarray) -> np.ndarray:
    """Map seven free slacks to q2,...,q8 while enforcing H811 inputs.

    We write e_n=1-q_n.  Monotonic q is e_{n+1}<=e_n and translated
    D3 is

      e_{n+1}^2 >= (1-e_{n+1})^2 e_n e_{n+2}.

    Hence choosing each next deficit below the displayed rational upper
    bound enforces every D3 constraint by construction.
    """

    batch = slacks.shape[0]
    deficits = np.empty((batch, 7), dtype=np.float64)
    deficits[:, 0] = 0.5 + 0.5 * slacks[:, 0]  # e2>=1/2, hence q2<=1/2.
    deficits[:, 1] = deficits[:, 0] * slacks[:, 1]
    tiny = np.finfo(np.float64).tiny
    for index in range(2, 7):
        previous = deficits[:, index - 2]
        current = deficits[:, index - 1]
        q_current = np.maximum(1.0 - current, tiny)
        d3_bound = current * current / (q_current * q_current * previous)
        upper = np.minimum(current, d3_bound)
        deficits[:, index] = upper * slacks[:, index]
    return 1.0 - deficits


def q_from_slacks_fraction(slacks: list[Fraction]) -> list[Fraction]:
    deficits = [Fraction(1, 2) + Fraction(1, 2) * slacks[0]]
    deficits.append(deficits[0] * slacks[1])
    for index in range(2, 7):
        previous = deficits[index - 2]
        current = deficits[index - 1]
        q_current = 1 - current
        if q_current <= 0:
            raise ValueError("the rationalized slack produced a zero q")
        d3_bound = current * current / (q_current * q_current * previous)
        upper = min(current, d3_bound)
        deficits.append(upper * slacks[index])
    return [1 - deficit for deficit in deficits]


def eval_poly_float(q_values: np.ndarray, poly: dict[tuple[int, ...], int]) -> np.ndarray:
    value = np.zeros(q_values.shape[0], dtype=np.float64)
    for monomial, coefficient in poly.items():
        term = np.ones(q_values.shape[0], dtype=np.float64)
        for index, exponent in enumerate(monomial):
            if exponent:
                term *= q_values[:, index] ** exponent
        value += coefficient * term
    return value


def eval_poly_fraction(q_values: list[Fraction], poly: dict[tuple[int, ...], int]) -> Fraction:
    value = Fraction(0)
    for monomial, coefficient in poly.items():
        term = Fraction(coefficient)
        for index, exponent in enumerate(monomial):
            if exponent:
                term *= q_values[index] ** exponent
        value += term
    return value


def normalized_d3_fraction(q_values: list[Fraction]) -> list[Fraction]:
    out = [1 - 2 * q_values[0] + q_values[0] * q_values[0] * q_values[1]]
    for index in range(0, len(q_values) - 2):
        x, y, z = q_values[index : index + 3]
        out.append((1 - y) ** 2 - y * y * (1 - x) * (1 - z))
    return out


def verify_exact(
    item: dict[str, Any], q_values: list[Fraction]
) -> dict[str, Any] | None:
    if not q_values or not (0 < q_values[0] <= Fraction(1, 2)):
        return None
    if any(not (0 < value <= 1) for value in q_values):
        return None
    if any(left > right for left, right in zip(q_values, q_values[1:])):
        return None
    d3 = normalized_d3_fraction(q_values)
    if any(value < 0 for value in d3):
        return None
    normalized = eval_poly_fraction(q_values, item["poly"])
    if normalized >= 0:
        return None
    sequence, ratios = build_sequence_fraction(q_values)
    determinant = toeplitz_minor_fraction(sequence, item["rows"], item["cols"])
    if determinant >= 0:
        raise AssertionError("normalized polynomial and exact determinant disagree")
    return {
        "rows": list(item["rows"]),
        "cols": list(item["cols"]),
        "family": item["analysis"]["family"],
        "polynomial": item["analysis"]["polynomial"],
        "q_values": [fraction_text(value) for value in q_values],
        "ratios": [fraction_text(value) for value in ratios],
        "sequence": [fraction_text(value) for value in sequence],
        "normalized_d3": [fraction_text(value) for value in d3],
        "normalized_sparse_polynomial": fraction_text(normalized),
        "sparse_toeplitz_determinant": fraction_text(determinant),
        "infinite_extension": "Set q_n=1 for n>=9; all later translated D3 constraints are then nonnegative.",
    }


def rationalize_candidate(
    item: dict[str, Any], slack_values: list[float], max_denominator: int
) -> dict[str, Any] | None:
    denominator = 32
    while denominator <= max_denominator:
        rational_slacks = [
            Fraction(value).limit_denominator(denominator) for value in slack_values
        ]
        if rational_slacks[0] >= 1:
            rational_slacks[0] = Fraction(denominator - 1, denominator)
        q_values = q_from_slacks_fraction(rational_slacks)
        witness = verify_exact(item, q_values)
        if witness is not None:
            witness["slacks"] = [fraction_text(value) for value in rational_slacks]
            witness["rationalization_denominator_limit"] = denominator
            return witness
        denominator *= 2
    return None


def run_search(args: argparse.Namespace) -> dict[str, Any]:
    pairs = remaining_pairs(args.length)
    rng = np.random.default_rng(args.seed)
    start = time.monotonic()
    best: dict[str, Any] | None = None
    witness = None
    evaluated = 0
    negative_float_candidates = 0

    for batch_index in range(args.batches):
        slacks = edge_biased_unit(rng, (args.batch_size, 7))
        # Deterministic boundary probes are included in every batch.  Values
        # remain strictly inside (0,1), which leaves room for rationalization.
        probes = min(16, args.batch_size)
        if probes:
            exponents = np.linspace(2.0, 14.0, probes)
            slacks[:probes] = 1.0 - 10.0 ** (-exponents[:, None])
            slacks[:probes, 0] = np.linspace(0.0, 0.999, probes)
        q_values = q_from_slacks_float(slacks)
        evaluated += q_values.shape[0]

        for pair_index, item in enumerate(pairs):
            values = eval_poly_float(q_values, item["poly"])
            local_index = int(np.argmin(values))
            local_value = float(values[local_index])
            if best is None or local_value < best["value"]:
                best = {
                    "value": local_value,
                    "pair_index": pair_index,
                    "rows": list(item["rows"]),
                    "cols": list(item["cols"]),
                    "family": item["analysis"]["family"],
                    "polynomial": item["analysis"]["polynomial"],
                    "q_values": q_values[local_index].tolist(),
                    "slacks": slacks[local_index].tolist(),
                    "batch_index": batch_index,
                }
            negative_indices = np.flatnonzero(values < -args.float_tolerance)
            if negative_indices.size:
                negative_float_candidates += int(negative_indices.size)
                for candidate_index in negative_indices[: args.rationalize_per_hit]:
                    witness = rationalize_candidate(
                        item,
                        slacks[int(candidate_index)].tolist(),
                        args.max_denominator,
                    )
                    if witness is not None:
                        break
            if witness is not None:
                break
        if witness is not None:
            break

    elapsed = time.monotonic() - start
    return {
        "parameters": {
            "length": args.length,
            "batch_size": args.batch_size,
            "batches": args.batches,
            "seed": args.seed,
            "float_tolerance": args.float_tolerance,
            "max_denominator": args.max_denominator,
        },
        "counts": {
            "h815_remaining_polynomials": len(pairs),
            "feasible_points_evaluated": evaluated,
            "negative_float_candidates": negative_float_candidates,
        },
        "elapsed_seconds": elapsed,
        "best_float_candidate": best,
        "exact_witness": witness,
    }


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    search = run_search(args)
    found = search["exact_witness"] is not None
    return {
        "schema": "rh_h817b_h811_sparse_nonlinear_countercheck.v0",
        "classification": (
            "h811_scalar_package_exact_sparse_counterexample_found"
            if found
            else "h811_scalar_package_nonlinear_search_no_counterexample"
        ),
        "question": (
            "Do positivity, PF2, q2<=1/2, nondecreasing q, and every translated "
            "contiguous D3 minor imply full PF3?"
        ),
        "method": (
            "Search the H815/H816 unresolved normalized determinant polynomials in a "
            "deficit/slack parametrization that enforces monotone q and D3 exactly by "
            "construction; rationalize every negative candidate and recheck with Fraction."
        ),
        "search": search,
        "decision": (
            "The abstract implication is false; any full-PF3 proof for Xi needs an additional Xi-specific input."
            if found
            else "No falsification was found in this nonlinear feasible-set search; this is evidence, not a proof."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--length", type=int, default=7)
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument("--batches", type=int, default=100)
    parser.add_argument("--seed", type=int, default=8172)
    parser.add_argument("--float-tolerance", type=float, default=1e-12)
    parser.add_argument("--rationalize-per-hit", type=int, default=8)
    parser.add_argument("--max-denominator", type=int, default=1 << 16)
    parser.add_argument("--out", default=DEFAULT_OUT)
    args = parser.parse_args()
    report = build_report(args)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({
        "classification": report["classification"],
        "counts": report["search"]["counts"],
        "best_float_value": None if report["search"]["best_float_candidate"] is None else report["search"]["best_float_candidate"]["value"],
        "exact_witness": report["search"]["exact_witness"],
    }, indent=2))


if __name__ == "__main__":
    main()
