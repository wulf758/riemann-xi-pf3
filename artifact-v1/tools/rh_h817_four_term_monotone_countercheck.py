import argparse
import json
import math
import sys
from fractions import Fraction
from itertools import combinations_with_replacement
from pathlib import Path
from typing import Any

import numpy as np

THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from rh_h813_sparse_row_symbolic_hierarchy import poly_to_string
from rh_h815_four_term_sparse_bridge import (
    certify_four_term_bridge,
    mixed_unclassified_pairs,
    monomial_text,
)


DEFAULT_OUT = "research/riemann/h817_four_term_monotone_countercheck.json"


def batch_iterator(denominator: int, length: int, batch_size: int):
    batch: list[tuple[int, ...]] = []
    for values in combinations_with_replacement(range(1, denominator + 1), length):
        batch.append(values)
        if len(batch) >= batch_size:
            yield np.array(batch, dtype=np.int64)
            batch = []
    if batch:
        yield np.array(batch, dtype=np.int64)


def h815_four_term_failures(length: int) -> list[dict[str, Any]]:
    out = []
    for item in mixed_unclassified_pairs(length):
        if item["analysis"]["term_count"] != 4:
            continue
        if certify_four_term_bridge(item["poly"]) is not None:
            continue
        out.append(item)
    return out


def eval_poly_float(q_values: np.ndarray, poly: dict[tuple[int, ...], int]) -> np.ndarray:
    values = np.zeros(q_values.shape[0], dtype=np.float64)
    for monomial, coeff in poly.items():
        term = np.ones(q_values.shape[0], dtype=np.float64)
        for index, exponent in enumerate(monomial):
            if exponent:
                term *= q_values[:, index] ** exponent
        values += coeff * term
    return values


def eval_poly_fraction(q_ints: list[int], denominator: int, poly: dict[tuple[int, ...], int]) -> Fraction:
    q_values = [Fraction(value, denominator) for value in q_ints]
    total = Fraction(0)
    for monomial, coeff in poly.items():
        term = Fraction(coeff)
        for index, exponent in enumerate(monomial):
            if exponent:
                term *= q_values[index] ** exponent
        total += term
    return total


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def serialize_counterexample(
    item: dict[str, Any],
    q_ints: list[int],
    denominator: int,
    determinant: Fraction,
) -> dict[str, Any]:
    return {
        "rows": list(item["rows"]),
        "cols": list(item["cols"]),
        "family": item["analysis"]["family"],
        "polynomial": item["analysis"]["polynomial"],
        "q_values": [
            {"name": f"q{index + 2}", "fraction": fraction_text(Fraction(value, denominator)), "float": value / denominator}
            for index, value in enumerate(q_ints)
        ],
        "determinant_polynomial_value": {
            "fraction": fraction_text(determinant),
            "float": float(determinant),
        },
    }


def run_grid(args: argparse.Namespace) -> dict[str, Any]:
    failures = h815_four_term_failures(args.length)
    total = math.comb(args.denominator + args.length - 1, args.length)
    checked = 0
    q2_half_candidates = 0
    negative_polynomial_indices: set[int] = set()
    counterexamples = []
    best = None

    for q_int_batch in batch_iterator(args.denominator, args.length, args.batch_size):
        checked += int(q_int_batch.shape[0])
        mask = q_int_batch[:, 0] <= args.denominator // 2
        if not np.any(mask):
            continue
        q_int_batch = q_int_batch[mask]
        q2_half_candidates += int(q_int_batch.shape[0])
        q_values = q_int_batch.astype(np.float64) / float(args.denominator)
        for index, item in enumerate(failures):
            values = eval_poly_float(q_values, item["poly"])
            local_index = int(np.argmin(values))
            local_value = float(values[local_index])
            if best is None or local_value < best["float_value"]:
                best = {
                    "float_value": local_value,
                    "failure_index": index,
                    "q_ints": [int(value) for value in q_int_batch[local_index]],
                    "rows": list(item["rows"]),
                    "cols": list(item["cols"]),
                    "polynomial": item["analysis"]["polynomial"],
                }
            if local_value < -args.float_tolerance and index not in negative_polynomial_indices:
                q_ints = [int(value) for value in q_int_batch[local_index]]
                exact = eval_poly_fraction(q_ints, args.denominator, item["poly"])
                if exact < 0:
                    negative_polynomial_indices.add(index)
                    if len(counterexamples) < args.example_limit:
                        counterexamples.append(serialize_counterexample(item, q_ints, args.denominator, exact))

    if best is not None:
        best_item = failures[best["failure_index"]]
        exact = eval_poly_fraction(best["q_ints"], args.denominator, best_item["poly"])
        best["exact_value"] = {
            "fraction": fraction_text(exact),
            "float": float(exact),
        }

    return {
        "parameters": {
            "length": args.length,
            "denominator": args.denominator,
            "batch_size": args.batch_size,
            "float_tolerance": args.float_tolerance,
        },
        "stats": {
            "grid_total": total,
            "checked": checked,
            "q2_le_half_candidates": q2_half_candidates,
            "h815_four_term_failures": len(failures),
            "negative_failure_polynomial_count": len(negative_polynomial_indices),
            "best_float_value": None if best is None else best["float_value"],
            "best_exact": best,
        },
        "counterexamples": counterexamples,
    }


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    grid = run_grid(args)
    negatives = grid["stats"]["negative_failure_polynomial_count"]
    if negatives:
        classification = "h817_monotone_only_four_term_counterexamples_found"
        decision = (
            "Some H815 four-term failures are genuinely negative under monotonicity and q2<=1/2 alone. "
            "Therefore Lemma B cannot be rescued as a monotone product lemma; D3, Plucker-Dodgson, "
            "or a Xi-specific analytic input is necessary."
        )
    else:
        classification = "h817_no_monotone_counterexample_on_grid"
        decision = (
            "No monotone-only counterexample was found on this grid. A stronger product lemma may exist, "
            "but this is finite evidence only."
        )
    return {
        "schema": "rh_h817_four_term_monotone_countercheck.v0",
        "classification": classification,
        "question": (
            "Are H815 four-term failures already nonnegative under 0<q2<=q3<=...<=1 and q2<=1/2, "
            "without using translated D3?"
        ),
        "grid": grid,
        "decision": decision,
        "next_target": {
            "name": "H818 Plucker-Dodgson or Xi-kernel pivot",
            "statement": (
                "If H817 finds monotone-only counterexamples, abandon monotone-only Lemma B and either "
                "derive a determinant identity using D3/H811 or pivot to the Xi-kernel unit-deficit route."
            ),
        },
    }


def write_markdown(report: dict[str, Any], out_path: Path) -> None:
    lines = [
        "# H817 Four-Term Monotone Countercheck",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Question",
        "",
        report["question"],
        "",
        "## Grid Stats",
        "",
        "```json",
        json.dumps(report["grid"]["stats"], indent=2),
        "```",
        "",
    ]
    if report["grid"]["counterexamples"]:
        lines.extend(["## Certified Counterexamples", ""])
        for item in report["grid"]["counterexamples"]:
            lines.extend(["```json", json.dumps(item, indent=2), "```", ""])
    lines.extend(["## Decision", "", report["decision"], "", "## Next Target", ""])
    lines.append(f"`{report['next_target']['name']}`: {report['next_target']['statement']}")
    lines.append("")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="H817 monotone-only countercheck for H815 failures.")
    parser.add_argument("--length", type=int, default=7)
    parser.add_argument("--denominator", type=int, default=16)
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument("--float-tolerance", type=float, default=1e-14)
    parser.add_argument("--example-limit", type=int, default=12)
    parser.add_argument("--out", default=DEFAULT_OUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = build_report(args)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    write_markdown(report, out_path.with_suffix(".md"))
    print(
        json.dumps(
            {
                "classification": report["classification"],
                "grid_stats": report["grid"]["stats"],
                "counterexample_count": len(report["grid"]["counterexamples"]),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
