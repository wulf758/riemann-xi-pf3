import argparse
import json
import math
import sys
import time
from fractions import Fraction
from itertools import combinations_with_replacement
from pathlib import Path
from typing import Any

import numpy as np

THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from rh_h807_essential_minor_hard_solver_audit import (
    all_order3_pairs,
    build_float_sequences,
    build_sequence_fraction,
    contiguous_d3_pairs,
    determinant_terms,
    eval_terms_float,
    serialize_fraction,
    toeplitz_minor_fraction,
)


DEFAULT_OUT = "research/riemann/h812_sparse_row_pf3_obstruction_map.json"


def is_consecutive_rows(rows: tuple[int, ...]) -> bool:
    return rows[1] == rows[0] + 1 and rows[2] == rows[1] + 1


def batch_iterator(denominator: int, length: int, batch_size: int):
    batch: list[tuple[int, ...]] = []
    for values in combinations_with_replacement(range(1, denominator + 1), length):
        batch.append(values)
        if len(batch) >= batch_size:
            yield np.array(batch, dtype=np.int64)
            batch = []
    if batch:
        yield np.array(batch, dtype=np.int64)


def q2_fraction(q_ints: list[int], denominator: int) -> Fraction:
    return Fraction(q_ints[0], denominator)


def exact_sparse_diagnosis(q_ints: list[int] | None, denominator: int) -> dict[str, Any] | None:
    if q_ints is None:
        return None
    q_values = [Fraction(value, denominator) for value in q_ints]
    sequence, ratios = build_sequence_fraction(q_values)
    contiguous_negative = []
    for rows, cols in contiguous_d3_pairs(len(sequence)):
        determinant = toeplitz_minor_fraction(sequence, rows, cols)
        if determinant < 0:
            contiguous_negative.append((rows, cols, determinant))

    consecutive_min: tuple[tuple[int, ...], tuple[int, ...], Fraction] | None = None
    sparse_min: tuple[tuple[int, ...], tuple[int, ...], Fraction] | None = None
    sparse_negative = []
    sparse_family_counts: dict[str, int] = {}
    sparse_structural_zero_count = 0
    for rows, cols in all_order3_pairs(len(sequence)):
        determinant = toeplitz_minor_fraction(sequence, rows, cols)
        if is_consecutive_rows(rows):
            if consecutive_min is None or determinant < consecutive_min[2]:
                consecutive_min = (rows, cols, determinant)
            continue
        if not determinant_terms(rows, cols, len(sequence)):
            sparse_structural_zero_count += 1
            continue
        if sparse_min is None or determinant < sparse_min[2]:
            sparse_min = (rows, cols, determinant)
        if determinant < 0:
            sparse_negative.append((rows, cols, determinant))
            family = f"row_gaps=({rows[1]-rows[0]},{rows[2]-rows[1]}), col_gaps=({cols[1]-cols[0]},{cols[2]-cols[1]})"
            sparse_family_counts[family] = sparse_family_counts.get(family, 0) + 1

    def serialize_minor(
        minor: tuple[tuple[int, ...], tuple[int, ...], Fraction] | None,
    ) -> dict[str, Any] | None:
        if minor is None:
            return None
        rows, cols, determinant = minor
        return {
            "rows": list(rows),
            "cols": list(cols),
            "determinant": serialize_fraction(determinant),
        }

    first_sparse_negative = sparse_negative[0] if sparse_negative else None
    return {
        "q_values": [serialize_fraction(value) for value in q_values],
        "q2": serialize_fraction(q_values[0]),
        "ratios": [serialize_fraction(value) for value in ratios],
        "sequence": [serialize_fraction(value) for value in sequence],
        "q2_le_half": bool(q2_fraction(q_ints, denominator) <= Fraction(1, 2)),
        "contiguous_d3_negative_count": len(contiguous_negative),
        "first_contiguous_d3_negative": serialize_minor(contiguous_negative[0])
        if contiguous_negative
        else None,
        "consecutive_row_minimum": serialize_minor(consecutive_min),
        "sparse_row_structural_zero_count": sparse_structural_zero_count,
        "sparse_row_negative_count": len(sparse_negative),
        "first_sparse_row_negative": serialize_minor(first_sparse_negative),
        "active_sparse_row_minimum": serialize_minor(sparse_min),
        "sparse_negative_family_counts": sparse_family_counts,
    }


def exact_verify_sparse(q_ints: list[int], denominator: int) -> dict[str, Any] | None:
    diagnosis = exact_sparse_diagnosis(q_ints, denominator)
    if diagnosis is None:
        return None
    if not diagnosis["q2_le_half"]:
        return None
    if diagnosis["contiguous_d3_negative_count"]:
        return None
    if diagnosis["sparse_row_negative_count"] <= 0:
        return None
    return diagnosis


def family_key(rows: tuple[int, ...], cols: tuple[int, ...]) -> str:
    return (
        f"row_gaps=({rows[1]-rows[0]},{rows[2]-rows[1]}), "
        f"col_gaps=({cols[1]-cols[0]},{cols[2]-cols[1]})"
    )


def run_grid(args: argparse.Namespace) -> dict[str, Any]:
    start = time.monotonic()
    sequence_len = args.length + 2
    all_pairs = all_order3_pairs(sequence_len)
    sparse_zero_pairs = []
    sparse_terms = []
    for rows, cols in all_pairs:
        if is_consecutive_rows(rows):
            continue
        terms = determinant_terms(rows, cols, sequence_len)
        if terms:
            sparse_terms.append(((rows, cols), terms))
        else:
            sparse_zero_pairs.append((rows, cols))
    sparse_pairs = [pair for pair, _ in sparse_terms]
    contiguous_pairs = contiguous_d3_pairs(sequence_len)
    contiguous_terms = [
        (pair, determinant_terms(pair[0], pair[1], sequence_len)) for pair in contiguous_pairs
    ]
    total = math.comb(args.denominator + args.length - 1, args.length)
    checked = 0
    q2_le_half = 0
    h811_float_candidates = 0
    verified_counterexample = None
    best_sparse_margin = math.inf
    best_sparse_pair = None
    best_sparse_q = None
    best_contiguous_margin = math.inf
    best_contiguous_q = None
    exact_candidate_checks = 0
    family_minima: dict[str, dict[str, Any]] = {}

    for q_int_batch in batch_iterator(args.denominator, args.length, args.batch_size):
        if args.max_seconds is not None and time.monotonic() - start > args.max_seconds:
            break

        checked += int(q_int_batch.shape[0])
        q2_mask = q_int_batch[:, 0] <= args.denominator // 2
        q2_le_half += int(np.count_nonzero(q2_mask))
        if not np.any(q2_mask):
            continue

        q_int_batch = q_int_batch[q2_mask]
        q_grid = q_int_batch.astype(np.float64) / float(args.denominator)
        sequences = build_float_sequences(q_grid)
        batch_size = q_int_batch.shape[0]
        contiguous_mask = np.ones(batch_size, dtype=bool)
        contiguous_min = np.full(batch_size, math.inf, dtype=np.float64)
        for _, terms in contiguous_terms:
            det_values = eval_terms_float(sequences, terms)
            contiguous_mask &= det_values >= -args.float_tolerance
            contiguous_min = np.minimum(contiguous_min, det_values)

        local_contiguous_best = int(np.argmin(contiguous_min))
        if float(contiguous_min[local_contiguous_best]) < best_contiguous_margin:
            best_contiguous_margin = float(contiguous_min[local_contiguous_best])
            best_contiguous_q = [int(value) for value in q_int_batch[local_contiguous_best]]

        contiguous_indices = np.nonzero(contiguous_mask)[0]
        h811_float_candidates += int(contiguous_indices.size)
        if not contiguous_indices.size:
            continue

        subset_sequences = sequences[contiguous_indices]
        subset_q_int = q_int_batch[contiguous_indices]
        min_values = np.full(contiguous_indices.size, math.inf, dtype=np.float64)
        min_pair_indices = np.zeros(contiguous_indices.size, dtype=np.int64)
        for pair_index, (pair, terms) in enumerate(sparse_terms):
            det_values = eval_terms_float(subset_sequences, terms)
            pair_min_index = int(np.argmin(det_values))
            pair_min_value = float(det_values[pair_min_index])
            key = family_key(pair[0], pair[1])
            current = family_minima.get(key)
            if current is None or pair_min_value < current["float_margin"]:
                family_minima[key] = {
                    "float_margin": pair_min_value,
                    "rows": list(pair[0]),
                    "cols": list(pair[1]),
                    "q_ints": [int(value) for value in subset_q_int[pair_min_index]],
                }
            update = det_values < min_values
            min_values[update] = det_values[update]
            min_pair_indices[update] = pair_index

        local_best_index = int(np.argmin(min_values))
        if float(min_values[local_best_index]) < best_sparse_margin:
            best_sparse_margin = float(min_values[local_best_index])
            best_sparse_pair = sparse_pairs[int(min_pair_indices[local_best_index])]
            best_sparse_q = [int(value) for value in subset_q_int[local_best_index]]

        negative_indices = np.nonzero(min_values < -args.float_tolerance)[0]
        for candidate_index in negative_indices[: args.exact_candidates_per_batch]:
            exact_candidate_checks += 1
            candidate = exact_verify_sparse(
                [int(value) for value in subset_q_int[int(candidate_index)]],
                args.denominator,
            )
            if candidate is not None:
                verified_counterexample = candidate
                break
        if verified_counterexample is not None:
            break

    sorted_families = sorted(family_minima.items(), key=lambda item: item[1]["float_margin"])
    return {
        "parameters": {
            "length": args.length,
            "denominator": args.denominator,
            "batch_size": args.batch_size,
            "float_tolerance": args.float_tolerance,
            "exact_candidates_per_batch": args.exact_candidates_per_batch,
            "max_seconds": args.max_seconds,
        },
        "shape_counts": {
            "all_order3_pairs": len(all_pairs),
            "consecutive_row_pairs_covered_by_h811": sum(1 for rows, _ in all_pairs if is_consecutive_rows(rows)),
            "sparse_row_pairs_total": len(sparse_pairs) + len(sparse_zero_pairs),
            "sparse_row_structural_zero_pairs": len(sparse_zero_pairs),
            "active_sparse_row_pairs_checked": len(sparse_pairs),
            "contiguous_d3_inputs_checked": len(contiguous_pairs),
        },
        "stats": {
            "grid_total": total,
            "checked": checked,
            "completed_full_grid": checked == total,
            "q2_le_half_float_candidates": q2_le_half,
            "h811_input_float_candidates": h811_float_candidates,
            "elapsed_seconds": time.monotonic() - start,
            "exact_candidate_checks": exact_candidate_checks,
            "best_sparse_float_margin": best_sparse_margin,
            "best_sparse_rows": None if best_sparse_pair is None else list(best_sparse_pair[0]),
            "best_sparse_cols": None if best_sparse_pair is None else list(best_sparse_pair[1]),
            "best_sparse_q_ints": best_sparse_q,
            "best_sparse_exact_diagnosis": exact_sparse_diagnosis(best_sparse_q, args.denominator),
            "best_contiguous_float_margin": best_contiguous_margin,
            "best_contiguous_q_ints": best_contiguous_q,
            "best_contiguous_exact_diagnosis": exact_sparse_diagnosis(best_contiguous_q, args.denominator),
            "top_sparse_family_minima": [
                {"family": key, **value} for key, value in sorted_families[: args.family_limit]
            ],
        },
        "verified_counterexample": verified_counterexample,
    }


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    grid = run_grid(args)
    found = grid["verified_counterexample"] is not None
    complete = grid["stats"]["completed_full_grid"]
    if found:
        classification = "h812_sparse_row_counterexample_found_under_h811_inputs"
        decision = (
            "The repaired H811 inputs are insufficient for PF3: a certified sparse-row "
            "order-3 Toeplitz counterexample was found. Pivot to adding an Xi-specific "
            "sparse-row lemma or abandon the PF2+q+D3 bridge."
        )
    elif complete:
        classification = "h812_sparse_row_grid_survival_under_h811_inputs"
        decision = (
            "No certified sparse-row obstruction was found on the full rational grid. "
            "Combined with H811, this makes the remaining PF3 bridge symbolically plausible "
            "at order 3, but still finite evidence only."
        )
    else:
        classification = "h812_sparse_row_grid_partial_no_counterexample"
        decision = (
            "The sparse-row scan stopped before completing the grid and did not certify a "
            "counterexample. Treat this only as partial orientation."
        )
    return {
        "schema": "rh_h812_sparse_row_pf3_obstruction_map.v0",
        "classification": classification,
        "candidate_statement": (
            "Under H811 inputs (PF2, nondecreasing q_n, q2<=1/2, and all translated "
            "contiguous order-3 Toeplitz minors nonnegative), the only remaining PF3 "
            "obstructions are sparse-row order-3 minors; search those directly."
        ),
        "h811_quotient": {
            "covered": "all consecutive-row order-3 Toeplitz minors",
            "remaining": "nonconsecutive-row order-3 Toeplitz minors",
        },
        "grid": grid,
        "decision": decision,
        "limitations": [
            "This is a finite rational-grid obstruction map, not a continuum proof.",
            "The scan assumes H811 as a proved conditional slice and only checks sparse-row minors.",
            "Even full PF3 would remain far weaker than PF-infinity/RH.",
        ],
        "next_target": {
            "name": "H813 sparse-row symbolic hierarchy",
            "statement": (
                "Use the top H812 sparse-row boundary families to seek an anchored or "
                "planar-orientation lemma in standard ratio variables."
            ),
        },
    }


def write_markdown(report: dict[str, Any], out_path: Path) -> None:
    grid = report["grid"]
    lines = [
        "# H812 Sparse-Row PF3 Obstruction Map",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Candidate",
        "",
        report["candidate_statement"],
        "",
        "## H811 Quotient",
        "",
        f"Covered: {report['h811_quotient']['covered']}",
        "",
        f"Remaining: {report['h811_quotient']['remaining']}",
        "",
        "## Grid",
        "",
        "Parameters:",
        "",
        "```json",
        json.dumps(grid["parameters"], indent=2),
        "```",
        "",
        "Shape counts:",
        "",
        "```json",
        json.dumps(grid["shape_counts"], indent=2),
        "```",
        "",
        "Stats:",
        "",
        "```json",
        json.dumps(grid["stats"], indent=2),
        "```",
    ]
    candidate = grid["verified_counterexample"]
    if candidate is None:
        lines.extend(["", "No certified sparse-row counterexample was found in this scan.", ""])
    else:
        lines.extend(
            [
                "",
                "Certified sparse-row counterexample:",
                "",
                "```json",
                json.dumps(candidate, indent=2),
                "```",
            ]
        )
    lines.extend(["", "## Decision", "", report["decision"], "", "## Limitations", ""])
    for limitation in report["limitations"]:
        lines.append(f"- {limitation}")
    lines.extend(
        [
            "",
            "## Next Target",
            "",
            f"`{report['next_target']['name']}`: {report['next_target']['statement']}",
            "",
        ]
    )
    out_path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="H812 sparse-row PF3 obstruction map.")
    parser.add_argument("--length", type=int, default=7)
    parser.add_argument("--denominator", type=int, default=16)
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument("--float-tolerance", type=float, default=0.0)
    parser.add_argument("--exact-candidates-per-batch", type=int, default=12)
    parser.add_argument("--family-limit", type=int, default=12)
    parser.add_argument("--max-seconds", type=float, default=90.0)
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
                "verified_counterexample": report["grid"]["verified_counterexample"] is not None,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

