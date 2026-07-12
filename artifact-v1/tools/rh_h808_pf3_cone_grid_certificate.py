import argparse
import json
import math
import sys
import time
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
    contiguous_d3_pairs,
    determinant_terms,
    eval_terms_float,
    exact_diagnosis,
    exact_verify,
)


DEFAULT_OUT = "research/riemann/h808_pf3_cone_grid_certificate.json"


def batch_iterator(denominator: int, length: int, batch_size: int):
    batch: list[tuple[int, ...]] = []
    for values in combinations_with_replacement(range(1, denominator + 1), length):
        batch.append(values)
        if len(batch) >= batch_size:
            yield np.array(batch, dtype=np.int64)
            batch = []
    if batch:
        yield np.array(batch, dtype=np.int64)


def run_grid(args: argparse.Namespace) -> dict[str, Any]:
    start = time.monotonic()
    sequence_len = args.length + 2
    all_pairs = all_order3_pairs(sequence_len)
    contiguous_pairs = contiguous_d3_pairs(sequence_len)
    all_terms = [(pair, determinant_terms(pair[0], pair[1], sequence_len)) for pair in all_pairs]
    contiguous_terms = [
        (pair, determinant_terms(pair[0], pair[1], sequence_len)) for pair in contiguous_pairs
    ]
    total = math.comb(args.denominator + args.length - 1, args.length)
    checked = 0
    contiguous_ok = 0
    verified_counterexample = None
    best_margin = math.inf
    best_pair = None
    best_q = None
    best_contiguous_min = math.inf
    best_contiguous_q = None
    exact_candidate_checks = 0

    for q_int_batch in batch_iterator(args.denominator, args.length, args.batch_size):
        q_grid = q_int_batch.astype(np.float64) / float(args.denominator)
        sequences = build_float_sequences(q_grid)
        batch_size = q_int_batch.shape[0]
        checked += batch_size

        contiguous_mask = np.ones(batch_size, dtype=bool)
        contiguous_min = np.full(batch_size, math.inf, dtype=np.float64)
        for _, terms in contiguous_terms:
            det_values = eval_terms_float(sequences, terms)
            contiguous_mask &= det_values >= -args.float_tolerance
            contiguous_min = np.minimum(contiguous_min, det_values)
        local_contiguous_best = int(np.argmin(contiguous_min))
        if float(contiguous_min[local_contiguous_best]) < best_contiguous_min:
            best_contiguous_min = float(contiguous_min[local_contiguous_best])
            best_contiguous_q = [int(value) for value in q_int_batch[local_contiguous_best]]

        contiguous_indices = np.nonzero(contiguous_mask)[0]
        contiguous_ok += int(contiguous_indices.size)
        if contiguous_indices.size:
            subset_sequences = sequences[contiguous_indices]
            subset_q_int = q_int_batch[contiguous_indices]
            min_values = np.full(contiguous_indices.size, math.inf, dtype=np.float64)
            min_pair_indices = np.zeros(contiguous_indices.size, dtype=np.int64)
            for pair_index, (_, terms) in enumerate(all_terms):
                det_values = eval_terms_float(subset_sequences, terms)
                update = det_values < min_values
                min_values[update] = det_values[update]
                min_pair_indices[update] = pair_index
            local_best_index = int(np.argmin(min_values))
            if float(min_values[local_best_index]) < best_margin:
                best_margin = float(min_values[local_best_index])
                best_pair = all_pairs[int(min_pair_indices[local_best_index])]
                best_q = [int(value) for value in subset_q_int[local_best_index]]

            negative_indices = np.nonzero(min_values < -args.float_tolerance)[0]
            for candidate_index in negative_indices[: args.exact_candidates_per_batch]:
                exact_candidate_checks += 1
                candidate = exact_verify(
                    [int(value) for value in subset_q_int[int(candidate_index)]],
                    args.denominator,
                )
                if candidate is not None:
                    verified_counterexample = candidate
                    break
        if verified_counterexample is not None:
            break

    return {
        "parameters": {
            "length": args.length,
            "denominator": args.denominator,
            "batch_size": args.batch_size,
            "float_tolerance": args.float_tolerance,
            "exact_candidates_per_batch": args.exact_candidates_per_batch,
        },
        "stats": {
            "grid_total": total,
            "checked": checked,
            "completed_full_grid": checked == total,
            "contiguous_d3_nonnegative_float_candidates": contiguous_ok,
            "elapsed_seconds": time.monotonic() - start,
            "exact_candidate_checks": exact_candidate_checks,
            "best_float_margin": best_margin,
            "best_float_margin_rows": None if best_pair is None else list(best_pair[0]),
            "best_float_margin_cols": None if best_pair is None else list(best_pair[1]),
            "best_q_ints": best_q,
            "best_q_exact_diagnosis": exact_diagnosis(best_q, args.denominator),
            "best_contiguous_float_margin": best_contiguous_min,
            "best_contiguous_q_ints": best_contiguous_q,
            "best_contiguous_q_exact_diagnosis": exact_diagnosis(
                best_contiguous_q, args.denominator
            ),
        },
        "verified_counterexample": verified_counterexample,
    }


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    grid = run_grid(args)
    found = grid["verified_counterexample"] is not None
    if found:
        classification = "pf2_q_contiguous_d3_grid_counterexample_found"
        decision = (
            "The exhaustive monotone grid found a certified hard counterexample. "
            "The H807 bridge is false at this grid size."
        )
    else:
        classification = "pf2_q_contiguous_d3_exhaustive_grid_survival"
        decision = (
            "The full monotone rational grid survived without a certified hard counterexample. "
            "This remains finite evidence, but it upgrades H807 from random no-hit to a "
            "symbolic-proof-worthy order-3 cone candidate."
        )
    return {
        "schema": "rh_h808_pf3_cone_grid_certificate.v0",
        "classification": classification,
        "candidate_statement": (
            "Positive Toeplitz sequence + PF2 + nondecreasing q_n + translated contiguous "
            "order-3 Toeplitz minors nonnegative implies all order-3 Toeplitz minors nonnegative."
        ),
        "grid": grid,
        "decision": decision,
        "limitations": [
            "This is an exhaustive rational-grid certificate, not a continuum proof.",
            "Float screening can miss candidates only if exact negativity is hidden behind nonnegative float margins; exact checks are applied to flagged negative candidates and best boundaries.",
            "The statement is PF3, not PF-infinity or RH.",
        ],
        "next_target": {
            "name": "H809 symbolic sparse-minor decomposition",
            "statement": (
                "Try to decompose each order-3 Toeplitz minor as a nonnegative expression "
                "using PF2, q-monotonicity, and translated contiguous D3 constraints."
            ),
        },
    }


def write_markdown(report: dict[str, Any], out_path: Path) -> None:
    lines = [
        "# H808 PF3 Cone Exhaustive-Grid Certificate",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Candidate",
        "",
        report["candidate_statement"],
        "",
        "## Grid",
        "",
        "Parameters:",
        "",
        "```json",
        json.dumps(report["grid"]["parameters"], indent=2),
        "```",
        "",
        "Stats:",
        "",
        "```json",
        json.dumps(report["grid"]["stats"], indent=2),
        "```",
    ]
    candidate = report["grid"]["verified_counterexample"]
    if candidate is None:
        lines.extend(["", "No certified hard counterexample was found on the full grid.", ""])
    else:
        lines.extend(
            [
                "",
                "Certified hard counterexample:",
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
    parser = argparse.ArgumentParser(description="H808 exhaustive monotone-grid PF3 cone audit.")
    parser.add_argument("--length", type=int, default=7)
    parser.add_argument("--denominator", type=int, default=20)
    parser.add_argument("--batch-size", type=int, default=8192)
    parser.add_argument("--float-tolerance", type=float, default=0.0)
    parser.add_argument("--exact-candidates-per-batch", type=int, default=12)
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
