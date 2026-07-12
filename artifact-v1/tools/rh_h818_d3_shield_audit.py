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

from rh_h813_sparse_row_symbolic_hierarchy import determinant_q_poly
from rh_h817_four_term_monotone_countercheck import (
    eval_poly_float,
    eval_poly_fraction,
    fraction_text,
    h815_four_term_failures,
)


DEFAULT_OUT = "research/riemann/h818_d3_shield_audit.json"


def batch_iterator(denominator: int, length: int, batch_size: int):
    batch: list[tuple[int, ...]] = []
    for values in combinations_with_replacement(range(1, denominator + 1), length):
        batch.append(values)
        if len(batch) >= batch_size:
            yield np.array(batch, dtype=np.int64)
            batch = []
    if batch:
        yield np.array(batch, dtype=np.int64)


def d3_polys(length: int) -> list[dict[str, Any]]:
    q_max = length + 1
    sequence_len = length + 2
    out = []
    for start in range(1, sequence_len - 2):
        poly, _ = determinant_q_poly((0, 1, 2), (start, start + 1, start + 2), q_max)
        out.append({"start": start, "poly": poly})
    return out


def eval_poly_fraction_named(q_ints: list[int], denominator: int, poly: dict[tuple[int, ...], int]) -> dict[str, Any]:
    value = eval_poly_fraction(q_ints, denominator, poly)
    return {"fraction": fraction_text(value), "float": float(value), "negative": value < 0}


def q_values_record(q_ints: list[int], denominator: int) -> list[dict[str, Any]]:
    return [
        {
            "name": f"q{index + 2}",
            "fraction": fraction_text(Fraction(value, denominator)),
            "float": value / denominator,
        }
        for index, value in enumerate(q_ints)
    ]


def run_grid(args: argparse.Namespace) -> dict[str, Any]:
    failures = h815_four_term_failures(args.length)
    d3s = d3_polys(args.length)
    total = math.comb(args.denominator + args.length - 1, args.length)
    checked = 0
    q2_half_candidates = 0

    records = [
        {
            "failure_index": index,
            "rows": list(item["rows"]),
            "cols": list(item["cols"]),
            "family": item["analysis"]["family"],
            "polynomial": item["analysis"]["polynomial"],
            "negative_count": 0,
            "negative_with_all_d3_nonnegative": 0,
            "single_d3_nonnegative_survivor_counts": {str(d3["start"]): 0 for d3 in d3s},
            "best_negative": None,
        }
        for index, item in enumerate(failures)
    ]

    for q_int_batch in batch_iterator(args.denominator, args.length, args.batch_size):
        checked += int(q_int_batch.shape[0])
        mask = q_int_batch[:, 0] <= args.denominator // 2
        if not np.any(mask):
            continue
        q_int_batch = q_int_batch[mask]
        q2_half_candidates += int(q_int_batch.shape[0])
        q_values = q_int_batch.astype(np.float64) / float(args.denominator)
        d3_values = np.column_stack([eval_poly_float(q_values, d3["poly"]) for d3 in d3s])
        all_d3_mask = np.all(d3_values >= -args.float_tolerance, axis=1)

        for index, item in enumerate(failures):
            values = eval_poly_float(q_values, item["poly"])
            negative_mask = values < -args.float_tolerance
            if not np.any(negative_mask):
                continue
            record = records[index]
            record["negative_count"] += int(np.count_nonzero(negative_mask))
            record["negative_with_all_d3_nonnegative"] += int(np.count_nonzero(negative_mask & all_d3_mask))
            for d3_index, d3 in enumerate(d3s):
                survivors = negative_mask & (d3_values[:, d3_index] >= -args.float_tolerance)
                record["single_d3_nonnegative_survivor_counts"][str(d3["start"])] += int(np.count_nonzero(survivors))

            local_indices = np.nonzero(negative_mask)[0]
            best_local = int(local_indices[np.argmin(values[negative_mask])])
            best_value = float(values[best_local])
            if record["best_negative"] is None or best_value < record["best_negative"]["float_value"]:
                q_ints = [int(value) for value in q_int_batch[best_local]]
                d3_exact = {
                    f"D3_start_{d3['start']}": eval_poly_fraction_named(q_ints, args.denominator, d3["poly"])
                    for d3 in d3s
                }
                exact = eval_poly_fraction(q_ints, args.denominator, item["poly"])
                record["best_negative"] = {
                    "float_value": best_value,
                    "exact_value": {"fraction": fraction_text(exact), "float": float(exact)},
                    "q_values": q_values_record(q_ints, args.denominator),
                    "d3_values": d3_exact,
                    "violated_d3_starts": [
                        d3["start"]
                        for d3 in d3s
                        if d3_exact[f"D3_start_{d3['start']}"]["negative"]
                    ],
                }

    negative_records = [record for record in records if record["negative_count"] > 0]
    for record in negative_records:
        shields = [
            int(start)
            for start, count in record["single_d3_nonnegative_survivor_counts"].items()
            if count == 0
        ]
        record["single_d3_shields"] = shields
        record["has_single_d3_shield"] = bool(shields)
        record["all_grid_negatives_blocked_by_all_d3"] = record["negative_with_all_d3_nonnegative"] == 0

    single_shield_counts: dict[str, int] = {str(d3["start"]): 0 for d3 in d3s}
    for record in negative_records:
        for start in record["single_d3_shields"]:
            single_shield_counts[str(start)] += 1

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
            "failures_with_monotone_negative": len(negative_records),
            "failures_with_negative_surviving_all_d3": sum(
                1 for record in negative_records if record["negative_with_all_d3_nonnegative"] > 0
            ),
            "failures_with_single_d3_shield": sum(1 for record in negative_records if record["has_single_d3_shield"]),
            "single_d3_shield_counts": single_shield_counts,
        },
        "negative_records": negative_records,
    }


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    grid = run_grid(args)
    stats = grid["stats"]
    if stats["failures_with_negative_surviving_all_d3"]:
        classification = "h818_d3_all_constraints_not_enough_on_grid"
        decision = (
            "Unexpectedly, some monotone negative four-term examples survive all D3 constraints on this grid. "
            "This would contradict the H812 no-hit picture and needs immediate audit."
        )
    elif stats["failures_with_single_d3_shield"] == stats["failures_with_monotone_negative"]:
        classification = "h818_single_d3_shields_all_monotone_negatives"
        decision = (
            "Every monotone negative H815 failure is blocked by at least one single translated D3 constraint. "
            "This suggests a local D3 jump-barrier lemma before trying full Plucker-Dodgson identities."
        )
    else:
        classification = "h818_mixed_d3_shield_requires_global_identity"
        decision = (
            "All monotone negatives are blocked by D3 collectively, but not always by a single D3. "
            "This points to a Plucker-Dodgson or multi-D3 determinant identity."
        )
    return {
        "schema": "rh_h818_d3_shield_audit.v0",
        "classification": classification,
        "question": (
            "For H815 four-term failures that are negative under monotonicity alone, are the negative grid "
            "points eliminated by one translated D3 constraint or only by a multi-D3 interaction?"
        ),
        "grid": {
            "parameters": grid["parameters"],
            "stats": stats,
            "examples": grid["negative_records"][: args.example_limit],
        },
        "decision": decision,
        "next_target": {
            "name": "H819 D3 jump-barrier lemma or Plucker-Dodgson",
            "statement": (
                "If single-D3 shields dominate, derive a local jump-barrier inequality from D3_k>=0. "
                "Otherwise search determinant identities combining several D3 constraints."
            ),
        },
    }


def write_markdown(report: dict[str, Any], out_path: Path) -> None:
    lines = [
        "# H818 D3 Shield Audit",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Question",
        "",
        report["question"],
        "",
        "## Stats",
        "",
        "```json",
        json.dumps(report["grid"]["stats"], indent=2),
        "```",
        "",
        "## Examples",
        "",
    ]
    for item in report["grid"]["examples"]:
        lines.extend(["```json", json.dumps(item, indent=2), "```", ""])
    lines.extend(["## Decision", "", report["decision"], "", "## Next Target", ""])
    lines.append(f"`{report['next_target']['name']}`: {report['next_target']['statement']}")
    lines.append("")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="H818 D3 shield audit.")
    parser.add_argument("--length", type=int, default=7)
    parser.add_argument("--denominator", type=int, default=16)
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument("--float-tolerance", type=float, default=1e-14)
    parser.add_argument("--example-limit", type=int, default=8)
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
                "stats": report["grid"]["stats"],
                "next_target": report["next_target"]["name"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
