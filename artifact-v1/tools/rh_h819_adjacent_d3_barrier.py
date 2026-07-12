import argparse
import json
import math
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations_with_replacement
from pathlib import Path
from typing import Any

import numpy as np

THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from rh_h813_sparse_row_symbolic_hierarchy import Poly, poly_key, poly_to_string
from rh_h817_four_term_monotone_countercheck import (
    eval_poly_float,
    eval_poly_fraction,
    fraction_text,
    h815_four_term_failures,
)
from rh_h818_d3_shield_audit import d3_polys, run_grid


DEFAULT_OUT = "research/riemann/h819_adjacent_d3_barrier.json"


def batch_iterator(denominator: int, length: int, batch_size: int):
    batch: list[tuple[int, ...]] = []
    for values in combinations_with_replacement(range(1, denominator + 1), length):
        batch.append(values)
        if len(batch) >= batch_size:
            yield np.array(batch, dtype=np.int64)
            batch = []
    if batch:
        yield np.array(batch, dtype=np.int64)


def parse_denominators(raw: str) -> list[int]:
    values = []
    for piece in raw.split(","):
        piece = piece.strip()
        if not piece:
            continue
        values.append(int(piece))
    return values


def q_values_record(q_ints: list[int], denominator: int) -> list[dict[str, Any]]:
    return [
        {
            "name": f"q{index + 2}",
            "fraction": fraction_text(Fraction(value, denominator)),
            "float": value / denominator,
        }
        for index, value in enumerate(q_ints)
    ]


def unique_h818_leftovers(args: argparse.Namespace) -> dict[str, Any]:
    h818_args = argparse.Namespace(
        length=args.length,
        denominator=args.h818_denominator,
        batch_size=args.batch_size,
        float_tolerance=args.float_tolerance,
    )
    h818_grid = run_grid(h818_args)
    failures = h815_four_term_failures(args.length)
    leftovers = [record for record in h818_grid["negative_records"] if not record["has_single_d3_shield"]]

    unique: dict[tuple[tuple[tuple[int, ...], int], ...], dict[str, Any]] = {}
    sources_by_key: dict[tuple[tuple[tuple[int, ...], int], ...], list[dict[str, Any]]] = defaultdict(list)
    for record in leftovers:
        failure = failures[record["failure_index"]]
        key = poly_key(failure["poly"])
        if key not in unique:
            unique[key] = {
                "key": len(unique),
                "poly": failure["poly"],
                "polynomial": record["polynomial"],
                "first_failure_index": record["failure_index"],
            }
        sources_by_key[key].append(
            {
                "failure_index": record["failure_index"],
                "rows": record["rows"],
                "cols": record["cols"],
                "family": record["family"],
                "negative_count": record["negative_count"],
                "single_d3_nonnegative_survivor_counts": record["single_d3_nonnegative_survivor_counts"],
            }
        )

    unique_items = []
    for key, item in unique.items():
        sources = sources_by_key[key]
        active_starts = sorted(
            {
                int(start)
                for source in sources
                for start, count in source["single_d3_nonnegative_survivor_counts"].items()
                if count < source["negative_count"]
            }
        )
        unique_items.append(
            {
                "key": item["key"],
                "poly": item["poly"],
                "polynomial": item["polynomial"],
                "multiplicity": len(sources),
                "active_d3_starts_among_negatives": active_starts,
                "sources": sources,
            }
        )

    return {
        "h818_stats": h818_grid["stats"],
        "leftover_count": len(leftovers),
        "unique_count": len(unique_items),
        "unique_polynomials": unique_items,
        "leftover_family_counts": Counter(record["family"] for record in leftovers).most_common(),
    }


def exact_value_record(poly: Poly, q_ints: list[int], denominator: int) -> dict[str, Any]:
    exact = eval_poly_fraction(q_ints, denominator, poly)
    return {"fraction": fraction_text(exact), "float": float(exact)}


def scan_denominator(
    unique_polynomials: list[dict[str, Any]],
    d3_index: list[dict[str, Any]],
    args: argparse.Namespace,
    denominator: int,
) -> dict[str, Any]:
    total = math.comb(denominator + args.length - 1, args.length)
    checked = 0
    q2_half_candidates = 0
    condition_counts = Counter()
    filters = [
        "monotone_q2_half",
        "d3_2_only",
        "d3_3_only",
        "d3_2_and_d3_3",
        "all_d3",
    ]
    results: dict[str, dict[str, Any]] = {
        str(item["key"]): {
            "polynomial": item["polynomial"],
            "by_filter": {
                name: {
                    "candidate_count": 0,
                    "negative_count": 0,
                    "best_float_value": None,
                    "best_q_ints": None,
                    "best_exact": None,
                    "best_q_values": None,
                }
                for name in filters
            },
        }
        for item in unique_polynomials
    }

    for q_int_batch in batch_iterator(denominator, args.length, args.batch_size):
        checked += int(q_int_batch.shape[0])
        q2_mask = q_int_batch[:, 0] <= denominator // 2
        if not np.any(q2_mask):
            continue
        q_int_batch = q_int_batch[q2_mask]
        q2_half_candidates += int(q_int_batch.shape[0])
        q_values = q_int_batch.astype(np.float64) / float(denominator)
        d3_values = {item["start"]: eval_poly_float(q_values, item["poly"]) for item in d3_index}

        masks = {
            "monotone_q2_half": np.ones(q_int_batch.shape[0], dtype=bool),
            "d3_2_only": d3_values[2] >= -args.float_tolerance,
            "d3_3_only": d3_values[3] >= -args.float_tolerance,
            "d3_2_and_d3_3": (d3_values[2] >= -args.float_tolerance)
            & (d3_values[3] >= -args.float_tolerance),
            "all_d3": np.logical_and.reduce(
                [values >= -args.float_tolerance for values in d3_values.values()]
            ),
        }
        for name, mask in masks.items():
            condition_counts[name] += int(np.count_nonzero(mask))

        for item in unique_polynomials:
            values = eval_poly_float(q_values, item["poly"])
            result = results[str(item["key"])]
            for name, mask in masks.items():
                if not np.any(mask):
                    continue
                slot = result["by_filter"][name]
                slot["candidate_count"] += int(np.count_nonzero(mask))
                masked_values = values[mask]
                negative_mask = masked_values < -args.float_tolerance
                slot["negative_count"] += int(np.count_nonzero(negative_mask))

                masked_indices = np.nonzero(mask)[0]
                local_offset = int(np.argmin(masked_values))
                batch_index = int(masked_indices[local_offset])
                best_value = float(values[batch_index])
                if slot["best_float_value"] is None or best_value < slot["best_float_value"]:
                    q_ints = [int(value) for value in q_int_batch[batch_index]]
                    slot["best_float_value"] = best_value
                    slot["best_q_ints"] = q_ints
                    slot["best_exact"] = exact_value_record(item["poly"], q_ints, denominator)
                    slot["best_q_values"] = q_values_record(q_ints, denominator)

    return {
        "denominator": denominator,
        "grid_total": total,
        "checked": checked,
        "q2_le_half_candidates": q2_half_candidates,
        "condition_counts": dict(condition_counts),
        "polynomial_results": results,
    }


def summarize_scans(scans: list[dict[str, Any]]) -> dict[str, Any]:
    pair_has_negative = False
    d2_only_has_negative = False
    d3_only_has_negative = False
    all_d3_has_negative = False
    worst_pair = None
    for scan in scans:
        for key, result in scan["polynomial_results"].items():
            by_filter = result["by_filter"]
            if by_filter["d3_2_and_d3_3"]["negative_count"]:
                pair_has_negative = True
            if by_filter["d3_2_only"]["negative_count"]:
                d2_only_has_negative = True
            if by_filter["d3_3_only"]["negative_count"]:
                d3_only_has_negative = True
            if by_filter["all_d3"]["negative_count"]:
                all_d3_has_negative = True
            pair_best = by_filter["d3_2_and_d3_3"]["best_float_value"]
            if pair_best is not None and (worst_pair is None or pair_best < worst_pair["best_float_value"]):
                worst_pair = {
                    "denominator": scan["denominator"],
                    "polynomial_key": key,
                    "polynomial": result["polynomial"],
                    "best_float_value": pair_best,
                    "best_exact": by_filter["d3_2_and_d3_3"]["best_exact"],
                    "best_q_values": by_filter["d3_2_and_d3_3"]["best_q_values"],
                }

    if pair_has_negative:
        classification = "h819_adjacent_two_d3_barrier_grid_counterexample"
        decision = (
            "The adjacent pair D3_2,D3_3 is not enough on the tested grids. "
            "Pivot to a wider Plucker-Dodgson identity or Xi-specific input."
        )
    elif d2_only_has_negative and d3_only_has_negative and not all_d3_has_negative:
        classification = "h819_adjacent_two_d3_barrier_supported_on_grids"
        decision = (
            "On the tested grids, every H818 leftover polynomial is nonnegative under "
            "monotonicity, q2<=1/2, and the adjacent constraints D3_2>=0,D3_3>=0. "
            "Each single D3 alone still permits negatives, so the natural lemma is a "
            "two-barrier statement rather than a one-D3 shield."
        )
    else:
        classification = "h819_adjacent_two_d3_barrier_inconclusive"
        decision = (
            "The adjacent pair did not fail, but the single-D3 contrast was not strong "
            "enough on these grids. More diagnostics are needed before promoting a lemma."
        )

    return {
        "classification": classification,
        "decision": decision,
        "pair_has_negative": pair_has_negative,
        "d3_2_only_has_negative": d2_only_has_negative,
        "d3_3_only_has_negative": d3_only_has_negative,
        "all_d3_has_negative": all_d3_has_negative,
        "worst_pair_case": worst_pair,
    }


def compact_scan_table(scans: list[dict[str, Any]]) -> list[dict[str, Any]]:
    table = []
    for scan in scans:
        row = {"denominator": scan["denominator"], "conditions": scan["condition_counts"], "polynomials": []}
        for key, result in scan["polynomial_results"].items():
            row["polynomials"].append(
                {
                    "key": int(key),
                    "polynomial": result["polynomial"],
                    "negative_counts": {
                        name: slot["negative_count"] for name, slot in result["by_filter"].items()
                    },
                    "best_d3_2_and_d3_3": result["by_filter"]["d3_2_and_d3_3"]["best_exact"],
                }
            )
        table.append(row)
    return table


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    leftovers = unique_h818_leftovers(args)
    d3_index = d3_polys(args.length)
    denominators = parse_denominators(args.denominators)
    scans = [scan_denominator(leftovers["unique_polynomials"], d3_index, args, denom) for denom in denominators]
    summary = summarize_scans(scans)
    return {
        "schema": "rh_h819_adjacent_d3_barrier.v0",
        "classification": summary["classification"],
        "question": (
            "Do the 20 H818 no-single-shield four-term failures reduce to an adjacent "
            "two-D3 barrier using only D3_start_2 and D3_start_3?"
        ),
        "leftovers": {
            "h818_stats": leftovers["h818_stats"],
            "leftover_count": leftovers["leftover_count"],
            "unique_count": leftovers["unique_count"],
            "leftover_family_counts": leftovers["leftover_family_counts"],
            "unique_polynomials": [
                {key: value for key, value in item.items() if key != "poly"}
                for item in leftovers["unique_polynomials"]
            ],
        },
        "scan_summary": summary,
        "scan_table": compact_scan_table(scans),
        "full_scans": scans,
        "underexploited_next_ideas": [
            {
                "name": "Adjacent two-D3 jump-barrier lemma",
                "status": "primary",
                "reason": "H818 leftovers only need starts 2 and 3 on the grid; H816 did not test this scalar two-barrier form.",
            },
            {
                "name": "Dodgson condensation identity for sparse rows",
                "status": "secondary",
                "reason": "If the scalar barrier proof stalls, express sparse minors through adjacent contiguous minors and PF2 factors.",
            },
            {
                "name": "Xi-specific unit-deficit input",
                "status": "fallback",
                "reason": "H817 showed monotone q alone is false; an analytic input from Xi coefficients may avoid proving a general q-cone theorem.",
            },
            {
                "name": "Cross-route Nyman-Beurling endpoint bridge",
                "status": "separate branch",
                "reason": "Older H520/H662 routes still contain formal candidates, especially endpoint Abel cancellation, not used by the PF3 corridor.",
            },
        ],
        "decision": summary["decision"],
        "next_target": {
            "name": "H820 formal adjacent D3 lemma",
            "statement": (
                "Try to prove the six unique normalized four-term forms from monotonicity, q2<=1/2, "
                "D3_2>=0, and D3_3>=0, or produce a continuous counterexample."
            ),
        },
    }


def write_markdown(report: dict[str, Any], out_path: Path) -> None:
    lines = [
        "# H819 Adjacent D3 Barrier",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Question",
        "",
        report["question"],
        "",
        "## H818 Leftovers",
        "",
        "```json",
        json.dumps(
            {
                "h818_stats": report["leftovers"]["h818_stats"],
                "leftover_count": report["leftovers"]["leftover_count"],
                "unique_count": report["leftovers"]["unique_count"],
                "leftover_family_counts": report["leftovers"]["leftover_family_counts"],
            },
            indent=2,
        ),
        "```",
        "",
        "## Unique Polynomials",
        "",
    ]
    for item in report["leftovers"]["unique_polynomials"]:
        lines.extend(
            [
                f"- P{item['key']}: `{item['polynomial']}`",
                f"  - multiplicity: `{item['multiplicity']}`",
                f"  - active D3 starts among negatives: `{item['active_d3_starts_among_negatives']}`",
            ]
        )
    lines.extend(
        [
            "",
            "## Scan Summary",
            "",
            "```json",
            json.dumps(report["scan_summary"], indent=2),
            "```",
            "",
            "## Compact Scan Table",
            "",
            "```json",
            json.dumps(report["scan_table"], indent=2),
            "```",
            "",
            "## Underexploited Next Ideas",
            "",
        ]
    )
    for item in report["underexploited_next_ideas"]:
        lines.append(f"- `{item['name']}` ({item['status']}): {item['reason']}")
    lines.extend(["", "## Decision", "", report["decision"], "", "## Next Target", ""])
    lines.append(f"`{report['next_target']['name']}`: {report['next_target']['statement']}")
    lines.append("")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="H819 adjacent two-D3 barrier probe.")
    parser.add_argument("--length", type=int, default=7)
    parser.add_argument("--h818-denominator", type=int, default=16)
    parser.add_argument("--denominators", default="16,24")
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument("--float-tolerance", type=float, default=1e-14)
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
                "leftover_count": report["leftovers"]["leftover_count"],
                "unique_count": report["leftovers"]["unique_count"],
                "scan_summary": report["scan_summary"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
