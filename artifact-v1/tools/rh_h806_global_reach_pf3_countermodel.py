import argparse
import json
import random
import time
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any


DEFAULT_OUT = "research/riemann/h806_global_reach_pf3_countermodel.json"


def det3(matrix: list[list[Fraction]]) -> Fraction:
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def toeplitz_value(sequence: list[Fraction], row: int, col: int) -> Fraction:
    index = col - row
    if index < 0 or index >= len(sequence):
        return Fraction(0)
    return sequence[index]


def toeplitz_minor(sequence: list[Fraction], rows: tuple[int, ...], cols: tuple[int, ...]) -> Fraction:
    matrix = [[toeplitz_value(sequence, row, col) for col in cols] for row in rows]
    return det3(matrix)


def build_sequence(q_values: list[Fraction], first_ratio: Fraction = Fraction(1)) -> tuple[list[Fraction], list[Fraction]]:
    ratios = [Fraction(0), first_ratio]
    for q_value in q_values:
        ratios.append(ratios[-1] * q_value)
    sequence = [Fraction(1)]
    for index in range(1, len(ratios)):
        sequence.append(sequence[-1] * ratios[index])
    return sequence, ratios[1:]


def contiguous_d3_rows(sequence: list[Fraction]) -> list[dict[str, Any]]:
    max_index = len(sequence) - 1
    rows = []
    for start in range(0, max_index - 1):
        row_tuple = (0, 1, 2)
        col_tuple = (start, start + 1, start + 2)
        determinant = toeplitz_minor(sequence, row_tuple, col_tuple)
        rows.append(
            {
                "rows": row_tuple,
                "cols": col_tuple,
                "determinant": determinant,
                "nonnegative": determinant >= 0,
            }
        )
    return rows


def negative_order3_minors(sequence: list[Fraction]) -> list[dict[str, Any]]:
    max_index = len(sequence) - 1
    negatives = []
    for rows in combinations(range(max_index + 1), 3):
        for cols in combinations(range(max_index + 1), 3):
            determinant = toeplitz_minor(sequence, rows, cols)
            if determinant < 0:
                negatives.append(
                    {
                        "rows": rows,
                        "cols": cols,
                        "determinant": determinant,
                    }
                )
    return negatives


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def fraction_float(value: Fraction) -> float:
    return float(value.numerator / value.denominator)


def serialize_fraction(value: Fraction) -> dict[str, Any]:
    return {
        "fraction": fraction_text(value),
        "float": fraction_float(value),
    }


def serialize_minor(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "rows": list(row["rows"]),
        "cols": list(row["cols"]),
        "determinant": serialize_fraction(row["determinant"]),
    }


def serialize_candidate(
    q_values: list[Fraction],
    sequence: list[Fraction],
    ratios: list[Fraction],
    negatives: list[dict[str, Any]],
    contiguous_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    return {
        "q_values": [serialize_fraction(value) for value in q_values],
        "ratios": [serialize_fraction(value) for value in ratios],
        "sequence": [serialize_fraction(value) for value in sequence],
        "first_negative_minor": serialize_minor(negatives[0]),
        "negative_minor_count": len(negatives),
        "contiguous_d3": [
            {
                "rows": list(row["rows"]),
                "cols": list(row["cols"]),
                "determinant": serialize_fraction(row["determinant"]),
                "nonnegative": row["nonnegative"],
            }
            for row in contiguous_rows
        ],
    }


def is_q_monotone(q_values: list[Fraction]) -> bool:
    return all(q_values[index] <= q_values[index + 1] for index in range(len(q_values) - 1))


def deterministic_weak_counterexample() -> dict[str, Any]:
    q_values = [Fraction(1, 12), Fraction(1, 12), Fraction(1, 12), Fraction(1, 12), Fraction(2, 3), Fraction(2, 3)]
    sequence, ratios = build_sequence(q_values)
    contiguous_rows = contiguous_d3_rows(sequence)
    negatives = negative_order3_minors(sequence)
    return serialize_candidate(q_values, sequence, ratios, negatives, contiguous_rows)


def random_q_values(length: int, denominator: int, rng: random.Random) -> list[Fraction]:
    values = [Fraction(rng.randint(1, denominator), denominator) for _ in range(length)]
    values.sort()
    return values


def search_hard_counterexample(
    *,
    length: int,
    denominator: int,
    attempts: int,
    seconds: float,
    seed: int,
) -> tuple[dict[str, Any] | None, dict[str, Any]]:
    rng = random.Random(seed)
    deadline = time.monotonic() + seconds
    checked = 0
    contiguous_ok = 0
    while checked < attempts and time.monotonic() < deadline:
        checked += 1
        q_values = random_q_values(length, denominator, rng)
        if not is_q_monotone(q_values):
            raise AssertionError("internal q generator failed")
        sequence, ratios = build_sequence(q_values)
        contiguous_rows = contiguous_d3_rows(sequence)
        if not all(row["nonnegative"] for row in contiguous_rows):
            continue
        contiguous_ok += 1
        negatives = negative_order3_minors(sequence)
        if negatives:
            return (
                serialize_candidate(q_values, sequence, ratios, negatives, contiguous_rows),
                {
                    "checked": checked,
                    "contiguous_d3_nonnegative_candidates": contiguous_ok,
                    "elapsed_seconds": time.monotonic() - (deadline - seconds),
                },
            )
    return (
        None,
        {
            "checked": checked,
            "contiguous_d3_nonnegative_candidates": contiguous_ok,
            "elapsed_seconds": time.monotonic() - (deadline - seconds),
        },
    )


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    weak = deterministic_weak_counterexample()
    hard, stats = search_hard_counterexample(
        length=args.length,
        denominator=args.denominator,
        attempts=args.attempts,
        seconds=args.seconds,
        seed=args.seed,
    )
    if hard is not None:
        classification = "pf2_q_contiguous_d3_does_not_imply_pf3_counterexample"
        decision = (
            "The current ratio conditions plus contiguous order-3 positivity still do not "
            "force PF3. H805/H804 cannot be a global JP2' route without additional "
            "Xi-specific sparse-minor structure or an essential-minor theorem."
        )
    else:
        classification = "pf2_q_contiguous_d3_counterexample_not_found_bounded_search"
        decision = (
            "A weak PF2+q counterexample exists, but the bounded hard search did not find "
            "a counterexample that also keeps all translated contiguous D3 minors nonnegative. "
            "This is not proof; it suggests a possible structural lemma worth formalizing or "
            "testing with a stronger solver."
        )
    return {
        "schema": "rh_h806_global_reach_pf3_countermodel.v0",
        "classification": classification,
        "definitions": {
            "PF2": "positive a_n with nonincreasing ratios R_n=a_n/a_{n-1}",
            "q_monotone": "q_n=R_n/R_{n-1} is nondecreasing",
            "contiguous_D3": "order-3 Toeplitz minors with rows (0,1,2), cols (k,k+1,k+2)",
            "PF3": "all order-3 Toeplitz minors are nonnegative",
        },
        "weak_counterexample": {
            "meaning": "PF2 + q_monotone alone does not imply PF3.",
            "candidate": weak,
        },
        "hard_search": {
            "meaning": "Search for PF2 + q_monotone + contiguous_D3_nonnegative but not PF3.",
            "parameters": {
                "length": args.length,
                "denominator": args.denominator,
                "attempts": args.attempts,
                "seconds": args.seconds,
                "seed": args.seed,
            },
            "stats": stats,
            "candidate": hard,
        },
        "decision": decision,
        "limitations": [
            "The hard search is bounded random evidence, not a theorem.",
            "The weak counterexample already blocks PF2+q_monotone as a global route.",
            "Failure to find the hard counterexample may reflect the search distribution.",
            "No statement here invalidates JP2' for the Xi coefficients.",
        ],
        "next_target": {
            "name": "H807 essential-minor or hard solver audit",
            "statement": (
                "Either prove a PF2+q+contiguous-D3-to-PF3 theorem from total-positivity "
                "literature, or use SAT/SMT/nonlinear search to find a hard counterexample."
            ),
        },
    }


def write_markdown(report: dict[str, Any], out_path: Path) -> None:
    lines = [
        "# H806 Global Reach PF3 Countermodel Audit",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Question",
        "",
        "Does the current H804/H805-style ratio control have global PF3 force, or only local force for selected sparse families?",
        "",
        "## Weak Counterexample",
        "",
        report["weak_counterexample"]["meaning"],
        "",
        "```json",
        json.dumps(report["weak_counterexample"]["candidate"], indent=2),
        "```",
        "",
        "## Hard Search",
        "",
        report["hard_search"]["meaning"],
        "",
        "Parameters:",
        "",
        "```json",
        json.dumps(report["hard_search"]["parameters"], indent=2),
        "```",
        "",
        "Stats:",
        "",
        "```json",
        json.dumps(report["hard_search"]["stats"], indent=2),
        "```",
    ]
    if report["hard_search"]["candidate"] is None:
        lines.extend(["", "No hard counterexample found in the bounded search.", ""])
    else:
        lines.extend(
            [
                "",
                "Hard counterexample:",
                "",
                "```json",
                json.dumps(report["hard_search"]["candidate"], indent=2),
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
    parser = argparse.ArgumentParser(description="H806 PF2/q-monotone global reach audit for PF3.")
    parser.add_argument("--length", type=int, default=7)
    parser.add_argument("--denominator", type=int, default=64)
    parser.add_argument("--attempts", type=int, default=25000)
    parser.add_argument("--seconds", type=float, default=20.0)
    parser.add_argument("--seed", type=int, default=806)
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
                "weak_negative_minor": report["weak_counterexample"]["candidate"]["first_negative_minor"],
                "hard_stats": report["hard_search"]["stats"],
                "hard_found": report["hard_search"]["candidate"] is not None,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
