import argparse
import json
import math
import time
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np


DEFAULT_OUT = "research/riemann/h807_essential_minor_hard_solver_audit.json"


def det3_fraction(matrix: list[list[Fraction]]) -> Fraction:
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def toeplitz_value_fraction(sequence: list[Fraction], row: int, col: int) -> Fraction:
    index = col - row
    if index < 0 or index >= len(sequence):
        return Fraction(0)
    return sequence[index]


def toeplitz_minor_fraction(
    sequence: list[Fraction], rows: tuple[int, ...], cols: tuple[int, ...]
) -> Fraction:
    matrix = [[toeplitz_value_fraction(sequence, row, col) for col in cols] for row in rows]
    return det3_fraction(matrix)


def build_sequence_fraction(q_values: list[Fraction]) -> tuple[list[Fraction], list[Fraction]]:
    ratios = [Fraction(1)]
    for q_value in q_values:
        ratios.append(ratios[-1] * q_value)
    sequence = [Fraction(1)]
    for ratio in ratios:
        sequence.append(sequence[-1] * ratio)
    return sequence, ratios


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def serialize_fraction(value: Fraction) -> dict[str, Any]:
    return {
        "fraction": fraction_text(value),
        "float": float(value),
    }


def all_order3_pairs(sequence_len: int) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    triples = list(combinations(range(sequence_len), 3))
    return [(rows, cols) for rows in triples for cols in triples]


def contiguous_d3_pairs(sequence_len: int) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    return [((0, 1, 2), (start, start + 1, start + 2)) for start in range(0, sequence_len - 2)]


def determinant_terms(rows: tuple[int, ...], cols: tuple[int, ...], sequence_len: int) -> list[tuple[int, int, int, int]]:
    permutations = [
        ((0, 1, 2), 1),
        ((0, 2, 1), -1),
        ((1, 0, 2), -1),
        ((1, 2, 0), 1),
        ((2, 0, 1), 1),
        ((2, 1, 0), -1),
    ]
    terms: list[tuple[int, int, int, int]] = []
    for perm, sign in permutations:
        indexes: list[int] = []
        valid = True
        for row_index, col_index in enumerate(perm):
            source_index = cols[col_index] - rows[row_index]
            if source_index < 0 or source_index >= sequence_len:
                valid = False
                break
            indexes.append(source_index)
        if valid:
            terms.append((sign, indexes[0], indexes[1], indexes[2]))
    return terms


def eval_terms_float(sequences: np.ndarray, terms: list[tuple[int, int, int, int]]) -> np.ndarray:
    values = np.zeros(sequences.shape[0], dtype=np.float64)
    for sign, first, second, third in terms:
        values += sign * sequences[:, first] * sequences[:, second] * sequences[:, third]
    return values


def build_float_sequences(q_grid: np.ndarray) -> np.ndarray:
    batch, length = q_grid.shape
    ratios = np.ones((batch, length + 1), dtype=np.float64)
    if length:
        ratios[:, 1:] = np.cumprod(q_grid, axis=1)
    sequences = np.ones((batch, length + 2), dtype=np.float64)
    sequences[:, 1:] = np.cumprod(ratios, axis=1)
    return sequences


def exact_verify(q_ints: list[int], denominator: int) -> dict[str, Any] | None:
    q_values = [Fraction(value, denominator) for value in q_ints]
    sequence, ratios = build_sequence_fraction(q_values)
    contiguous_rows = []
    for rows, cols in contiguous_d3_pairs(len(sequence)):
        determinant = toeplitz_minor_fraction(sequence, rows, cols)
        contiguous_rows.append((rows, cols, determinant))
        if determinant < 0:
            return None
    negatives = []
    for rows, cols in all_order3_pairs(len(sequence)):
        determinant = toeplitz_minor_fraction(sequence, rows, cols)
        if determinant < 0:
            negatives.append((rows, cols, determinant))
    if not negatives:
        return None
    return {
        "q_values": [serialize_fraction(value) for value in q_values],
        "ratios": [serialize_fraction(value) for value in ratios],
        "sequence": [serialize_fraction(value) for value in sequence],
        "contiguous_d3": [
            {
                "rows": list(rows),
                "cols": list(cols),
                "determinant": serialize_fraction(determinant),
            }
            for rows, cols, determinant in contiguous_rows
        ],
        "first_negative_minor": {
            "rows": list(negatives[0][0]),
            "cols": list(negatives[0][1]),
            "determinant": serialize_fraction(negatives[0][2]),
        },
        "negative_minor_count": len(negatives),
    }



def exact_diagnosis(q_ints: list[int] | None, denominator: int) -> dict[str, Any] | None:
    if q_ints is None:
        return None
    q_values = [Fraction(value, denominator) for value in q_ints]
    sequence, _ = build_sequence_fraction(q_values)
    contiguous_negative = []
    for rows, cols in contiguous_d3_pairs(len(sequence)):
        determinant = toeplitz_minor_fraction(sequence, rows, cols)
        if determinant < 0:
            contiguous_negative.append((rows, cols, determinant))
    min_rows = None
    min_cols = None
    min_det = None
    negative_count = 0
    for rows, cols in all_order3_pairs(len(sequence)):
        determinant = toeplitz_minor_fraction(sequence, rows, cols)
        if min_det is None or determinant < min_det:
            min_rows = rows
            min_cols = cols
            min_det = determinant
        if determinant < 0:
            negative_count += 1
    return {
        "q_values": [serialize_fraction(value) for value in q_values],
        "contiguous_d3_negative_count": len(contiguous_negative),
        "first_contiguous_d3_negative": None
        if not contiguous_negative
        else {
            "rows": list(contiguous_negative[0][0]),
            "cols": list(contiguous_negative[0][1]),
            "determinant": serialize_fraction(contiguous_negative[0][2]),
        },
        "all_order3_negative_count": negative_count,
        "minimum_order3_minor": {
            "rows": [] if min_rows is None else list(min_rows),
            "cols": [] if min_cols is None else list(min_cols),
            "determinant": serialize_fraction(Fraction(0) if min_det is None else min_det),
        },
    }

def hard_search(args: argparse.Namespace) -> dict[str, Any]:
    rng = np.random.default_rng(args.seed)
    deadline = time.monotonic() + args.seconds
    sequence_len = args.length + 2
    all_pairs = all_order3_pairs(sequence_len)
    contiguous_pairs = contiguous_d3_pairs(sequence_len)
    all_terms = [(pair, determinant_terms(pair[0], pair[1], sequence_len)) for pair in all_pairs]
    contiguous_terms = [
        (pair, determinant_terms(pair[0], pair[1], sequence_len)) for pair in contiguous_pairs
    ]
    checked = 0
    contiguous_ok = 0
    best_margin = math.inf
    best_margin_pair: tuple[tuple[int, ...], tuple[int, ...]] | None = None
    best_q: list[int] | None = None
    verified_candidate = None
    while checked < args.attempts and time.monotonic() < deadline:
        batch_size = min(args.batch_size, args.attempts - checked)
        q_int = rng.integers(1, args.denominator + 1, size=(batch_size, args.length))
        q_int.sort(axis=1)
        q_grid = q_int.astype(np.float64) / float(args.denominator)
        sequences = build_float_sequences(q_grid)
        contiguous_mask = np.ones(batch_size, dtype=bool)
        for _, terms in contiguous_terms:
            det_values = eval_terms_float(sequences, terms)
            contiguous_mask &= det_values >= -args.float_tolerance
        contiguous_indices = np.nonzero(contiguous_mask)[0]
        contiguous_ok += int(contiguous_indices.size)
        if contiguous_indices.size:
            subset_sequences = sequences[contiguous_indices]
            subset_q_int = q_int[contiguous_indices]
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
                best_margin_pair = all_pairs[int(min_pair_indices[local_best_index])]
                best_q = [int(value) for value in subset_q_int[local_best_index]]
            negative_indices = np.nonzero(min_values < -args.float_tolerance)[0]
            for candidate_index in negative_indices[: args.exact_candidates_per_batch]:
                candidate = exact_verify(
                    [int(value) for value in subset_q_int[int(candidate_index)]],
                    args.denominator,
                )
                if candidate is not None:
                    verified_candidate = candidate
                    break
        checked += batch_size
        if verified_candidate is not None:
            break
    return {
        "parameters": {
            "length": args.length,
            "denominator": args.denominator,
            "attempts": args.attempts,
            "seconds": args.seconds,
            "batch_size": args.batch_size,
            "seed": args.seed,
            "float_tolerance": args.float_tolerance,
        },
        "stats": {
            "checked": checked,
            "contiguous_d3_nonnegative_float_candidates": contiguous_ok,
            "elapsed_seconds": time.monotonic() - (deadline - args.seconds),
            "best_float_margin": best_margin,
            "best_float_margin_rows": None if best_margin_pair is None else list(best_margin_pair[0]),
            "best_float_margin_cols": None if best_margin_pair is None else list(best_margin_pair[1]),
            "best_q_ints": best_q,
            "best_q_exact_diagnosis": exact_diagnosis(best_q, args.denominator),
        },
        "verified_counterexample": verified_candidate,
    }


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    search = hard_search(args)
    found = search["verified_counterexample"] is not None
    source_audit = [
        {
            "source": "Gasca and Pena / Neville-elimination criterion as summarized in Launois-Lenagan-Rigal",
            "url": "https://arxiv.org/abs/1207.3613",
            "finding": (
                "Finite matrix total positivity can be tested by special initial minors, "
                "but the criterion asks for a broad family of initial minors, not just "
                "PF2, ratio-log-convexity, and translated contiguous D3 minors."
            ),
            "route_impact": "Does not close H807 candidate; points to missing essential minors.",
        },
        {
            "source": "Katkova and Vishnyakova, A sufficient condition for a sequence to be a Pólya frequency sequence",
            "url": "https://arxiv.org/abs/math/0504183",
            "finding": (
                "Strong adjacent dominance conditions can imply PF_m, but for m=3 the "
                "available sufficient condition is far stronger than the Xi tail data, "
                "where q_n tends to 1."
            ),
            "route_impact": "Useful benchmark, too strong for current Xi route.",
        },
        {
            "source": "Pólya/de Bruijn-Newman kernel PF-order warning",
            "url": "https://arxiv.org/abs/2602.20313",
            "finding": (
                "The continuous de Bruijn-Newman kernel is reported not to be PF5. "
                "This warns against continuous-kernel total-positivity shortcuts, but "
                "does not directly address the discrete Xi coefficient sequence."
            ),
            "route_impact": "Prevents overclaiming kernel TP; not an invalidation of JP2'.",
        },
    ]
    if found:
        classification = "pf2_q_contiguous_d3_not_global_pf3_counterexample"
        decision = (
            "The strengthened PF2+q+contiguous-D3 candidate is false. H805/H804 cannot "
            "be promoted to global PF3 without a stronger Xi-specific sparse-minor lemma."
        )
    else:
        classification = "essential_minor_theorem_not_found_hard_counterexample_not_found"
        decision = (
            "H807 did not find a primary theorem proving PF2+q+contiguous-D3 => PF3, "
            "and the faster bounded search did not find a certified hard counterexample. "
            "Treat the candidate as open but unsupported; next progress needs either a "
            "formal order-3 symbolic proof attempt or a stronger nonlinear solver."
        )
    return {
        "schema": "rh_h807_essential_minor_hard_solver_audit.v0",
        "classification": classification,
        "candidate_statement": (
            "Positive Toeplitz sequence + PF2 + nondecreasing q_n=R_n/R_{n-1} + "
            "all translated contiguous order-3 Toeplitz minors nonnegative implies PF3."
        ),
        "source_audit": source_audit,
        "hard_search": search,
        "decision": decision,
        "limitations": [
            "The source audit is a targeted primary-source scout, not a complete literature review.",
            "The hard search is bounded and randomized over rational grids.",
            "No-hit evidence is not proof of the candidate statement.",
            "The candidate statement is about PF3 only, not PF-infinity or RH.",
        ],
        "next_target": {
            "name": "H808 symbolic PF3 cone or nonlinear certificate",
            "statement": (
                "For order 3, symbolically reduce sparse Toeplitz minors under PF2+q+D3, "
                "or use a stronger nonlinear solver to find a certified hard counterexample."
            ),
        },
    }


def write_markdown(report: dict[str, Any], out_path: Path) -> None:
    lines = [
        "# H807 Essential-Minor And Hard-Solver Audit",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Candidate",
        "",
        report["candidate_statement"],
        "",
        "## Source Audit",
        "",
        "| source | finding | impact |",
        "| --- | --- | --- |",
    ]
    for item in report["source_audit"]:
        lines.append(f"| [{item['source']}]({item['url']}) | {item['finding']} | {item['route_impact']} |")
    lines.extend(
        [
            "",
            "## Hard Search",
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
    )
    candidate = report["hard_search"]["verified_counterexample"]
    if candidate is None:
        lines.extend(["", "No certified hard counterexample was found in this bounded run.", ""])
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
    parser = argparse.ArgumentParser(description="H807 essential-minor and hard-solver audit.")
    parser.add_argument("--length", type=int, default=7)
    parser.add_argument("--denominator", type=int, default=512)
    parser.add_argument("--attempts", type=int, default=120000)
    parser.add_argument("--seconds", type=float, default=25.0)
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument("--seed", type=int, default=807)
    parser.add_argument("--float-tolerance", type=float, default=1e-18)
    parser.add_argument("--exact-candidates-per-batch", type=int, default=8)
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
                "hard_stats": report["hard_search"]["stats"],
                "verified_counterexample": report["hard_search"]["verified_counterexample"] is not None,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

