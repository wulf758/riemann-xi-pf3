import argparse
import json
import random
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any


DEFAULT_OUT = "research/riemann/h809_consecutive_row_pf3_decomposition.json"


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


def toeplitz_minor(sequence: list[Fraction], rows: tuple[int, int, int], cols: tuple[int, int, int]) -> Fraction:
    return det3([[toeplitz_value(sequence, row, col) for col in cols] for row in rows])


def ratios(sequence: list[Fraction]) -> list[Fraction | None]:
    return [None] + [sequence[index] / sequence[index - 1] for index in range(1, len(sequence))]


def random_log_concave_sequence(length: int, rng: random.Random) -> list[Fraction]:
    q_values = sorted(Fraction(rng.randint(1, 11), 12) for _ in range(length - 2))
    ratio_values = [Fraction(1)]
    for q_value in q_values:
        ratio_values.append(ratio_values[-1] * q_value)
    sequence = [Fraction(1)]
    for ratio_value in ratio_values:
        sequence.append(sequence[-1] * ratio_value)
    return sequence


def column_point(sequence: list[Fraction], col: int) -> tuple[Fraction, Fraction]:
    r_values = ratios(sequence)
    return (r_values[col - 1], r_values[col - 1] * r_values[col])


def orientation(points: list[tuple[Fraction, Fraction]]) -> Fraction:
    matrix = [[point[1] for point in points], [point[0] for point in points], [Fraction(1), Fraction(1), Fraction(1)]]
    return det3(matrix)


def formula_col0(sequence: list[Fraction], j: int, k: int) -> Fraction:
    r_values = ratios(sequence)
    return sequence[0] * sequence[j - 2] * sequence[k - 2] * (r_values[j - 1] - r_values[k - 1])


def g_value(sequence: list[Fraction], col: int) -> Fraction:
    r_values = ratios(sequence)
    return r_values[col - 1] * (r_values[1] - r_values[col])


def formula_col1(sequence: list[Fraction], j: int, k: int) -> Fraction:
    return sequence[0] * sequence[j - 2] * sequence[k - 2] * (g_value(sequence, j) - g_value(sequence, k))


def formula_col_ge2(sequence: list[Fraction], i: int, j: int, k: int) -> Fraction:
    points = [column_point(sequence, col) for col in (i, j, k)]
    return sequence[i - 2] * sequence[j - 2] * sequence[k - 2] * orientation(points)


def contiguous_d3(sequence: list[Fraction], start: int) -> Fraction:
    return toeplitz_minor(sequence, (0, 1, 2), (start, start + 1, start + 2))


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def serialize_fraction(value: Fraction) -> dict[str, Any]:
    return {"fraction": fraction_text(value), "float": float(value)}


def verify_identities(max_col: int, trials: int, seed: int) -> dict[str, Any]:
    rng = random.Random(seed)
    failures: list[dict[str, Any]] = []
    for _ in range(trials):
        sequence = random_log_concave_sequence(max_col + 1, rng)
        for i, j, k in combinations(range(0, max_col + 1), 3):
            direct = toeplitz_minor(sequence, (0, 1, 2), (i, j, k))
            if i == 0 and j >= 2:
                formula = formula_col0(sequence, j, k)
                family = "first_col_0"
            elif i == 1 and j >= 2:
                formula = formula_col1(sequence, j, k)
                family = "first_col_1"
            elif i >= 2:
                formula = formula_col_ge2(sequence, i, j, k)
                family = "first_col_ge2"
            else:
                # cols (0,1,k) is triangular and directly positive.
                formula = direct
                family = "triangular_01"
            if direct != formula:
                failures.append(
                    {
                        "family": family,
                        "cols": [i, j, k],
                        "direct": fraction_text(direct),
                        "formula": fraction_text(formula),
                    }
                )
                if len(failures) >= 5:
                    break
        if len(failures) >= 5:
            break
    return {"trials": trials, "max_col": max_col, "failure_count": len(failures), "failures": failures}


def classify_order3_shapes(max_index: int) -> dict[str, Any]:
    counts = {
        "consecutive_rows_total": 0,
        "consecutive_rows_covered": 0,
        "consecutive_rows_zero_by_negative_shift": 0,
        "consecutive_rows_formula_covered": 0,
        "other_rows_uncovered": 0,
        "total_order3_minors": 0,
    }
    examples: dict[str, Any] = {}
    triples = list(combinations(range(max_index + 1), 3))
    for rows in triples:
        for cols in triples:
            counts["total_order3_minors"] += 1
            consecutive_rows = rows[1] == rows[0] + 1 and rows[2] == rows[1] + 1
            if consecutive_rows:
                counts["consecutive_rows_total"] += 1
                shifted_cols = tuple(col - rows[0] for col in cols)
                counts["consecutive_rows_covered"] += 1
                if shifted_cols[0] < 0:
                    counts["consecutive_rows_zero_by_negative_shift"] += 1
                    examples.setdefault("covered_zero_shift", {"rows": list(rows), "cols": list(cols)})
                else:
                    counts["consecutive_rows_formula_covered"] += 1
                    examples.setdefault("covered_formula", {"rows": list(rows), "cols": list(cols)})
            else:
                counts["other_rows_uncovered"] += 1
                examples.setdefault("uncovered_other_rows", {"rows": list(rows), "cols": list(cols)})
    return {"max_index": max_index, "counts": counts, "examples": examples}


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    identity_checks = verify_identities(args.max_col, args.trials, args.seed)
    shape_coverage = classify_order3_shapes(args.max_index)
    classification = (
        "consecutive_row_identity_decomposition_verified_but_pf2_d3_proof_invalidated"
        if identity_checks["failure_count"] == 0
        else "consecutive_row_order3_symbolic_decomposition_identity_failure"
    )
    decision = (
        "All tested exact identities passed, but the proposed proof from PF2 plus "
        "translated contiguous D3 is invalidated by H796. The first-column-1 family "
        "needs monotonicity of G_c=R_{c-1}(R1-R_c); adjacent differences of G_c are "
        "fixed-left minors with columns (1,c,c+1), not the translated contiguous D3 "
        "family (c-1,c,c+1). H809 is therefore an identity decomposition and a map of "
        "the missing anchored-minor hierarchy, not a proof slice."
    )
    return {
        "schema": "rh_h809_consecutive_row_pf3_decomposition.v0",
        "classification": classification,
        "statement": (
            "For one-sided Toeplitz order-3 minors with consecutive rows, the minors "
            "split into exact ratio/geometry identities. The earlier claim that PF2 "
            "plus translated contiguous D3 nonnegativity implies all sparse-column "
            "minors is false; extra anchored minors are needed."
        ),
        "normalization": "Signs are invariant under a_n -> a0 * R1^n * b_n, so formulas may set a0=R1=1 for intuition.",
        "families": [
            {
                "name": "first_col_0",
                "cols": "(0,j,k), 2<=j<k",
                "formula": "det = a0*a_{j-2}*a_{k-2}*(R_{j-1}-R_{k-1})",
                "nonnegative_source": "PF2 gives R_{j-1}>=R_{k-1}.",
            },
            {
                "name": "first_col_1",
                "cols": "(1,j,k), 2<=j<k",
                "formula": "det = a0*a_{j-2}*a_{k-2}*(G_j-G_k), G_c=R_{c-1}(R1-R_c)",
                "nonnegative_source": "Requires anchored fixed-left minors with cols (1,c,c+1); translated contiguous D3 alone is not enough.",
            },
            {
                "name": "first_col_ge2",
                "cols": "(i,j,k), 2<=i<j<k",
                "formula": "det = a_{i-2}*a_{j-2}*a_{k-2}*orient((X_i,Y_i),(X_j,Y_j),(X_k,Y_k))",
                "nonnegative_source": "X_c=R_{c-1} is monotone by PF2; contiguous D3 gives local convexity, hence all triple orientations.",
            },
        ],
        "identity_checks": identity_checks,
        "invalidating_counterexample": {
            "source": "H796",
            "sequence": [1, 7, 42, 252, 1512, 7560],
            "ratios": [7, 6, 6, 6, 5],
            "translated_contiguous_d3": {
                "D0": "1",
                "D1": "7",
                "D2": "0",
                "D3": "0",
            },
            "negative_sparse_consecutive_row_minor": {
                "rows": [0, 1, 2],
                "cols": [1, 2, 5],
                "determinant": "-1260",
            },
            "lesson": "G-monotonicity is controlled by anchored minors (1,c,c+1), not by translated contiguous D3.",
        },        "shape_coverage": shape_coverage,
        "limitations": [
            "This is an identity decomposition, not a proof from the H807 assumptions.",
            "PF2 plus translated contiguous D3 is insufficient, by H796.",
            "A viable route needs an anchored-minor hierarchy such as fixed-left minors (1,c,c+1), or a different Xi-specific constraint.",
        ],
        "decision": decision,
        "next_target": {
            "name": "H810 anchored-minor hierarchy audit",
            "statement": (
                "Try to transform sparse-row/consecutive-column and fully sparse order-3 "
                "Toeplitz minors into a dual convexity or monotone-slope problem."
            ),
        },
    }


def write_markdown(report: dict[str, Any], out_path: Path) -> None:
    lines = [
        "# H809 Consecutive-Row PF3 Decomposition",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Statement",
        "",
        report["statement"],
        "",
        "## Families",
        "",
        "| family | columns | formula | source of nonnegativity |",
        "| --- | --- | --- | --- |",
    ]
    for family in report["families"]:
        lines.append(
            f"| `{family['name']}` | `{family['cols']}` | `{family['formula']}` | {family['nonnegative_source']} |"
        )
    lines.extend(
        [
            "",
            "## Identity Checks",
            "",
            "```json",
            json.dumps(report["identity_checks"], indent=2),
            "```",
            "",
            "## Invalidating Counterexample",
            "",
            "```json",
            json.dumps(report["invalidating_counterexample"], indent=2),
            "```",
            "",
            "## Shape Coverage",
            "",
            "```json",
            json.dumps(report["shape_coverage"], indent=2),
            "```",
            "",
            "## Decision",
            "",
            report["decision"],
            "",
            "## Limitations",
            "",
        ]
    )
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
    parser = argparse.ArgumentParser(description="H809 symbolic decomposition for consecutive-row order-3 Toeplitz minors.")
    parser.add_argument("--max-col", type=int, default=10)
    parser.add_argument("--max-index", type=int, default=8)
    parser.add_argument("--trials", type=int, default=100)
    parser.add_argument("--seed", type=int, default=809)
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
                "identity_failure_count": report["identity_checks"]["failure_count"],
                "covered_consecutive_rows": report["shape_coverage"]["counts"]["consecutive_rows_covered"],
                "other_rows_uncovered": report["shape_coverage"]["counts"]["other_rows_uncovered"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()



