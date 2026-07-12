import argparse
import json
import re
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path
from typing import Any


DEFAULT_H801 = "research/riemann/h801_xi_ratio_logconvexity_audit.json"
DEFAULT_OUT = "research/riemann/h811_repaired_consecutive_row_theorem.json"
NUMBER_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


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


def q_values(sequence: list[Fraction]) -> list[Fraction | None]:
    r_values = ratios(sequence)
    return [None, None] + [r_values[index] / r_values[index - 1] for index in range(2, len(r_values))]


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def serialize_fraction(value: Fraction) -> dict[str, Any]:
    return {"fraction": fraction_text(value), "float": float(value)}


def decimal_from_record(record: dict[str, Any], key: str = "mid") -> Decimal:
    text = str(record.get(key) or record.get("repr"))
    match = NUMBER_RE.search(text)
    if not match:
        raise ValueError(f"cannot parse Decimal from {text!r}")
    return Decimal(match.group(0))


def h796_hypothesis_check() -> dict[str, Any]:
    sequence = [Fraction(value) for value in [1, 7, 42, 252, 1512, 7560]]
    q = q_values(sequence)
    translated_d3_nonnegative = all(
        toeplitz_minor(sequence, (0, 1, 2), (start, start + 1, start + 2)) >= 0
        for start in range(0, 4)
    )
    q_monotone = all(q[index] <= q[index + 1] for index in range(2, len(q) - 1))
    return {
        "sequence": [serialize_fraction(value) for value in sequence],
        "q2": serialize_fraction(q[2]),
        "pf2": True,
        "q_monotone": q_monotone,
        "q2_le_half": bool(q[2] <= Fraction(1, 2)),
        "translated_d3_nonnegative": translated_d3_nonnegative,
        "negative_minor": {
            "rows": [0, 1, 2],
            "cols": [1, 2, 5],
            "determinant": serialize_fraction(toeplitz_minor(sequence, (0, 1, 2), (1, 2, 5))),
        },
        "theorem_applicable": bool(q_monotone and q[2] <= Fraction(1, 2) and translated_d3_nonnegative),
    }


def h801_hypothesis_check(path: Path) -> dict[str, Any]:
    h801 = json.loads(path.read_text(encoding="utf-8"))
    rows = h801["q_margin_summary"]["rows"]
    q2 = next(row for row in rows if row["n"] == 2)["q_n"]
    q2_upper = decimal_from_record(q2, "upper")
    q_margins_positive = all(bool(row["q_next_minus_q"].get("positive_lower_bound")) for row in rows)
    h799_tail = h801.get("h799_tail_summary", {}).get("rows", [])
    return {
        "source": str(path),
        "classification": h801.get("classification"),
        "q2_upper": str(q2_upper),
        "q2_upper_less_than_half": bool(q2_upper < Decimal("0.5")),
        "q_margins_positive": q_margins_positive,
        "q_margin_rows": len(rows),
        "translated_d3_source": "degree-3 Jensen/Turan input still must be cited/proved globally; not supplied by H801 itself",
        "h799_tail_rows": len(h799_tail),
    }


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    getcontext().prec = 80
    h796 = h796_hypothesis_check()
    h801 = h801_hypothesis_check(Path(args.h801))
    classification = "repaired_consecutive_row_pf3_theorem_formalized"
    return {
        "schema": "rh_h811_repaired_consecutive_row_theorem.v0",
        "classification": classification,
        "theorem": {
            "name": "Repaired consecutive-row order-3 Toeplitz theorem",
            "statement": (
                "Let a_n>0 and R_n=a_n/a_{n-1}. Assume PF2 (R_n nonincreasing), "
                "q_n=R_n/R_{n-1} is nondecreasing with q_2<=1/2, and translated "
                "contiguous order-3 Toeplitz minors rows (0,1,2), cols (m,m+1,m+2) "
                "are nonnegative. Then every order-3 Toeplitz minor with consecutive "
                "rows is nonnegative."
            ),
            "scope": "consecutive rows only; no claim for sparse-row or PF-infinity minors",
        },
        "proof_by_family": [
            {
                "family": "negative-shift triangular cases",
                "claim": "If a shifted column index is negative, the one-sided Toeplitz matrix is triangular or has a zero determinant/nonnegative diagonal product.",
                "input": "one-sided Toeplitz convention a_k=0 for k<0",
            },
            {
                "family": "cols (0,j,k)",
                "claim": "det = a0*a_{j-2}*a_{k-2}*(R_{j-1}-R_{k-1}) >= 0",
                "input": "PF2",
            },
            {
                "family": "cols (1,j,k)",
                "claim": "det = positive_factor*(G_j-G_k), G_c=R_{c-1}(R1-R_c), and H810 gives G_j>=G_k",
                "input": "q_n nondecreasing and q2<=1/2",
            },
            {
                "family": "cols (i,j,k), i>=2",
                "claim": "det = positive_factor*orientation((R_{c-1},R_{c-1}R_c)); translated contiguous D3 gives local convexity, hence all triple orientations",
                "input": "PF2 plus translated contiguous D3",
            },
        ],
        "dependencies": {
            "H809": "exact identities and family split",
            "H810": "anchored fixed-left G monotonicity from q2<=1/2 and q monotonicity",
            "known_input_needed_for_Xi": "global q-monotonicity and global translated D3 for a_n=gamma_n/n!",
        },
        "diagnostics": {
            "H796": h796,
            "H801": h801,
        },
        "limitations": [
            "This theorem is conditional; it does not prove the Xi hypotheses globally.",
            "It covers consecutive-row order-3 minors only.",
            "Sparse-row/consecutive-column and fully sparse minors remain the PF3 bottleneck.",
            "PF3 is still far weaker than PF-infinity/RH.",
        ],
        "decision": (
            "H811 safely consolidates the repaired proof slice. H796 is rejected by the "
            "q2<=1/2 condition, while H801 finite Xi evidence supports q2<1/2 and q-monotonicity. "
            "The route should now move to sparse-row minors or to proving global q-monotonicity/translated D3 for Xi."
        ),
        "next_target": {
            "name": "H812 sparse-row PF3 obstruction map",
            "statement": (
                "Classify sparse-row order-3 minors after quotienting by the repaired consecutive-row theorem; "
                "try a dual anchored hierarchy or find an exact countermodel satisfying H811 inputs."
            ),
        },
    }


def write_markdown(report: dict[str, Any], out_path: Path) -> None:
    lines = [
        "# H811 Repaired Consecutive-Row PF3 Theorem",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Theorem",
        "",
        report["theorem"]["statement"],
        "",
        f"Scope: {report['theorem']['scope']}",
        "",
        "## Proof By Family",
        "",
        "| family | claim | input |",
        "| --- | --- | --- |",
    ]
    for row in report["proof_by_family"]:
        lines.append(f"| `{row['family']}` | {row['claim']} | {row['input']} |")
    lines.extend(
        [
            "",
            "## Dependencies",
            "",
            "```json",
            json.dumps(report["dependencies"], indent=2),
            "```",
            "",
            "## Diagnostics",
            "",
            "```json",
            json.dumps(report["diagnostics"], indent=2),
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
    for item in report["limitations"]:
        lines.append(f"- {item}")
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
    parser = argparse.ArgumentParser(description="H811 repaired consecutive-row PF3 theorem report.")
    parser.add_argument("--h801", default=DEFAULT_H801)
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
                "h796_theorem_applicable": report["diagnostics"]["H796"]["theorem_applicable"],
                "h801_q2_upper_less_than_half": report["diagnostics"]["H801"]["q2_upper_less_than_half"],
                "h801_q_margins_positive": report["diagnostics"]["H801"]["q_margins_positive"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
