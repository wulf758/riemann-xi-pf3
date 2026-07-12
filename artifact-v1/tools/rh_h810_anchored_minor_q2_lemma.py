import argparse
import json
import re
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path
from typing import Any


DEFAULT_H801 = "research/riemann/h801_xi_ratio_logconvexity_audit.json"
DEFAULT_OUT = "research/riemann/h810_anchored_minor_q2_lemma.json"
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


def q_values_from_ratios(r_values: list[Fraction | None]) -> list[Fraction | None]:
    return [None, None] + [
        r_values[index] / r_values[index - 1] for index in range(2, len(r_values))
    ]


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


def h796_diagnosis() -> dict[str, Any]:
    sequence = [Fraction(value) for value in [1, 7, 42, 252, 1512, 7560]]
    r_values = ratios(sequence)
    q_values = q_values_from_ratios(r_values)
    fixed_left = []
    for c in range(2, 5):
        fixed_left.append(
            {
                "c": c,
                "cols": [1, c, c + 1],
                "determinant": serialize_fraction(toeplitz_minor(sequence, (0, 1, 2), (1, c, c + 1))),
            }
        )
    translated = []
    for start in range(0, 4):
        translated.append(
            {
                "start": start,
                "cols": [start, start + 1, start + 2],
                "determinant": serialize_fraction(
                    toeplitz_minor(sequence, (0, 1, 2), (start, start + 1, start + 2))
                ),
            }
        )
    return {
        "sequence": [serialize_fraction(value) for value in sequence],
        "ratios": [None if value is None else serialize_fraction(value) for value in r_values],
        "q_values": [None if value is None else serialize_fraction(value) for value in q_values],
        "q2_le_half": bool(q_values[2] <= Fraction(1, 2)),
        "translated_contiguous_d3": translated,
        "fixed_left_anchored_minors": fixed_left,
        "negative_sparse_minor": {
            "rows": [0, 1, 2],
            "cols": [1, 2, 5],
            "determinant": serialize_fraction(toeplitz_minor(sequence, (0, 1, 2), (1, 2, 5))),
        },
    }


def h801_diagnosis(path: Path) -> dict[str, Any]:
    h801 = json.loads(path.read_text(encoding="utf-8"))
    rows = h801["q_margin_summary"]["rows"]
    q2_row = next(row for row in rows if row["n"] == 2)
    q2_upper = decimal_from_record(q2_row["q_n"], "upper")
    all_q_margins_positive = all(bool(row["q_next_minus_q"].get("positive_lower_bound")) for row in rows)
    return {
        "source": str(path),
        "classification": h801.get("classification"),
        "q2_upper": str(q2_upper),
        "q2_upper_less_than_half": bool(q2_upper < Decimal("0.5")),
        "q_margin_rows": len(rows),
        "all_q_next_minus_q_certified_positive": all_q_margins_positive,
        "q2_record": q2_row["q_n"],
    }


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    getcontext().prec = 80
    h796 = h796_diagnosis()
    h801 = h801_diagnosis(Path(args.h801))
    classification = (
        "anchored_minor_q2_lemma_supported_and_repairs_h809_gap"
        if h801["q2_upper_less_than_half"] and h801["all_q_next_minus_q_certified_positive"]
        else "anchored_minor_q2_lemma_needs_xi_input"
    )
    return {
        "schema": "rh_h810_anchored_minor_q2_lemma.v0",
        "classification": classification,
        "lemma": {
            "name": "anchored fixed-left minor q2 lemma",
            "statement": (
                "Let x_n=R_n/R1, x_1=1, q_n=x_n/x_{n-1}. If 0<q_n<=q_{n+1}<=1 "
                "for n>=2 and q_2<=1/2, then G_c=x_{c-1}(1-x_c) is nonincreasing "
                "for c>=2. Hence all fixed-left anchored minors rows (0,1,2), "
                "cols (1,c,c+1) are nonnegative."
            ),
            "proof_sketch": [
                "For c=2, G_2-G_3=(1-q2)-q2(1-q2*q3) >= 1-2*q2 >= 0.",
                "For c>=3, x_{c-1}<=x_2=q2<=1/2.",
                "G_c-G_{c+1}=x_{c-1}*((1-q_c)-x_{c-1}*q_c*(1-q_c*q_{c+1})).",
                "It is enough that x_{c-1}<=(1-q_c)/(q_c*(1-q_c*q_{c+1})).",
                "Since q_{c+1}>=q_c and q_c<=1, the right side is >=1/(q_c*(1+q_c))>=1/2.",
            ],
        },
        "h809_repair": {
            "invalidated_step": "translated contiguous D3 was incorrectly used to prove G_c monotone",
            "replacement_input": "q monotonicity plus q2<=1/2 proves the anchored fixed-left minors",
            "remaining_inputs_for_consecutive_rows": [
                "PF2 for first_col_0",
                "H810 anchored lemma for first_col_1",
                "translated contiguous D3/local convexity for first_col_ge2",
            ],
        },
        "h796_diagnosis": h796,
        "h801_xi_diagnosis": h801,
        "limitations": [
            "This repairs only the consecutive-row sparse-column part of order-3 PF.",
            "The Xi conclusion is finite to the H801 audited window unless global q-monotonicity is proved.",
            "Translated contiguous D3 and sparse-row minors remain separate bottlenecks.",
        ],
        "decision": (
            "H810 is a useful repair: H796 invalidates PF2+translated-D3, but it fails q2<=1/2. "
            "The Xi finite window satisfies q2<1/2 and q-monotonicity, so the anchored-minor "
            "gap found in H809 is aligned with the H801/H804 route rather than a dead end."
        ),
        "next_target": {
            "name": "H811 repaired consecutive-row PF3 lemma",
            "statement": (
                "Combine H809 identities with H810 anchored lemma and translated contiguous D3 "
                "to formalize a repaired consecutive-row order-3 theorem; then move to sparse rows."
            ),
        },
    }


def write_markdown(report: dict[str, Any], out_path: Path) -> None:
    lines = [
        "# H810 Anchored-Minor q2 Lemma",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Lemma",
        "",
        report["lemma"]["statement"],
        "",
        "## Proof Sketch",
        "",
    ]
    for item in report["lemma"]["proof_sketch"]:
        lines.append(f"- {item}")
    lines.extend(
        [
            "",
            "## H809 Repair",
            "",
            "```json",
            json.dumps(report["h809_repair"], indent=2),
            "```",
            "",
            "## H796 Diagnosis",
            "",
            "```json",
            json.dumps(report["h796_diagnosis"], indent=2),
            "```",
            "",
            "## H801 Xi Diagnosis",
            "",
            "```json",
            json.dumps(report["h801_xi_diagnosis"], indent=2),
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
    parser = argparse.ArgumentParser(description="H810 anchored fixed-left minor q2 lemma audit.")
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
                "h796_q2_le_half": report["h796_diagnosis"]["q2_le_half"],
                "h801_q2_upper_less_than_half": report["h801_xi_diagnosis"]["q2_upper_less_than_half"],
                "h801_q_margins_positive": report["h801_xi_diagnosis"]["all_q_next_minus_q_certified_positive"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
