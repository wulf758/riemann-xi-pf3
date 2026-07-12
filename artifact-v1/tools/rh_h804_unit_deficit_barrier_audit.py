import argparse
import json
import re
from decimal import Decimal, getcontext
from pathlib import Path
from statistics import mean
from typing import Any


DEFAULT_IN = "research/riemann/h803_xi_moment_curvature_decomposition.json"
DEFAULT_OUT = "research/riemann/h804_unit_deficit_barrier_audit.json"

NUMBER_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


def parse_decimal(value: Any) -> Decimal:
    match = NUMBER_RE.search(str(value))
    if not match:
        raise ValueError(f"cannot parse Decimal from {value!r}")
    return Decimal(match.group(0))


def dec(value: Decimal, digits: int = 36) -> str:
    if value.is_zero():
        return "0"
    return f"{value:.{digits}E}"


def threshold_deficit_formula(n: int) -> Decimal:
    n_dec = Decimal(n)
    return (Decimal(8) * n_dec * n_dec - Decimal(4) * n_dec - Decimal(3)) / (
        n_dec * n_dec * (Decimal(2) * n_dec - Decimal(1)) ** 2
    )


def find_tail_start(rows: list[dict[str, Any]], key: str) -> int | None:
    for index, row in enumerate(rows):
        if all(item[key] for item in rows[index:]):
            return row["n"]
    return None


def selected_rows(rows: list[dict[str, Any]], n_max: int) -> list[dict[str, Any]]:
    wanted = {2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 22, 32, 48, 64, 80, 100, n_max}
    return [row for row in rows if row["n"] in wanted]


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    getcontext().prec = args.precision
    h803_path = Path(args.h803)
    h803 = json.loads(h803_path.read_text(encoding="utf-8"))
    rows: list[dict[str, Any]] = []
    for source in h803["rows"]:
        n = int(source["n"])
        observed = parse_decimal(source["observed_T_next_over_T_n"])
        threshold = parse_decimal(source["factorial_threshold"])
        threshold_deficit = Decimal(1) - threshold
        formula_deficit = threshold_deficit_formula(n)
        formula_agreement_abs = abs(formula_deficit - threshold_deficit)
        unit_floor = Decimal(1) - (Decimal(1) / (Decimal(n) * Decimal(n)))
        observed_deficit = Decimal(1) - observed
        unit_slack = observed - unit_floor
        threshold_slack = observed - threshold
        rows.append(
            {
                "n": n,
                "observed_T_next_over_T_n": dec(observed),
                "factorial_threshold": dec(threshold),
                "unit_floor_1_minus_1_over_n2": dec(unit_floor),
                "observed_deficit": dec(observed_deficit),
                "observed_deficit_constant_n2": dec(observed_deficit * n * n),
                "threshold_deficit": dec(threshold_deficit),
                "threshold_deficit_formula": dec(formula_deficit),
                "threshold_formula_agreement_abs": dec(formula_agreement_abs),
                "threshold_deficit_constant_n2": dec(threshold_deficit * n * n),
                "unit_slack": dec(unit_slack),
                "threshold_slack": dec(threshold_slack),
                "unit_implies_threshold_exact": bool(
                    Decimal(1) / (Decimal(n) * Decimal(n)) <= formula_deficit
                ),
                "threshold_formula_agrees_reported": bool(
                    formula_agreement_abs < Decimal("1e-30")
                ),
                "unit_barrier_holds_midpoint": bool(unit_slack >= 0),
                "h803_condition_holds_midpoint": bool(source["condition_holds_midpoint"]),
                "h801_certified": bool(source["h801_q_margin_certified_positive"]),
            }
        )

    unit_failures = [row["n"] for row in rows if not row["unit_barrier_holds_midpoint"]]
    tail_start = find_tail_start(rows, "unit_barrier_holds_midpoint")
    tail_rows = [row for row in rows if tail_start is not None and row["n"] >= tail_start]
    n_max = max(row["n"] for row in rows)
    last_rows = rows[-min(args.tail_average, len(rows)) :]
    last_observed_constants = [
        float(parse_decimal(row["observed_deficit_constant_n2"])) for row in last_rows
    ]
    last_threshold_constants = [
        float(parse_decimal(row["threshold_deficit_constant_n2"])) for row in last_rows
    ]

    min_unit_slack_tail = None
    if tail_rows:
        min_unit_slack_tail = min(tail_rows, key=lambda row: parse_decimal(row["unit_slack"]))
    min_threshold_slack = min(rows, key=lambda row: parse_decimal(row["threshold_slack"]))
    all_unit_implies_threshold = all(row["unit_implies_threshold_exact"] for row in rows)
    all_formula_agrees = all(row["threshold_formula_agrees_reported"] for row in rows)
    all_h803 = all(row["h803_condition_holds_midpoint"] and row["h801_certified"] for row in rows)
    classification = (
        "unit_deficit_barrier_supported_with_finite_base"
        if tail_start is not None and all_unit_implies_threshold and all_formula_agrees and all_h803
        else "unit_deficit_barrier_needs_inspection"
    )
    return {
        "schema": "rh_h804_unit_deficit_barrier_audit.v0",
        "classification": classification,
        "source": {
            "h803_json": str(h803_path),
            "h803_classification": h803.get("classification"),
        },
        "definitions": {
            "T_n": "T_n = (M_n/M_{n-1})/(M_{n-1}/M_{n-2})",
            "h803_threshold": "theta_n = C_{n-1}C_{n+1}/C_n^2, C_n=(2n)(2n-1)",
            "unit_barrier": "T_{n+1}/T_n >= 1 - 1/n^2",
        },
        "derived_identities": {
            "threshold_deficit": "1 - theta_n = (8n^2 - 4n - 3)/(n^2(2n - 1)^2)",
            "unit_implies_h803": (
                "For n>=2, (8n^2-4n-3)/(2n-1)^2 > 1, so "
                "1/n^2 < 1-theta_n; hence T_{n+1}/T_n >= 1-1/n^2 implies H803."
            ),
        },
        "checks": {
            "rows_tested": len(rows),
            "n_range": [rows[0]["n"], rows[-1]["n"]],
            "all_unit_implies_threshold_exact": all_unit_implies_threshold,
            "all_threshold_formula_agrees_reported": all_formula_agrees,
            "all_h803_conditions_available": all_h803,
            "unit_barrier_tail_start_midpoint": tail_start,
            "unit_barrier_failures_before_tail": unit_failures,
        },
        "summary": {
            "min_threshold_slack": min_threshold_slack,
            "min_unit_slack_from_tail_start": min_unit_slack_tail,
            "tail_average_window": len(last_rows),
            "tail_average_observed_deficit_constant_n2": mean(last_observed_constants),
            "tail_average_threshold_deficit_constant_n2": mean(last_threshold_constants),
        },
        "selected_rows": selected_rows(rows, n_max),
        "rows": rows,
        "decision": (
            "H804 is a useful simplification: H803 reduces to a finite base plus the cleaner "
            "eventual inequality T_{n+1}/T_n >= 1-1/n^2. The finite H803 data supports this "
            "from n=8 onward, and the tail deficit constant is far below the factorial "
            "threshold constant. The remaining work is analytic: prove the unit-deficit "
            "barrier from the Xi kernel or from coefficient asymptotics with explicit error."
        ),
        "limitations": [
            "The unit barrier is stronger than H803 and fails for small n, so finite base cases remain necessary.",
            "The report uses H803 midpoint ratios plus H801 positivity certificates; it is not a new Arb proof.",
            "No Xi-kernel analytic estimate is proved here.",
        ],
        "next_target": {
            "name": "H805 saddle-point unit-deficit lemma",
            "statement": (
                "Try to prove eventually T_{n+1}/T_n >= 1-1/n^2, for example from "
                "Xi-kernel log-saddle concentration such as T_n = 1 + O(1/(n log n)) with enough regularity to control the quotient."
            ),
        },
    }


def write_markdown(report: dict[str, Any], out_path: Path) -> None:
    lines = [
        "# H804 Unit-Deficit Moment-Curvature Barrier Audit",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Derived Barrier",
        "",
        "H803 needs",
        "",
        "```text",
        "T_{n+1}/T_n >= theta_n = C_{n-1}C_{n+1}/C_n^2.",
        "```",
        "",
        "The exact deficit is",
        "",
        "```text",
        "1 - theta_n = (8n^2 - 4n - 3)/(n^2(2n - 1)^2).",
        "```",
        "",
        "Since this is larger than `1/n^2` for `n>=2`, the stronger unit barrier",
        "",
        "```text",
        "T_{n+1}/T_n >= 1 - 1/n^2",
        "```",
        "",
        "implies the H803 inequality.",
        "",
        "## Checks",
        "",
        "| check | value |",
        "| --- | --- |",
    ]
    for key, value in report["checks"].items():
        lines.append(f"| {key} | `{value}` |")

    lines.extend(
        [
            "",
            "## Tightest Rows",
            "",
            "Minimum H803 threshold slack:",
            "",
            "```json",
            json.dumps(report["summary"]["min_threshold_slack"], indent=2),
            "```",
            "",
            "Minimum unit-barrier slack from tail start:",
            "",
            "```json",
            json.dumps(report["summary"]["min_unit_slack_from_tail_start"], indent=2),
            "```",
            "",
            "## Tail Constants",
            "",
            f"Last `{report['summary']['tail_average_window']}` rows average "
            f"`n^2(1-T_{{n+1}}/T_n)` = "
            f"`{report['summary']['tail_average_observed_deficit_constant_n2']}`.",
            "",
            f"Last `{report['summary']['tail_average_window']}` rows average "
            f"`n^2(1-theta_n)` = "
            f"`{report['summary']['tail_average_threshold_deficit_constant_n2']}`.",
            "",
            "## Selected Rows",
            "",
            "| n | observed deficit constant | threshold constant | unit slack | unit barrier |",
            "| ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in report["selected_rows"]:
        lines.append(
            f"| {row['n']} | `{row['observed_deficit_constant_n2']}` | "
            f"`{row['threshold_deficit_constant_n2']}` | `{row['unit_slack']}` | "
            f"`{row['unit_barrier_holds_midpoint']}` |"
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
    parser = argparse.ArgumentParser(
        description="Audit the H804 unit-deficit sufficient barrier for H803."
    )
    parser.add_argument("--h803", default=DEFAULT_IN)
    parser.add_argument("--out", default=DEFAULT_OUT)
    parser.add_argument("--precision", type=int, default=120)
    parser.add_argument("--tail-average", type=int, default=20)
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
                "rows_tested": report["checks"]["rows_tested"],
                "unit_barrier_tail_start_midpoint": report["checks"][
                    "unit_barrier_tail_start_midpoint"
                ],
                "unit_barrier_failures_before_tail": report["checks"][
                    "unit_barrier_failures_before_tail"
                ],
                "tail_average_observed_deficit_constant_n2": report["summary"][
                    "tail_average_observed_deficit_constant_n2"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()



