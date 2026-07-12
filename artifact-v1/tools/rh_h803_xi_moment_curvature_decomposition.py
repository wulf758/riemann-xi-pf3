import argparse
import json
import re
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Any


DEFAULT_IN = "research/riemann/h801_xi_ratio_logconvexity_audit.json"
DEFAULT_OUT = "research/riemann/h803_xi_moment_curvature_decomposition.json"

NUMBER_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


def decimal_from_record(record: dict[str, Any], key: str = "mid") -> Decimal:
    text = str(record.get(key) or record.get("repr"))
    match = NUMBER_RE.search(text)
    if not match:
        raise ValueError(f"cannot parse Decimal from {text!r}")
    return Decimal(match.group(0))


def c_factor(n: int) -> Decimal:
    return Decimal(2 * n * (2 * n - 1))


def dec(value: Decimal, digits: int = 36) -> str:
    if value.is_zero():
        return "0"
    return f"{value:.{digits}E}"


def selected_rows(rows: list[dict[str, Any]], n_max: int) -> list[dict[str, Any]]:
    wanted = {2, 3, 4, 5, 6, 7, 8, 12, 16, 22, 32, 48, 64, 80, 100, n_max}
    return [row for row in rows if row["n"] in wanted]


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    getcontext().prec = args.precision
    h801_path = Path(args.h801)
    h801 = json.loads(h801_path.read_text(encoding="utf-8"))
    q_rows_raw = h801["q_margin_summary"]["rows"]
    q_by_n = {row["n"]: decimal_from_record(row["q_n"]) for row in q_rows_raw}
    q_margin_certified = {
        row["n"]: bool(row["q_next_minus_q"].get("positive_lower_bound"))
        for row in q_rows_raw
    }

    rows: list[dict[str, Any]] = []
    for n in sorted(q_by_n):
        if n + 1 not in q_by_n:
            continue
        c_prev = c_factor(n - 1)
        c_n = c_factor(n)
        c_next = c_factor(n + 1)
        q_n = q_by_n[n]
        q_next = q_by_n[n + 1]
        t_n = q_n * c_n / c_prev
        t_next = q_next * c_next / c_n
        observed = t_next / t_n
        threshold = c_prev * c_next / (c_n * c_n)
        additive_slack = observed - threshold
        normalized_slack = observed / threshold - Decimal(1)
        q_ratio_minus_one = q_next / q_n - Decimal(1)
        rows.append(
            {
                "n": n,
                "q_n": dec(q_n),
                "q_next": dec(q_next),
                "T_n": dec(t_n),
                "T_next": dec(t_next),
                "observed_T_next_over_T_n": dec(observed),
                "factorial_threshold": dec(threshold),
                "additive_slack": dec(additive_slack),
                "normalized_slack": dec(normalized_slack),
                "q_next_over_q_n_minus_1": dec(q_ratio_minus_one),
                "h801_q_margin_certified_positive": q_margin_certified[n],
                "condition_holds_midpoint": bool(additive_slack > 0),
            }
        )

    min_additive = min(rows, key=lambda row: Decimal(row["additive_slack"]))
    min_normalized = min(rows, key=lambda row: Decimal(row["normalized_slack"]))
    all_midpoint = all(row["condition_holds_midpoint"] for row in rows)
    all_certified = all(row["h801_q_margin_certified_positive"] for row in rows)
    n_max = max(row["n"] for row in rows)
    classification = (
        "xi_moment_curvature_reformulation_finite_positive"
        if all_midpoint and all_certified
        else "xi_moment_curvature_reformulation_needs_inspection"
    )

    return {
        "schema": "rh_h803_xi_moment_curvature_decomposition.v0",
        "classification": classification,
        "source": {
            "h801_json": str(h801_path),
            "h801_classification": h801.get("classification"),
            "h801_parameters": h801.get("parameters"),
        },
        "definitions": {
            "a_n": "a_n = gamma_n/n!",
            "M_n": "M_n = (2n)! a_n, the even cosine-moment normalization used in H789",
            "S_n": "S_n = M_n/M_{n-1}",
            "T_n": "T_n = S_n/S_{n-1} = M_n M_{n-2}/M_{n-1}^2",
            "C_n": "C_n = (2n)(2n-1)",
            "R_n": "R_n = a_n/a_{n-1} = S_n/C_n",
            "q_n": "q_n = R_n/R_{n-1} = T_n C_{n-1}/C_n",
        },
        "corrected_equivalence": {
            "target": "q_{n+1} >= q_n",
            "equivalent_moment_curvature_condition": (
                "T_{n+1}/T_n >= C_{n-1} C_{n+1}/C_n^2"
            ),
            "factorial_threshold": (
                "C_{n-1} C_{n+1}/C_n^2 = "
                "((2n-2)(2n-3)(2n+2)(2n+1))/((2n)(2n-1))^2"
            ),
            "invalidated_initial_orientation": (
                "The reciprocal C_n^2/(C_{n-1} C_{n+1}) was the wrong threshold."
            ),
        },
        "checks": {
            "rows_tested": len(rows),
            "n_range": [rows[0]["n"], rows[-1]["n"]],
            "all_h801_q_margins_certified_positive": all_certified,
            "all_corrected_moment_curvature_conditions_hold_midpoint": all_midpoint,
        },
        "summary": {
            "min_additive_slack": min_additive,
            "min_normalized_slack": min_normalized,
        },
        "selected_rows": selected_rows(rows, n_max),
        "rows": rows,
        "limitations": [
            "This is an exact algebraic reformulation plus a finite audit of H801 midpoint data.",
            "Rigorous positivity is inherited only from H801's Arb-certified q_{n+1}-q_n rows.",
            "No analytic proof of the Xi moment-curvature condition is supplied here.",
            "The generic moment/log-concavity route remains invalidated by H802.",
        ],
        "decision": (
            "H803 isolates the surviving Xi-specific sublemma: the even moments M_n must have "
            "enough second-ratio persistence to beat the exact factorial threshold "
            "C_{n-1}C_{n+1}/C_n^2. This is cleaner than treating q-monotonicity as a black box, "
            "but it is still a target lemma, not a proof of JP2' or RH."
        ),
        "next_target": {
            "name": "H804 Xi kernel moment-curvature proof or counterexample",
            "statement": (
                "Use the explicit Riemann Xi kernel or saddle-point asymptotics to prove, "
                "disprove, or finitely reduce T_{n+1}/T_n >= C_{n-1}C_{n+1}/C_n^2."
            ),
        },
    }


def write_markdown(report: dict[str, Any], out_path: Path) -> None:
    lines = [
        "# H803 Xi Moment-Curvature Decomposition",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Corrected Equivalence",
        "",
        "With",
        "",
        "```text",
        "M_n = (2n)! a_n",
        "S_n = M_n/M_{n-1}",
        "T_n = S_n/S_{n-1}",
        "C_n = (2n)(2n-1)",
        "R_n = S_n/C_n",
        "q_n = T_n C_{n-1}/C_n",
        "```",
        "",
        "the H801 bottleneck `q_{n+1} >= q_n` is exactly",
        "",
        "```text",
        "T_{n+1}/T_n >= C_{n-1} C_{n+1}/C_n^2.",
        "```",
        "",
        "The reciprocal threshold from the first H803 hypothesis card is invalidated.",
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
            "Minimum additive slack:",
            "",
            "```json",
            json.dumps(report["summary"]["min_additive_slack"], indent=2),
            "```",
            "",
            "Minimum normalized slack:",
            "",
            "```json",
            json.dumps(report["summary"]["min_normalized_slack"], indent=2),
            "```",
            "",
            "## Selected Rows",
            "",
            "| n | T_next/T_n | threshold | additive slack | normalized slack | H801 certified |",
            "| ---: | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in report["selected_rows"]:
        lines.append(
            f"| {row['n']} | `{row['observed_T_next_over_T_n']}` | "
            f"`{row['factorial_threshold']}` | `{row['additive_slack']}` | "
            f"`{row['normalized_slack']}` | `{row['h801_q_margin_certified_positive']}` |"
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
        description="H803 exact moment-curvature decomposition of Xi q-log-convexity."
    )
    parser.add_argument("--h801", default=DEFAULT_IN)
    parser.add_argument("--out", default=DEFAULT_OUT)
    parser.add_argument("--precision", type=int, default=120)
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
                "n_range": report["checks"]["n_range"],
                "min_additive_slack_n": report["summary"]["min_additive_slack"]["n"],
                "min_normalized_slack_n": report["summary"]["min_normalized_slack"]["n"],
                "all_h801_certified": report["checks"]["all_h801_q_margins_certified_positive"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
