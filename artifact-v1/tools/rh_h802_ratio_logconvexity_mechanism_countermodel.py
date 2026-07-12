import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any


DEFAULT_XI_SOURCE = "research/riemann/h801_xi_ratio_logconvexity_audit.json"
DEFAULT_OUT = "research/riemann/h802_ratio_logconvexity_mechanism_countermodel.json"


def fraction_record(value: Fraction) -> dict[str, Any]:
    return {
        "fraction": f"{value.numerator}/{value.denominator}",
        "decimal": format(float(value), ".17g"),
        "positive": value > 0,
        "negative": value < 0,
    }


def build_sequence_from_ratios(ratios: list[Fraction]) -> list[Fraction]:
    values = [Fraction(1)]
    for ratio in ratios[1:]:
        values.append(values[-1] * ratio)
    return values


def finite_profile(ratios: list[Fraction]) -> dict[str, Any]:
    values = build_sequence_from_ratios(ratios)
    q_values: list[Fraction | None] = [None, None]
    for index in range(2, len(ratios)):
        q_values.append(ratios[index] / ratios[index - 1])
    q_margins: list[Fraction | None] = [None, None]
    for index in range(2, len(q_values) - 1):
        q_margins.append(q_values[index + 1] - q_values[index])
    ratio_decrease = [
        ratios[index - 1] - ratios[index]
        for index in range(2, len(ratios))
    ]
    h799_tail = []
    b = ratios[2]
    for m in range(3, len(ratios) - 1):
        q = q_values[m]
        p = q_values[m + 1]
        assert q is not None and p is not None
        h799_tail.append(
            {
                "m": m,
                "tail_margin": b - q * ratios[m - 1] * (1 + q),
                "p_minus_q": p - q,
                "one_minus_p": 1 - p,
            }
        )
    return {
        "values": values,
        "ratios": ratios,
        "q_values": q_values,
        "q_margins": q_margins,
        "ratio_decrease": ratio_decrease,
        "h799_tail": h799_tail,
    }


def exact_counterexample() -> dict[str, Any]:
    # R_0 is unused. Ratios are positive, nonincreasing, and can be continued
    # by R_n=R_5/2^(n-5), forcing infinite radius, while q_4<q_3.
    ratios = [
        Fraction(0),
        Fraction(1, 2),
        Fraction(1, 4),
        Fraction(1, 5),
        Fraction(1, 20),
        Fraction(1, 50),
    ]
    profile = finite_profile(ratios)
    q3 = profile["q_values"][3]
    q4 = profile["q_values"][4]
    assert q3 is not None and q4 is not None
    return {
        "definition": {
            "a0": "1",
            "ratios": {
                f"R_{index}": fraction_record(value)
                for index, value in enumerate(ratios)
                if index >= 1
            },
            "tail_continuation": "For n>=5, set R_{n+1}=R_n/2. Then R_n decreases to 0 geometrically and the ordinary generating function has infinite radius.",
        },
        "a_values": {
            f"a_{index}": fraction_record(value)
            for index, value in enumerate(profile["values"])
        },
        "checks": {
            "positive_values": all(value > 0 for value in profile["values"]),
            "positive_ratios": all(value > 0 for value in ratios[1:]),
            "ratios_nonincreasing": all(value >= 0 for value in profile["ratio_decrease"]),
            "finite_log_concavity": all(value >= 0 for value in profile["ratio_decrease"]),
            "entire_tail_possible": True,
            "q_monotonicity_fails": q4 - q3 < 0,
        },
        "q_values": {
            f"q_{index}": None if value is None else fraction_record(value)
            for index, value in enumerate(profile["q_values"])
        },
        "q_failure": {
            "q_3": fraction_record(q3),
            "q_4": fraction_record(q4),
            "q_4_minus_q_3": fraction_record(q4 - q3),
        },
        "h799_tail_rows": [
            {
                "m": row["m"],
                "tail_margin": fraction_record(row["tail_margin"]),
                "p_minus_q": fraction_record(row["p_minus_q"]),
                "one_minus_p": fraction_record(row["one_minus_p"]),
                "condition_pass": row["tail_margin"] > 0 and row["p_minus_q"] >= 0 and row["one_minus_p"] > 0,
            }
            for row in profile["h799_tail"]
        ],
    }


def xi_profile(source_path: Path) -> dict[str, Any]:
    source = json.loads(source_path.read_text(encoding="utf-8"))
    q_min = source["q_margin_summary"]["min_q_next_minus_q"]
    tail_rows = source["h799_tail_summary"]["rows"]
    certified_from_m4 = [
        row for row in tail_rows if row["m"] >= 4 and row["condition_certified"]
    ]
    return {
        "source": str(source_path),
        "classification": source["classification"],
        "parameters": source["parameters"],
        "checks": source["checks"],
        "min_q_margin": q_min,
        "certified_h799_tail_from_m4_count": len(certified_from_m4),
        "total_h799_tail_from_m4_count": len([row for row in tail_rows if row["m"] >= 4]),
    }


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    counterexample = exact_counterexample()
    xi = xi_profile(Path(args.xi_source))
    return {
        "schema": "rh_h802_ratio_logconvexity_mechanism_countermodel.v0",
        "classification": "generic_ratio_logconvexity_mechanism_invalidated_xi_special_needed",
        "question": (
            "Does q_{n+1}>=q_n follow from positivity, adjacent log-concavity, "
            "decreasing ratios with R_n->0, or generic H799-style tail smallness?"
        ),
        "counterexample": counterexample,
        "xi_comparison": xi,
        "decision": (
            "The generic mechanism is false: positive log-concave sequences with decreasing ratios "
            "and an entire tail can have q_{n+1}<q_n. Therefore H801's Xi behavior is not explained "
            "by the broad properties isolated in H800. H802 must pivot to a Xi-specific mechanism, "
            "such as an integral representation with a special kernel or explicit coefficient asymptotics."
        ),
        "not_proved": [
            "Ratio log-convexity for Xi coefficients.",
            "A Xi-specific analytic mechanism.",
            "The H798 family for all m.",
            "PF3, JP2', PF-infinity, or RH.",
        ],
        "next_target": {
            "name": "H803 Xi-specific ratio-log-convexity mechanism",
            "statement": (
                "Use the Riemann Xi integral representation or coefficient asymptotics to seek a "
                "special proof of q_{n+1}>=q_n, since generic log-concavity/entire-tail assumptions are insufficient."
            ),
        },
    }


def write_markdown(report: dict[str, Any], out_path: Path) -> None:
    counterexample = report["counterexample"]
    lines = [
        "# H802 Ratio-Log-Convexity Mechanism Countermodel",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Question",
        "",
        report["question"],
        "",
        "## Exact Counterexample",
        "",
        counterexample["definition"]["tail_continuation"],
        "",
        "| item | value |",
        "| --- | ---: |",
    ]
    for key, value in counterexample["definition"]["ratios"].items():
        lines.append(f"| `{key}` | `{value['fraction']}` |")
    lines.extend(
        [
            "",
            "Checks:",
            "",
            "| check | value |",
            "| --- | --- |",
        ]
    )
    for key, value in counterexample["checks"].items():
        lines.append(f"| {key} | `{value}` |")
    lines.extend(
        [
            "",
            "Ratio-log-convexity failure:",
            "",
            "```json",
            json.dumps(counterexample["q_failure"], indent=2),
            "```",
            "",
            "H799 rows in the counterexample:",
            "",
            "| m | tail margin | p-q | 1-p | pass |",
            "| ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in counterexample["h799_tail_rows"]:
        lines.append(
            f"| {row['m']} | `{row['tail_margin']['fraction']}` | `{row['p_minus_q']['fraction']}` | "
            f"`{row['one_minus_p']['fraction']}` | `{row['condition_pass']}` |"
        )
    xi = report["xi_comparison"]
    lines.extend(
        [
            "",
            "## Xi Comparison",
            "",
            f"H801 classification: `{xi['classification']}`",
            "",
            f"H801 parameters: `{xi['parameters']}`",
            "",
            f"Minimum H801 q margin lower bound: `{xi['min_q_margin']['q_next_minus_q']['lower']}`",
            "",
            (
                "H799 tail certified from m=4 count: "
                f"`{xi['certified_h799_tail_from_m4_count']}/{xi['total_h799_tail_from_m4_count']}`"
            ),
            "",
            "## Decision",
            "",
            report["decision"],
            "",
            "## Not Proved",
            "",
        ]
    )
    for item in report["not_proved"]:
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
    parser = argparse.ArgumentParser(description="H802 generic mechanism countermodel for q ratio log-convexity.")
    parser.add_argument("--xi-source", default=DEFAULT_XI_SOURCE)
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
                "q_failure": report["counterexample"]["q_failure"],
                "xi_classification": report["xi_comparison"]["classification"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
