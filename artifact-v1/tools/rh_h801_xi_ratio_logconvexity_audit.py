import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from rh_h218_jensen_arb_discriminants import gamma_coefficients, load_flint


DEFAULT_OUT = "research/riemann/h801_xi_ratio_logconvexity_audit.json"
DEFAULT_FLINT_PATH = "C:/tmp/h784_pydeps_copy"


def arb_record(value: Any) -> dict[str, Any]:
    return {
        "repr": str(value),
        "lower": str(value.lower()),
        "upper": str(value.upper()),
        "lower_float": float(value.lower()),
        "upper_float": float(value.upper()),
        "mid": str(value.mid()),
        "rad": str(value.rad()),
        "positive_lower_bound": bool(value.lower() > 0),
        "negative_upper_bound": bool(value.upper() < 0),
    }


def acb_radius_record(value: Any) -> dict[str, str]:
    return {
        "repr": str(value),
        "real": arb_record(value.real),
        "imag": arb_record(value.imag),
    }


def midpoint_float(value: Any) -> float:
    return float(value.mid())


def make_ordinary_sequence(arb: Any, gammas: list[Any]) -> list[Any]:
    return [gamma / arb(str(math.factorial(index))) for index, gamma in enumerate(gammas)]


def ratios(values: list[Any]) -> list[Any]:
    return [None] + [values[index] / values[index - 1] for index in range(1, len(values))]


def quotient_ratios(r_values: list[Any]) -> list[Any]:
    return [None, None] + [r_values[index] / r_values[index - 1] for index in range(2, len(r_values))]


def build_config(
    acb: Any,
    arb: Any,
    n_max: int,
    radius: str,
    samples: int,
) -> dict[str, Any]:
    gammas = gamma_coefficients(acb, arb, n_max, radius, samples)
    ordinary = make_ordinary_sequence(arb, gammas)
    r_values = ratios(ordinary)
    q_values = quotient_ratios(r_values)
    q_margins = [None, None] + [
        q_values[index + 1] - q_values[index]
        for index in range(2, len(q_values) - 1)
    ]
    r_decrease_margins = [None, None] + [
        r_values[index - 1] - r_values[index]
        for index in range(2, len(r_values))
    ]
    b = r_values[2]
    h799_tail_margins = [None, None, None] + [
        b - q_values[m] * r_values[m - 1] * (1 + q_values[m])
        for m in range(3, len(q_values) - 1)
    ]
    h799_p_minus_q = [None, None, None] + [
        q_values[m + 1] - q_values[m]
        for m in range(3, len(q_values) - 1)
    ]
    h799_one_minus_p = [None, None, None] + [
        1 - q_values[m + 1]
        for m in range(3, len(q_values) - 1)
    ]
    return {
        "radius": radius,
        "gammas": gammas,
        "ordinary": ordinary,
        "ratios": r_values,
        "q": q_values,
        "q_margins": q_margins,
        "r_decrease_margins": r_decrease_margins,
        "h799_tail_margins": h799_tail_margins,
        "h799_p_minus_q": h799_p_minus_q,
        "h799_one_minus_p": h799_one_minus_p,
    }


def selected_indices(n_max: int) -> list[int]:
    wanted = [0, 1, 2, 3, 4, 5, 8, 12, 16, 22, 24, 32, 40, 48, 56, 64, 80]
    return [index for index in wanted if index <= n_max]


def min_row(records: list[dict[str, Any]], key: str) -> dict[str, Any] | None:
    usable = [row for row in records if row[key] is not None]
    if not usable:
        return None
    return min(usable, key=lambda row: row[key]["lower_float"])


def max_mid_spread(values_by_radius: dict[str, list[Any]], index: int) -> dict[str, Any]:
    mids = {radius: midpoint_float(values[index]) for radius, values in values_by_radius.items()}
    low = min(mids.values())
    high = max(mids.values())
    center = (abs(low) + abs(high)) / 2
    spread = high - low
    relative = None if center == 0 else abs(spread) / center
    return {
        "index": index,
        "midpoints": mids,
        "spread": spread,
        "relative_spread": relative,
    }


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    flint, acb, arb, ctx = load_flint(args.flint_path)
    ctx.dps = args.dps
    radii = [item.strip() for item in args.radii.split(",") if item.strip()]
    if args.reference_radius not in radii:
        raise ValueError("--reference-radius must be included in --radii")

    configs = {
        radius: build_config(acb, arb, args.n_max, radius, args.samples)
        for radius in radii
    }
    reference = configs[args.reference_radius]
    selected = selected_indices(args.n_max)
    selected_gamma_rows = [
        {
            "n": index,
            "gamma": arb_record(reference["gammas"][index]),
            "a_n": arb_record(reference["ordinary"][index]),
        }
        for index in selected
    ]
    selected_ratio_rows = [
        {
            "n": index,
            "R_n": None if index < 1 else arb_record(reference["ratios"][index]),
            "q_n": None if index < 2 else arb_record(reference["q"][index]),
            "R_decrease_margin": None if index < 2 else arb_record(reference["r_decrease_margins"][index]),
            "q_logconvex_margin": None
            if index < 2 or index + 1 > args.n_max
            else arb_record(reference["q_margins"][index]),
        }
        for index in selected
    ]
    q_rows = [
        {
            "n": index,
            "q_n": arb_record(reference["q"][index]),
            "q_next_minus_q": arb_record(reference["q_margins"][index]),
        }
        for index in range(2, args.n_max)
    ]
    h799_rows = [
        {
            "m": m,
            "tail_margin": arb_record(reference["h799_tail_margins"][m]),
            "p_minus_q": arb_record(reference["h799_p_minus_q"][m]),
            "one_minus_p": arb_record(reference["h799_one_minus_p"][m]),
            "condition_certified": bool(
                reference["h799_tail_margins"][m].lower() > 0
                and reference["h799_p_minus_q"][m].lower() >= 0
                and reference["h799_one_minus_p"][m].lower() > 0
            ),
        }
        for m in range(3, args.n_max - 1)
    ]

    q_by_radius = {radius: config["q"] for radius, config in configs.items()}
    r_by_radius = {radius: config["ratios"] for radius, config in configs.items()}
    q_stability_rows = [max_mid_spread(q_by_radius, index) for index in range(2, args.n_max + 1)]
    r_stability_rows = [max_mid_spread(r_by_radius, index) for index in range(1, args.n_max + 1)]
    worst_q_stability = max(q_stability_rows, key=lambda row: row["relative_spread"] or 0.0)
    worst_r_stability = max(r_stability_rows, key=lambda row: row["relative_spread"] or 0.0)

    min_q_margin = min_row(q_rows, "q_next_minus_q")
    min_r_decrease = min_row(selected_ratio_rows, "R_decrease_margin")
    min_tail = min_row(h799_rows, "tail_margin")
    h799_base_cases = [row["m"] for row in h799_rows if not row["condition_certified"]]
    worst_q_relative = worst_q_stability["relative_spread"]
    worst_r_relative = worst_r_stability["relative_spread"]
    checks = {
        "all_gamma_positive_reference": all(item.lower() > 0 for item in reference["gammas"]),
        "all_a_positive_reference": all(item.lower() > 0 for item in reference["ordinary"]),
        "all_R_decrease_margins_positive_reference": all(
            reference["r_decrease_margins"][index].lower() > 0
            for index in range(2, args.n_max + 1)
        ),
        "all_q_logconvex_margins_positive_reference": all(
            reference["q_margins"][index].lower() > 0
            for index in range(2, args.n_max)
        ),
        "all_h799_tail_conditions_certified_reference": all(row["condition_certified"] for row in h799_rows),
        "all_h799_tail_conditions_certified_from_m4_reference": all(row["condition_certified"] for row in h799_rows if row["m"] >= 4),
        "h799_base_cases_not_tail_certified": h799_base_cases,
        "q_cross_radius_relative_spread_small": bool(
            worst_q_relative is not None and worst_q_relative < float(args.stability_threshold)
        ),
        "R_cross_radius_relative_spread_small": bool(
            worst_r_relative is not None and worst_r_relative < float(args.stability_threshold)
        ),
    }
    if not checks["all_q_logconvex_margins_positive_reference"]:
        classification = "xi_ratio_logconvexity_candidate_failure"
        decision = (
            "The reference extraction contains a nonpositive Arb lower bound for q_{n+1}-q_n. "
            "This is a candidate invalidation of H799's ratio-log-convexity route; inspect cross-radius stability before treating it as real."
        )
    elif checks["q_cross_radius_relative_spread_small"]:
        classification = "xi_ratio_logconvexity_survives_extended_stable_window"
        decision = (
            "The H801 window certifies positive q_{n+1}-q_n margins in the reference Arb extraction, "
            "with small cross-radius midpoint spread. This strengthens H799 as a proof target but remains finite Cauchy/DFT evidence."
        )
    else:
        classification = "xi_ratio_logconvexity_survives_reference_window_stability_uncertain"
        decision = (
            "The reference Arb extraction has positive q_{n+1}-q_n margins, but cross-radius spread is not small enough. "
            "Improve coefficient extraction before using the margins as evidence."
        )

    report = {
        "schema": "rh_h801_xi_ratio_logconvexity_audit.v0",
        "classification": classification,
        "parameters": {
            "n_max": args.n_max,
            "samples": args.samples,
            "radii": radii,
            "reference_radius": args.reference_radius,
            "dps": args.dps,
            "flint_path": args.flint_path,
            "stability_threshold": args.stability_threshold,
            "python_flint_version": getattr(flint, "__version__", "unknown"),
        },
        "definitions": {
            "a_n": "a_n=gamma_n/n!",
            "R_n": "R_n=a_n/a_{n-1}",
            "q_n": "q_n=R_n/R_{n-1}",
            "target": "q_{n+1}-q_n >= 0 for ratio-log-convexity of R_n",
            "h799_tail": "B-q_m*R_{m-1}*(1+q_m)>0 plus q_m<=q_{m+1}<1",
        },
        "checks": checks,
        "selected_gamma_rows": selected_gamma_rows,
        "selected_ratio_rows": selected_ratio_rows,
        "q_margin_summary": {
            "min_q_next_minus_q": min_q_margin,
            "rows": q_rows,
        },
        "h799_tail_summary": {
            "min_tail_margin": min_tail,
            "rows": h799_rows,
        },
        "cross_radius_stability": {
            "worst_q": worst_q_stability,
            "worst_R": worst_r_stability,
            "selected_q": [row for row in q_stability_rows if row["index"] in selected],
            "selected_R": [row for row in r_stability_rows if row["index"] in selected],
        },
        "limitations": [
            "The coefficient extraction is finite Cauchy/DFT Arb arithmetic without a fresh rigorous outer-circle alias bound.",
            "Cross-radius agreement is a stability test, not a theorem.",
            "This audits only a finite coefficient window.",
            "Positive H801 margins do not prove H798 globally, PF3, JP2', or RH.",
        ],
        "decision": decision,
        "next_target": {
            "name": "H802 ratio-log-convexity mechanism",
            "statement": (
                "If H801 remains stable, derive q_{n+1}>=q_n from an integral representation "
                "or from explicit coefficient asymptotics with a finite certified cutoff."
            ),
        },
    }
    return report


def write_markdown(report: dict[str, Any], out_path: Path) -> None:
    lines = [
        "# H801 Xi Ratio-Log-Convexity Audit",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Parameters",
        "",
        "| parameter | value |",
        "| --- | ---: |",
    ]
    for key, value in report["parameters"].items():
        lines.append(f"| {key} | `{value}` |")
    lines.extend(
        [
            "",
            "## Checks",
            "",
            "| check | value |",
            "| --- | --- |",
        ]
    )
    for key, value in report["checks"].items():
        lines.append(f"| {key} | `{value}` |")

    min_q = report["q_margin_summary"]["min_q_next_minus_q"]
    min_tail = report["h799_tail_summary"]["min_tail_margin"]
    lines.extend(
        [
            "",
            "## Margin Summary",
            "",
            "Minimum `q_{n+1}-q_n`:",
            "",
            "```json",
            json.dumps(min_q, indent=2),
            "```",
            "",
            "Minimum H799 tail margin:",
            "",
            "```json",
            json.dumps(min_tail, indent=2),
            "```",
            "",
            "Worst cross-radius q stability:",
            "",
            "```json",
            json.dumps(report["cross_radius_stability"]["worst_q"], indent=2),
            "```",
            "",
            "## Selected Ratios",
            "",
            "| n | R_n | q_n | R decrease margin | q next minus q |",
            "| ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in report["selected_ratio_rows"]:
        r = None if row["R_n"] is None else row["R_n"]["repr"]
        q = None if row["q_n"] is None else row["q_n"]["repr"]
        rd = None if row["R_decrease_margin"] is None else row["R_decrease_margin"]["repr"]
        qm = None if row["q_logconvex_margin"] is None else row["q_logconvex_margin"]["repr"]
        lines.append(f"| {row['n']} | `{r}` | `{q}` | `{rd}` | `{qm}` |")

    lines.extend(
        [
            "",
            "## H799 Tail Conditions",
            "",
            "| m | tail margin | p-q | 1-p | certified |",
            "| ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in report["h799_tail_summary"]["rows"]:
        if row["m"] <= 8 or row["m"] in {12, 16, 22, 32, 48, report["parameters"]["n_max"] - 2}:
            lines.append(
                f"| {row['m']} | `{row['tail_margin']['repr']}` | `{row['p_minus_q']['repr']}` | "
                f"`{row['one_minus_p']['repr']}` | `{row['condition_certified']}` |"
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
    parser = argparse.ArgumentParser(description="H801 audit for Xi ratio log-convexity q_{n+1}>=q_n.")
    parser.add_argument("--n-max", type=int, default=48)
    parser.add_argument("--samples", type=int, default=2048)
    parser.add_argument("--radii", default="4,5,6")
    parser.add_argument("--reference-radius", default="5")
    parser.add_argument("--dps", type=int, default=130)
    parser.add_argument("--flint-path", default=DEFAULT_FLINT_PATH)
    parser.add_argument("--stability-threshold", default="1e-40")
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
                "n_max": report["parameters"]["n_max"],
                "all_q_margins_positive": report["checks"]["all_q_logconvex_margins_positive_reference"],
                "all_h799_tail_conditions_from_m4": report["checks"]["all_h799_tail_conditions_certified_from_m4_reference"],
                "worst_q_relative_spread": report["cross_radius_stability"]["worst_q"]["relative_spread"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()






