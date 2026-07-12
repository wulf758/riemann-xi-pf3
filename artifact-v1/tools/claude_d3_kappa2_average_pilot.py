"""Rigorous one-row translated-D3 pilot through kappa2 averages.

For M(theta)=int_R exp(theta*s)V(s)ds and Lambda=log(M),

    A_n = (Lambda(2n)-2 Lambda(2n-2)+Lambda(2n-4))/4
        = int_[0,1]^2 kappa2(2(n-2)+2(x+y)) dx dy,

and

    q_n = exp(4 A_n) * (2n-2)(2n-3)/((2n)(2n-1)).

This pilot rigorously integrates the five moments needed for n=127,128,129,
certifies D_127 through both the kappa2-average and raw moment-ratio routes,
and checks q_127,q_128 against the independent H1321 DFT/alias seam.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any


THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

import claude_kappa3_window_certificate as window  # noqa: E402


DEFAULT_OUT = "research/riemann/claude_d3_kappa2_average_pilot.json"
DEFAULT_FINITE_SOURCE = (
    "research/riemann/claude_d3_alias_repaired_finite_certificate.json"
)


def arb_record(value: Any) -> dict[str, Any]:
    return {
        "repr": str(value),
        "lower": str(value.lower()),
        "upper": str(value.upper()),
        "mid": str(value.mid()),
        "rad": str(value.rad()),
        "positive_lower_bound": bool(value.lower() > 0),
        "negative_upper_bound": bool(value.upper() < 0),
        "contains_zero": bool(value.lower() <= 0 <= value.upper()),
    }


def c(index: int) -> int:
    return (2 * index) * (2 * index - 1)


def positive_tail_interval(bound: Any) -> Any:
    return window.arb(0).union(bound.upper())


def rigorous_moment(theta: int) -> tuple[Any, dict[str, Any]]:
    center, sigma, s_lo, s_hi = window.window_params(theta)
    right = window.right_tail_bound(theta, center, s_hi, 0)
    left = window.left_tail_bound(theta, center, s_lo, 0)
    if right is None:
        raise RuntimeError(f"right-tail condition failed at theta={theta}")
    truncated = window.rig_integral(
        theta,
        center,
        0,
        s_lo,
        s_hi,
        sig=sigma,
    )
    tail_bound = right + left
    moment = truncated + positive_tail_interval(tail_bound)
    if not moment.lower() > 0:
        raise RuntimeError(f"moment lower bound is not positive at theta={theta}")
    metadata = {
        "theta": theta,
        "center": center,
        "sigma": sigma,
        "s_lo": s_lo,
        "s_hi": s_hi,
        "truncated_integral": arb_record(truncated),
        "right_tail_upper": str(right.upper()),
        "left_tail_upper": str(left.upper()),
        "total_tail_upper": str(tail_bound.upper()),
        "relative_tail_upper": str(tail_bound.upper() / moment.lower()),
        "moment": arb_record(moment),
    }
    return moment, metadata


def q_from_kappa2_average(
    moments: dict[int, Any],
    n: int,
) -> tuple[Any, Any, Any]:
    log_high = moments[2 * n].log()
    log_mid = moments[2 * n - 2].log()
    log_low = moments[2 * n - 4].log()
    average = (log_high - 2 * log_mid + log_low) / 4
    factorial_ratio = window.arb(c(n - 1)) / window.arb(c(n))
    q_average = (4 * average).exp() * factorial_ratio
    q_direct = (
        moments[2 * n]
        * moments[2 * n - 4]
        / moments[2 * n - 2] ** 2
        * factorial_ratio
    )
    return average, q_average, q_direct


def translated_d3(q_values: dict[int, Any], n: int) -> Any:
    one = q_values[n + 1] * 0 + 1
    return (
        (one - q_values[n + 1]) ** 2
        - q_values[n + 1] ** 2
        * (one - q_values[n])
        * (one - q_values[n + 2])
    )


def parse_factor(text: str) -> tuple[int, int]:
    pieces = text.split(":")
    if len(pieces) != 2:
        raise ValueError(f"invalid rational factor {text!r}")
    numerator, denominator = int(pieces[0]), int(pieces[1])
    if numerator <= 0 or denominator <= 0:
        raise ValueError("mutation factors must be positive")
    return numerator, denominator


def mutated_d3_rows(
    moments: dict[int, Any],
    target_n: int,
    factor_texts: list[str],
) -> list[dict[str, Any]]:
    rows = []
    mutated_theta = 2 * (target_n + 1)
    for factor_text in factor_texts:
        numerator, denominator = parse_factor(factor_text)
        mutated = dict(moments)
        mutated[mutated_theta] *= window.arb(numerator) / window.arb(denominator)
        q_average: dict[int, Any] = {}
        q_direct: dict[int, Any] = {}
        for n in range(target_n, target_n + 3):
            _, q_average[n], q_direct[n] = q_from_kappa2_average(mutated, n)
        d3_average = translated_d3(q_average, target_n)
        d3_direct = translated_d3(q_direct, target_n)
        rows.append(
            {
                "target_n": target_n,
                "mutated_theta": mutated_theta,
                "factor": factor_text,
                "d3_from_kappa2_average": arb_record(d3_average),
                "d3_from_moment_ratio": arb_record(d3_direct),
                "breaks_positive_certificate": bool(
                    d3_average.upper() < 0 and d3_direct.upper() < 0
                ),
            }
        )
    return rows


def finite_overlap_balls(
    finite_source: Path,
    arb: Any,
) -> dict[int, Any]:
    report = json.loads(finite_source.read_text(encoding="utf-8"))
    rows = report["finite_d3"]["rows"]
    row126 = next(row for row in rows if int(row["n"]) == 126)
    return {
        127: arb(row126["q_center"]["repr"]),
        128: arb(row126["q_right"]["repr"]),
    }


def write_markdown(report: dict[str, Any], out_path: Path) -> None:
    result = report.get("result", {})
    lines = [
        "# Rigorous Kappa2-Average D3 Pilot",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Result",
        "",
        report["decision"],
        "",
        "The exact bridge used is",
        "",
        "```text",
        "A_n = [Lambda(2n)-2Lambda(2n-2)+Lambda(2n-4)]/4",
        "    = int_[0,1]^2 kappa2(2(n-2)+2(x+y)) dx dy,",
        "q_n = exp(4 A_n) * (2n-2)(2n-3)/((2n)(2n-1)).",
        "```",
        "",
        f"- target row: `D_{result.get('target_n')}`",
        f"- kappa2-average lower bound: `{result.get('d3_from_kappa2_average', {}).get('lower')}`",
        f"- direct moment-ratio lower bound: `{result.get('d3_from_moment_ratio', {}).get('lower')}`",
        f"- overlap with H1321 for q127/q128: `{result.get('finite_overlap_pass')}`",
        "",
        "## Checks",
        "",
    ]
    for name, value in report.get("checks", {}).items():
        lines.append(f"- `{name}`: `{value}`")
    lines.extend(
        [
            "",
            "## Scope",
            "",
            "This is one certified row beyond the H1321 finite seam. It does not",
            "supply a uniform tail lemma for all `n>=127`, PF-infinity, Jensen",
            "hyperbolicity in all degrees, or RH.",
            "",
        ]
    )
    out_path.write_text("\n".join(lines), encoding="utf-8")


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    window.ctx.prec = args.prec
    target_n = args.target_n
    theta_values = list(range(2 * target_n - 4, 2 * target_n + 5, 2))
    moments: dict[int, Any] = {}
    moment_rows = []
    for theta in theta_values:
        moments[theta], metadata = rigorous_moment(theta)
        moment_rows.append(metadata)

    averages: dict[int, Any] = {}
    q_average: dict[int, Any] = {}
    q_direct: dict[int, Any] = {}
    q_rows = []
    for n in range(target_n, target_n + 3):
        averages[n], q_average[n], q_direct[n] = q_from_kappa2_average(moments, n)
        difference = q_average[n] - q_direct[n]
        q_rows.append(
            {
                "n": n,
                "kappa2_square_average": arb_record(averages[n]),
                "q_from_kappa2_average": arb_record(q_average[n]),
                "q_from_moment_ratio": arb_record(q_direct[n]),
                "cross_evaluation_difference": arb_record(difference),
            }
        )

    d3_average = translated_d3(q_average, target_n)
    d3_direct = translated_d3(q_direct, target_n)
    d3_difference = d3_average - d3_direct

    finite_source = Path(args.finite_source)
    finite_q = finite_overlap_balls(finite_source, window.arb)
    overlap_rows = []
    for n in (target_n, target_n + 1):
        for route_name, value in (
            ("kappa2_average", q_average[n]),
            ("moment_ratio", q_direct[n]),
        ):
            difference = value - finite_q[n]
            overlap_rows.append(
                {
                    "n": n,
                    "route": route_name,
                    "pilot_q": arb_record(value),
                    "finite_q": arb_record(finite_q[n]),
                    "difference": arb_record(difference),
                    "overlaps": bool(difference.lower() <= 0 <= difference.upper()),
                }
            )

    factor_texts = [
        text.strip() for text in args.canary_factors.split(",") if text.strip()
    ]
    canary_rows = mutated_d3_rows(moments, target_n, factor_texts)
    sign_canary = -d3_average
    checks = {
        "five_even_theta_moments_computed": theta_values
        == list(range(2 * target_n - 4, 2 * target_n + 5, 2)),
        "all_moment_lower_bounds_positive": all(
            moments[theta].lower() > 0 for theta in theta_values
        ),
        "all_q_cross_evaluations_overlap": all(
            row["cross_evaluation_difference"]["contains_zero"] for row in q_rows
        ),
        "d3_kappa2_average_lower_bound_positive": bool(d3_average.lower() > 0),
        "d3_moment_ratio_lower_bound_positive": bool(d3_direct.lower() > 0),
        "d3_cross_evaluations_overlap": bool(
            d3_difference.lower() <= 0 <= d3_difference.upper()
        ),
        "all_h1321_overlap_checks_pass": all(row["overlaps"] for row in overlap_rows),
        "moment_mutation_canary_breaks_certificate": any(
            row["breaks_positive_certificate"] for row in canary_rows
        ),
        "sign_canary_breaks_certificate": bool(sign_canary.upper() < 0),
    }
    all_checks_pass = all(checks.values())
    return {
        "schema": "claude_d3_kappa2_average_pilot.v1",
        "classification": (
            "xi_d3_kappa2_average_pilot_certified"
            if all_checks_pass
            else "xi_d3_kappa2_average_pilot_failed"
        ),
        "parameters": {
            "target_n": target_n,
            "theta_values": theta_values,
            "arb_precision_bits": args.prec,
            "finite_source": str(finite_source),
            "canary_factors": factor_texts,
            "integration_module": "tools/claude_kappa3_window_certificate.py",
        },
        "identity": {
            "Lambda": "Lambda(theta)=log int_R exp(theta*s)V(s) ds",
            "kappa2": "kappa2=Lambda''",
            "average": (
                "A_n=(Lambda(2n)-2Lambda(2n-2)+Lambda(2n-4))/4 "
                "=int_[0,1]^2 kappa2(2(n-2)+2(x+y)) dxdy"
            ),
            "q": "q_n=exp(4A_n)*C_(n-1)/C_n, C_n=(2n)(2n-1)",
            "D": "D_n=(1-q_(n+1))^2-q_(n+1)^2(1-q_n)(1-q_(n+2))",
        },
        "proof_contract": {
            "integration": (
                "acb.integral on the conditioned core plus real-ball hulls outside; "
                "elementary positive left/right tail intervals are added"
            ),
            "kernel_truncation": (
                "reuses the m=15 Phi truncation and explicit m-tail enclosure from "
                "claude_kappa3_window_certificate.py"
            ),
            "double_evaluation": (
                "q and D are evaluated both via exp of the kappa2 square average "
                "and via raw moment ratios"
            ),
            "seam_check": (
                "q_127 and q_128 must overlap the independent H1321 DFT/alias balls"
            ),
        },
        "moments": moment_rows,
        "q_rows": q_rows,
        "finite_overlap": overlap_rows,
        "result": {
            "target_n": target_n,
            "d3_from_kappa2_average": arb_record(d3_average),
            "d3_from_moment_ratio": arb_record(d3_direct),
            "d3_cross_evaluation_difference": arb_record(d3_difference),
            "finite_overlap_pass": all(row["overlaps"] for row in overlap_rows),
        },
        "canaries": {
            "moment_mutations": canary_rows,
            "moment_mutation_kills_certificate": any(
                row["breaks_positive_certificate"] for row in canary_rows
            ),
            "sign_canary": arb_record(sign_canary),
            "sign_canary_kills_certificate": bool(sign_canary.upper() < 0),
        },
        "checks": checks,
        "all_checks_pass": all_checks_pass,
        "decision": (
            f"The exact kappa2-average route certifies D_{target_n}>0 and "
            "overlaps H1321 on q_127,q_128. This is a rigorous pilot row, not "
            "a uniform tail proof."
            if all_checks_pass
            else "At least one moment, overlap, D3, or canary check failed."
        ),
        "scope": {
            "proved_if_pass": f"D_{target_n}>0 for the true Xi moment sequence",
            "does_not_prove": [
                f"D_n>0 uniformly for every n>{target_n}",
                "PF-infinity",
                "all-degree Jensen hyperbolicity",
                "the Riemann Hypothesis",
            ],
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Certify a one-row translated-D3 pilot through kappa2 averages."
    )
    parser.add_argument("--out", default=DEFAULT_OUT)
    parser.add_argument("--finite-source", default=DEFAULT_FINITE_SOURCE)
    parser.add_argument("--target-n", type=int, default=127)
    parser.add_argument("--prec", type=int, default=256)
    parser.add_argument(
        "--canary-factors",
        default="10001:10000,1001:1000,101:100",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        report = build_report(args)
    except Exception as exc:
        report = {
            "schema": "claude_d3_kappa2_average_pilot.v1",
            "classification": "xi_d3_kappa2_average_pilot_replay_unavailable",
            "all_checks_pass": False,
            "error": f"{type(exc).__name__}: {exc}",
            "decision": "The rigorous kappa2-average pilot did not complete.",
        }

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if "result" in report:
        write_markdown(report, out_path.with_suffix(".md"))
    print(
        json.dumps(
            {
                "classification": report["classification"],
                "all_checks_pass": report["all_checks_pass"],
                "out": str(out_path),
                "error": report.get("error"),
            },
            indent=2,
        )
    )
    return 0 if report["all_checks_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
