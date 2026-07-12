"""Rigorous finite translated-D3 certificate for the Xi coefficients.

This replays the H1321 Cauchy/DFT coefficient extraction, enlarges every
coefficient ball by the proved outer-circle alias bound, and only then forms

    q_n = a_n a_(n-2) / a_(n-1)^2,
    D_n = (1-q_(n+1))^2 - q_(n+1)^2 (1-q_n)(1-q_(n+2)).

For each 2 <= n <= 126, D_n is also evaluated as the normalized contiguous
Toeplitz determinant

    det [[a_n,   a_(n+1), a_(n+2)],
         [a_(n-1), a_n,   a_(n+1)],
         [a_(n-2), a_(n-1), a_n  ]] / a_n^3.

All arithmetic used for the certificate is python-flint ball arithmetic.
The coefficient-mutation and sign canaries must fail before the report can
receive a certified classification.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any


THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from rh_h218_jensen_arb_discriminants import (  # noqa: E402
    gamma_coefficients,
    load_flint,
    xi_function,
)


DEFAULT_OUT = "research/riemann/claude_d3_alias_repaired_finite_certificate.json"
DEFAULT_FLINT_PATH = "C:/tmp/h784_pydeps_copy"


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


def ratios(values: list[Any]) -> list[Any | None]:
    out: list[Any | None] = [None]
    out.extend(values[n] / values[n - 1] for n in range(1, len(values)))
    return out


def q_ratios(values: list[Any]) -> list[Any | None]:
    r_values = ratios(values)
    out: list[Any | None] = [None, None]
    out.extend(r_values[n] / r_values[n - 1] for n in range(2, len(values)))
    return out


def translated_d3(q_values: list[Any | None], n: int) -> Any:
    one = q_values[n].parent()(1)
    q_left = q_values[n]
    q_center = q_values[n + 1]
    q_right = q_values[n + 2]
    return (one - q_center) ** 2 - q_center**2 * (one - q_left) * (one - q_right)


def determinant3(matrix: list[list[Any]]) -> Any:
    a, b, c = matrix[0]
    d, e, f = matrix[1]
    g, h, i = matrix[2]
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def normalized_contiguous_d3(values: list[Any], n: int) -> Any:
    matrix = [
        [values[n], values[n + 1], values[n + 2]],
        [values[n - 1], values[n], values[n + 1]],
        [values[n - 2], values[n - 1], values[n]],
    ]
    return determinant3(matrix) / values[n] ** 3


def parse_factor(text: str) -> tuple[int, int]:
    pieces = text.split(":")
    if len(pieces) != 2:
        raise ValueError(f"canary factor must be numerator:denominator, got {text!r}")
    numerator, denominator = int(pieces[0]), int(pieces[1])
    if numerator <= 0 or denominator <= 0:
        raise ValueError("canary factors must be positive")
    return numerator, denominator


def mutation_canaries(
    arb: Any,
    repaired_a: list[Any],
    target_n: int,
    factor_texts: list[str],
) -> list[dict[str, Any]]:
    rows = []
    mutation_index = target_n + 1
    for factor_text in factor_texts:
        numerator, denominator = parse_factor(factor_text)
        factor = arb(numerator) / arb(denominator)
        mutated = list(repaired_a)
        mutated[mutation_index] *= factor
        q_values = q_ratios(mutated)
        value = translated_d3(q_values, target_n)
        direct = normalized_contiguous_d3(mutated, target_n)
        rows.append(
            {
                "target_n": target_n,
                "mutated_a_index": mutation_index,
                "factor": factor_text,
                "d3_from_q": arb_record(value),
                "d3_from_determinant": arb_record(direct),
                "breaks_positive_certificate": bool(value.upper() < 0 and direct.upper() < 0),
            }
        )
    return rows


def write_markdown(report: dict[str, Any], out_path: Path) -> None:
    summary = report.get("finite_d3", {})
    checks = report.get("checks", {})
    lines = [
        "# Alias-Repaired Finite Translated-D3 Certificate",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Result",
        "",
        report["decision"],
        "",
        "The certified inequality is",
        "",
        "```text",
        "D_n = (1-q_(n+1))^2 - q_(n+1)^2 (1-q_n)(1-q_(n+2)) > 0,",
        "q_n = a_n a_(n-2) / a_(n-1)^2.",
        "```",
        "",
        f"- certified integer range: `{summary.get('n_range')}`",
        f"- row count: `{summary.get('row_count')}`",
        f"- worst lower bound: `{summary.get('worst_d3_from_q', {}).get('interval', {}).get('lower')}` at `n={summary.get('worst_d3_from_q', {}).get('n')}`",
        f"- determinant cross-evaluation worst lower bound: `{summary.get('worst_d3_from_determinant', {}).get('interval', {}).get('lower')}`",
        "",
        "## Proof Contract",
        "",
        "The H1321 root-of-unity coefficient balls are recomputed. Each is enlarged by",
        "",
        "```text",
        "Xi(iR) R^(-2n) (r/R)^N / (1-(r/R)^N),  r=11, R=12, N=8192,",
        "```",
        "",
        "before any quotient or determinant is formed. Every row is then evaluated both",
        "from the quotient formula and from the contiguous 3x3 Toeplitz determinant",
        "divided by `a_n^3`.",
        "",
        "## Checks",
        "",
    ]
    for name, value in checks.items():
        lines.append(f"- `{name}`: `{value}`")
    lines.extend(
        [
            "",
            "## Canaries",
            "",
            f"- coefficient mutation killed the target row: `{report.get('canaries', {}).get('coefficient_mutation_kills_certificate')}`",
            f"- sign/comparator mutation killed the target row: `{report.get('canaries', {}).get('sign_canary_kills_certificate')}`",
            "",
            "## Scope",
            "",
        ]
    )
    for item in report["scope"]["does_not_prove"]:
        lines.append(f"- Does not prove {item}.")
    lines.append("")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    if args.n_max < 128:
        raise ValueError("--n-max must be at least 128 to certify D_n through n=126")
    if args.samples <= 2 * args.n_max:
        raise ValueError("--samples must exceed 2*n_max")

    flint, acb, arb, ctx = load_flint(args.flint_path)
    ctx.dps = args.dps
    radius = arb(args.radius)
    outer_radius = arb(args.outer_radius)
    if not radius < outer_radius:
        raise ValueError("--outer-radius must exceed --radius")

    raw_gammas = gamma_coefficients(acb, arb, args.n_max, args.radius, args.samples)
    raw_a = [
        raw_gammas[n] / arb(str(math.factorial(n)))
        for n in range(args.n_max + 1)
    ]

    outer_s = acb(arb("0.5") + outer_radius, arb(0))
    outer_xi = xi_function(acb, arb, outer_s)
    outer_abs_upper = outer_xi.abs_upper()
    rho = (radius / outer_radius) ** args.samples
    alias_geometric = rho / (1 - rho)
    alias_bounds = [
        outer_abs_upper * outer_radius ** (-2 * n) * alias_geometric
        for n in range(args.n_max + 1)
    ]
    repaired_a = [
        raw_a[n] + arb(0, alias_bounds[n].upper())
        for n in range(args.n_max + 1)
    ]

    r_values = ratios(repaired_a)
    q_values = q_ratios(repaired_a)
    runtime_rows: list[dict[str, Any]] = []
    for n in range(2, 127):
        d3_q = translated_d3(q_values, n)
        d3_det = normalized_contiguous_d3(repaired_a, n)
        difference = d3_q - d3_det
        runtime_rows.append(
            {
                "n": n,
                "q_n": q_values[n],
                "q_center": q_values[n + 1],
                "q_right": q_values[n + 2],
                "d3_q": d3_q,
                "d3_det": d3_det,
                "difference": difference,
            }
        )

    worst_q = min(runtime_rows, key=lambda row: row["d3_q"].lower())
    worst_det = min(runtime_rows, key=lambda row: row["d3_det"].lower())
    rows = [
        {
            "n": row["n"],
            "q_n": arb_record(row["q_n"]),
            "q_center": arb_record(row["q_center"]),
            "q_right": arb_record(row["q_right"]),
            "d3_from_q": arb_record(row["d3_q"]),
            "d3_from_determinant": arb_record(row["d3_det"]),
            "cross_evaluation_difference": arb_record(row["difference"]),
        }
        for row in runtime_rows
    ]

    canary_factor_texts = [
        text.strip() for text in args.canary_factors.split(",") if text.strip()
    ]
    canary_rows = mutation_canaries(
        arb,
        repaired_a,
        int(worst_q["n"]),
        canary_factor_texts,
    )
    sign_canary = -worst_q["d3_q"]

    pf2_margins = [r_values[n] - r_values[n + 1] for n in range(1, 128)]
    q_monotone_margins = [q_values[n + 1] - q_values[n] for n in range(2, 128)]
    checks = {
        "python_flint_loaded": True,
        "reference_radius_is_11": args.radius == "11",
        "outer_radius_is_12": args.outer_radius == "12",
        "sample_count_is_8192": args.samples == 8192,
        "n_max_reaches_128": args.n_max >= 128,
        "outer_xi_real_part_positive": bool(outer_xi.real.lower() > 0),
        "outer_xi_imaginary_part_contains_zero": bool(
            outer_xi.imag.lower() <= 0 <= outer_xi.imag.upper()
        ),
        "alias_ratio_strictly_between_zero_and_one": bool(0 < rho < 1),
        "all_alias_bounds_positive": all(bound.lower() > 0 for bound in alias_bounds),
        "all_repaired_coefficients_positive": all(value.lower() > 0 for value in repaired_a),
        "pf2_ratios_strictly_decreasing_through_R128": all(
            margin.lower() > 0 for margin in pf2_margins
        ),
        "q2_strictly_below_one_half": bool(q_values[2].upper() < arb(1) / 2),
        "q_strictly_increasing_through_q128": all(
            margin.lower() > 0 for margin in q_monotone_margins
        ),
        "row_indices_exactly_2_through_126": [row["n"] for row in runtime_rows]
        == list(range(2, 127)),
        "all_d3_q_lower_bounds_positive": all(
            row["d3_q"].lower() > 0 for row in runtime_rows
        ),
        "all_d3_determinant_lower_bounds_positive": all(
            row["d3_det"].lower() > 0 for row in runtime_rows
        ),
        "all_cross_evaluation_differences_contain_zero": all(
            row["difference"].lower() <= 0 <= row["difference"].upper()
            for row in runtime_rows
        ),
        "coefficient_mutation_canary_breaks_certificate": any(
            row["breaks_positive_certificate"] for row in canary_rows
        ),
        "sign_canary_breaks_certificate": bool(sign_canary.upper() < 0),
    }
    all_checks_pass = all(checks.values())

    report = {
        "schema": "claude_d3_alias_repaired_finite_certificate.v1",
        "classification": (
            "xi_translated_d3_alias_repaired_finite_certified"
            if all_checks_pass
            else "xi_translated_d3_alias_repaired_finite_failed"
        ),
        "parameters": {
            "n_max": args.n_max,
            "samples": args.samples,
            "reference_radius": args.radius,
            "outer_radius": args.outer_radius,
            "dps": args.dps,
            "flint_path": args.flint_path,
            "python_flint_version": getattr(flint, "__version__", "unknown"),
            "canary_factors": canary_factor_texts,
        },
        "definitions": {
            "Xi": "Xi(z)=xi(1/2+i*z)=sum_n (-1)^n a_n z^(2n)",
            "q_n": "q_n=a_n*a_(n-2)/a_(n-1)^2",
            "e_n": "e_n=1-q_n",
            "D_n": "D_n=e_(n+1)^2-q_(n+1)^2*e_n*e_(n+2)",
            "determinant_identity": (
                "D_n=det(T_n)/a_n^3 with T_n rows "
                "[a_n,a_(n+1),a_(n+2)], [a_(n-1),a_n,a_(n+1)], "
                "[a_(n-2),a_(n-1),a_n]"
            ),
        },
        "proof_contract": {
            "kernel_majorant": "sup_|z|=R |Xi(z)| <= Xi(iR)=xi(1/2+R)",
            "dft_alias_identity": (
                "the N-point DFT coefficient at exponent 2n aliases only "
                "exponents 2n+lN"
            ),
            "cauchy_alias_bound": (
                "|a_tilde_n-a_n| <= M_R R^(-2n) "
                "(r/R)^N/(1-(r/R)^N)"
            ),
            "interval_policy": (
                "every raw coefficient ball is enlarged symmetrically by the "
                "analytic alias bound before quotients and determinants"
            ),
            "double_evaluation": (
                "D_n is evaluated through q-ratios and through a direct 3x3 "
                "Toeplitz determinant"
            ),
        },
        "outer_majorant": {
            "xi_ball": {
                "real": arb_record(outer_xi.real),
                "imag": arb_record(outer_xi.imag),
            },
            "rho": arb_record(rho),
            "alias_geometric_factor": arb_record(alias_geometric),
        },
        "selected_alias_bounds": [
            {"n": n, "raw_a": arb_record(raw_a[n]), "alias_bound": arb_record(alias_bounds[n])}
            for n in (0, 2, 64, 126, 127, 128)
        ],
        "finite_d3": {
            "n_range": [2, 126],
            "row_count": len(rows),
            "worst_d3_from_q": {
                "n": worst_q["n"],
                "interval": arb_record(worst_q["d3_q"]),
            },
            "worst_d3_from_determinant": {
                "n": worst_det["n"],
                "interval": arb_record(worst_det["d3_det"]),
            },
            "rows": rows,
        },
        "canaries": {
            "target_n": worst_q["n"],
            "coefficient_mutations": canary_rows,
            "coefficient_mutation_kills_certificate": any(
                row["breaks_positive_certificate"] for row in canary_rows
            ),
            "sign_canary": arb_record(sign_canary),
            "sign_canary_kills_certificate": bool(sign_canary.upper() < 0),
        },
        "checks": checks,
        "all_checks_pass": all_checks_pass,
        "decision": (
            "The alias-repaired Xi coefficient balls rigorously certify translated "
            "D3 for every integer 2<=n<=126. This closes the finite seam only; "
            "an analytic or independently certified tail for n>=127 is still required."
            if all_checks_pass
            else "At least one finite-D3, alias, cross-evaluation, or canary check failed."
        ),
        "scope": {
            "proved_if_pass": "D_n>0 for the true Xi coefficients for all 2<=n<=126",
            "does_not_prove": [
                "translated D3 for n>=127",
                "full PF3 (only H811 consecutive-row consequences are relevant)",
                "PF-infinity",
                "all-degree Jensen hyperbolicity",
                "the Riemann Hypothesis",
            ],
        },
    }
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Certify the finite translated-D3 seam using H1321 alias-repaired Xi balls."
    )
    parser.add_argument("--out", default=DEFAULT_OUT)
    parser.add_argument("--flint-path", default=DEFAULT_FLINT_PATH)
    parser.add_argument("--dps", type=int, default=280)
    parser.add_argument("--n-max", type=int, default=128)
    parser.add_argument("--samples", type=int, default=8192)
    parser.add_argument("--radius", default="11")
    parser.add_argument("--outer-radius", default="12")
    parser.add_argument(
        "--canary-factors",
        default="10001:10000,1001:1000,101:100",
        help="comma-separated positive rational factors numerator:denominator",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        report = build_report(args)
    except Exception as exc:
        report = {
            "schema": "claude_d3_alias_repaired_finite_certificate.v1",
            "classification": "xi_translated_d3_alias_repaired_finite_replay_unavailable",
            "all_checks_pass": False,
            "error": f"{type(exc).__name__}: {exc}",
            "decision": "The rigorous finite-D3 replay did not complete.",
            "scope": {"does_not_prove": ["any new translated-D3 statement"]},
        }

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
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
